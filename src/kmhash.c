#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "kmhash.h"
#include "semaphore_wrapper.h"
#include "utils.h"
#include "verbose.h"

#define KMHASH_IDLE			0
#define KMHASH_BUSY			1

#define KMMASK_EMPTY			0
#define KMMASK_OLD			1
#define KMMASK_NEW			2
#define KMMASK_LOADING			3

#define HM_MAGIC_1			UINT64_C(0xbf58476d1ce4e5b9)
#define HM_MAGIC_2			UINT64_C(0x94d049bb133111eb)

#define atomic_add_and_fetch_kmint	atomic_add_and_fetch64
#define atomic_val_CAS_kmkey		atomic_val_CAS64

#define rs_set_old(x, i) ((x)[(i) >> 4] = ((x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1)))) |       \
				((uint32_t)KMMASK_OLD << (((i) & 15) << 1)))

#define rs_set_new(x, i) ((x)[(i) >> 4] = ((x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1)))) |       \
				((uint32_t)KMMASK_NEW << (((i) & 15) << 1)))

#define rs_set_empty(x, i) ((x)[(i) >> 4] = (x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1))))

#define rs_get_flag(x, i) (((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)

#define rs_is_old(x, i) ((((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)  \
							== (uint32_t)KMMASK_OLD)

struct kmresize_bundle_t {
	struct kmhash_t *h;
	int n_threads;
	int thread_no;
	int adj_included;
	pthread_barrier_t *barrier;
};

static inline kmkey_t __hash_int2(kmkey_t k)
{
	kmkey_t x = k;
	x = (x ^ (x >> 30)) * HM_MAGIC_1;
	x = (x ^ (x >> 27)) * HM_MAGIC_2;
	x ^= (x > 31);
	return x;
}

static inline kmint_t estimate_probe_3(kmint_t size)
{
	kmint_t s, i;
	i = s = 0;
	while (s < size) {
		++i;
		s += i * i * i * 64;
	}
	return i;
}

/*
 * put an item to hash table
 * return position if put ok or already on table
 * otherwise, return KMHASH_MAX_SIZE
 */
static kmint_t internal_kmhash_put(struct kmhash_t *h, kmkey_t key)
{
	kmint_t mask, step, i, n_probe;
	kmkey_t cur_key, k;

	mask = h->size - 1;
	k = __hash_int2(key);
	n_probe = h->n_probe;
	i = k & mask;
	{
		cur_key = atomic_val_CAS_kmkey(&(h->keys[i]), TOMB_STONE, key);
		if (cur_key == TOMB_STONE || cur_key == key) {
			if (cur_key == TOMB_STONE)
				atomic_add_and_fetch_kmint(&(h->n_item), 1);
			return i;
		}
	}
	step = 0;
	do {
		++step;
		i = (i + (step * (step + 1)) / 2) & mask;
		cur_key = atomic_val_CAS_kmkey(&(h->keys[i]), TOMB_STONE, key);
	} while (step < n_probe && cur_key != key && cur_key != TOMB_STONE);
	if (cur_key == TOMB_STONE || cur_key == key) {
		if (cur_key == TOMB_STONE)
			atomic_add_and_fetch_kmint(&(h->n_item), 1);
		return i;
	}
	return KMHASH_MAX_SIZE;
}

/*
 * put an item to hash table,
 * and check if its count greater than 1
 */
static kmint_t internal_kmhash_singleton_put(struct kmhash_t *h, kmkey_t key)
{
	kmint_t mask, step, i, n_probe;
	kmkey_t cur_key, k;

	kmkey_t *keys;
	uint32_t *sgts;

	mask = h->size - 1;
	keys = h->keys;
	sgts = h->sgts;
	n_probe = h->n_probe;

	k = __hash_int2(key);
	i = k & mask;
	{
		cur_key = atomic_val_CAS_kmkey(keys + i, TOMB_STONE, key);
		if (cur_key == TOMB_STONE || cur_key == key) {
			if (cur_key == TOMB_STONE)
				atomic_add_and_fetch_kmint(&(h->n_item), 1);
			else
				atomic_set_bit32(sgts + (i >> 5), i & 31);
			return i;
		}
	}
	step = 0;
	do {
		++step;
		i = (i + (step * (step + 1)) / 2) & mask;
		cur_key = atomic_val_CAS_kmkey(keys + i, TOMB_STONE, key);
	} while (step < n_probe && cur_key != key && cur_key != TOMB_STONE);
	if (cur_key == TOMB_STONE || cur_key == key) {
		if (cur_key == TOMB_STONE)
			atomic_add_and_fetch_kmint(&(h->n_item), 1);
		else
			atomic_set_bit32(sgts + (i >> 5), i & 31);
		return i;
	}
	return KMHASH_MAX_SIZE;
}

kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key)
{
	kmint_t step, i, n_probe, mask;
	kmkey_t k;
	mask = h->size - 1;
	k = __hash_int2(key);
	i = k & mask;
	n_probe = h->n_probe;
	if (h->keys[i] == key)
		return i;
	if (h->keys[i] == TOMB_STONE)
		return KMHASH_MAX_SIZE;
	step = 0;
	do {
		++step;
		i = (i + step * (step + 1) / 2) & mask;
		if (h->keys[i] == key)
			return i;
	} while (step < n_probe && h->keys[i] != TOMB_STONE);
	return KMHASH_MAX_SIZE;
}

void *kmhash_resize_worker(void *data)
{
	struct kmresize_bundle_t *bundle = (struct kmresize_bundle_t *)data;
	struct kmhash_t *h;

	kmkey_t  *keys;
	uint32_t *sgts;
	uint8_t  *adjs, *flag;

	kmint_t i, j, step, l, r, cap, size, old_size, mask, n_probe;
	kmkey_t  x, xt, k;
	uint32_t b, bt;
	uint8_t  a, at, current_flag;
	int adj_included;

	adj_included = bundle->adj_included;
	h = bundle->h;
	size = h->size;
	old_size = h->old_size;
	mask = size - 1;
	n_probe = h->n_probe;

	keys = h->keys;

	sgts = h->sgts;
	adjs = h->adjs;
	flag = h->flag;

	/* Init all keys */
	cap = (size - old_size) / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), size - old_size);
	for (i = l; i < r; ++i) {
		keys[old_size + i] = TOMB_STONE;
		if (adj_included)
			adjs[old_size + i] = 0;
	}

	cap = old_size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), old_size);
	for (i = l; i < r; ++i) {
		if (keys[i] != TOMB_STONE)
			flag[i] = KMMASK_OLD;
	}

	pthread_barrier_wait(bundle->barrier);

	/* Re-positioning items */
	for (i = l; i < r; ++i) {
		if (atomic_val_CAS8(flag + i, KMMASK_OLD, KMMASK_LOADING)
								== KMMASK_OLD) {
			x = keys[i];
			keys[i] = TOMB_STONE;
			b = atomic_get_bit32(sgts + (i >> 5), i & 31);
			atomic_clear_bit32(sgts + (i >> 5), i & 31);
			if (adj_included) {
				a = adjs[i];
				adjs[i] = 0;
			}
			flag[i] = KMMASK_EMPTY;

			while (1) { /* kick-out process */
				k = __hash_int2(x);
				j = k & mask;
				step = 0;
				current_flag = KMMASK_NEW;
				while (step <= n_probe) {
					j = (j + step * (step + 1) / 2) & mask;
					if ((current_flag =
						atomic_val_CAS8(flag + j,
								KMMASK_EMPTY,
								KMMASK_NEW))
							== KMMASK_EMPTY) {
						keys[j] = x;
						atomic_set_bit_val32(sgts +
								(j >> 5),
								j & 31, b);
						if (adj_included)
							adjs[j] = a;
						break;
					} else if ((current_flag =
						atomic_val_CAS8(flag + j,
								KMMASK_OLD,
								KMMASK_NEW))
							== KMMASK_OLD) {
						xt = keys[j];
						keys[j] = x;
						x = xt;
						bt = atomic_get_bit32(sgts +
								(j >> 5),
								j & 31);
						atomic_set_bit_var32(sgts +
								(j >> 5),
								j & 31, b);
						b = bt;
						if (adj_included) {
							at = adjs[j];
							adjs[j] = a;
							a = at;
						}
						break;
					}
					++step;
				}
				if (current_flag == KMMASK_EMPTY)
					break;
				else if (current_flag == KMMASK_NEW)
					__ERROR("[Multi-thread] Resize kmhash failed");
			}
		}
	}
	pthread_exit(NULL);
}

void kmhash_resize_multi(struct kmhash_t *h, int adj_included)
{
	int n_threads, i;
	n_threads = h->n_worker;

	h->old_size = h->size;
	h->size <<= 1;
	h->n_probe = estimate_probe_3(h->size);

	h->keys = realloc(h->keys, h->size * sizeof(kmkey_t));
	h->sgts = realloc(h->sgts, (h->size >> 5) * sizeof(uint32_t));
	if (adj_included)
		h->adjs = realloc(h->adjs, h->size * sizeof(uint8_t));
	h->flag = calloc(h->size, sizeof(uint8_t));
	memset(h->sgts + (h->old_size >> 5), 0,
		((h->size >> 5) - (h->old_size >> 5)) * sizeof(uint32_t));

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_t *t;
	t = calloc(n_threads, sizeof(pthread_t));

	pthread_barrier_t barrier;
	pthread_barrier_init(&barrier, NULL, n_threads);

	struct kmresize_bundle_t *bundles;
	bundles = calloc(n_threads, sizeof(struct kmresize_bundle_t));

	for (i = 0; i < n_threads; ++i) {
		bundles[i].n_threads = n_threads;
		bundles[i].thread_no = i;
		bundles[i].h = h;
		bundles[i].barrier = &barrier;
		bundles[i].adj_included = adj_included;
		pthread_create(t + i, &attr, kmhash_resize_worker, bundles + i);
	}

	for (i = 0; i < n_threads; ++i)
		pthread_join(t[i], NULL);

	pthread_attr_destroy(&attr);
	pthread_barrier_destroy(&barrier);

	free(t);
	free(bundles);
	free(h->flag);
	h->flag = NULL;
}

void kmhash_resize_single(struct kmhash_t *h, int adj_included)
{
	kmkey_t  *keys;
	uint32_t *sgts;
	uint8_t  *adjs, *flag;

	kmint_t i, j, mask, n_probe, step, size, old_size;
	kmkey_t x, xt, k;
	uint32_t b, bt;
	uint8_t current_flag, a, at;

	old_size = h->size;
	h->size <<= 1;
	size = h->size;
	mask = size - 1;
	h->n_probe = estimate_probe_3(h->size);
	n_probe = h->n_probe;

	h->keys = realloc(h->keys, h->size * sizeof(kmkey_t));
	h->sgts = realloc(h->sgts, (h->size >> 5) * sizeof(uint32_t));
	if (adj_included) {
		adjs = h->adjs;
		h->adjs = realloc(h->adjs, h->size * sizeof(uint8_t));
		memset(adjs + old_size, 0, (size - old_size) * sizeof(uint8_t));
	}
	flag = calloc(h->size, sizeof(uint8_t));

	keys = h->keys;
	sgts = h->sgts;

	memset(sgts + (old_size >> 5), 0,
			((size >> 5) - (old_size >> 5)) * sizeof(uint32_t));

	for (i = old_size; i < size; ++i)
		keys[i] = TOMB_STONE;

	for (i = 0; i < old_size; ++i)
		if (keys[i] == TOMB_STONE)
			flag[i] = KMMASK_OLD;

	for (i = 0; i < old_size; ++i) {
		if (flag[i] == KMMASK_OLD) {
			x = keys[i];
			keys[i] = TOMB_STONE;
			b = get_bit32(sgts[i >> 5], i & 31);
			sgts[i >> 5] = clear_bit32(sgts[i >> 5], i & 31);
			if (adj_included) {
				a = adjs[i];
				adjs[i] = 0;
			}
			flag[i] = KMMASK_EMPTY;
			while (1) { /* kick-out process */
				k = __hash_int2(x);
				j = k & mask;
				step = 0;
				current_flag = KMMASK_NEW;
				while (step <= n_probe) {
					j = (j + step * (step + 1) / 2) & mask;
					current_flag = flag[j];
					if (current_flag == KMMASK_EMPTY) {
						flag[j] = KMMASK_NEW;
						keys[j] = x;
						sgts[j >> 5] = set_bit_val32(sgts[j >> 5],
									j & 31, b);
						if (adj_included)
							adjs[j] = a;
						break;
					} else if (current_flag == KMMASK_OLD) {
						xt = keys[j];
						keys[j] = x;
						x = xt;
						bt = get_bit32(sgts[j >> 5], j & 31);
						sgts[j >> 5] = set_bit_var32(sgts[j >> 5],
									j & 31, b);
						b = bt;
						if (adj_included) {
							at = adjs[j];
							adjs[j] = a;
							a = at;
						}
						flag[j] = KMMASK_NEW;
						break;
					}
					++step;
				}
				if (current_flag == KMMASK_EMPTY)
					break;
				else if (current_flag == KMMASK_NEW)
					__ERROR("[Single-thread] Resizing kmhash failed");
			}
		}
	}
	free(flag);
}

void kmhash_resize(struct kmhash_t *h, int adj_included)
{
	/* Resize process including adjs array */
	int i;
	for (i = 0; i < h->n_worker; ++i)
		pthread_mutex_lock(h->locks + i);

	if (h->size == KMHASH_MAX_SIZE)
		__ERROR("Hash table size limit (number of items = %llu, maximum size = %llu)",
			(unsigned long long)h->n_item,
			(unsigned long long)KMHASH_MAX_SIZE);

	if (h->size <= KMHASH_SINGLE_RESIZE)
		kmhash_resize_single(h, adj_included);
	else
		kmhash_resize_multi(h, adj_included);

	for (i = 0; i < h->n_worker; ++i)
		pthread_mutex_unlock(h->locks + i);
}

void kmhash_add_edge(struct kmhash_t *h, kmkey_t key, int c, pthread_mutex_t *lock)
{
	kmint_t k;
	pthread_mutex_lock(lock);
	k = kmhash_get(h, key);
	assert(k != KMHASH_MAX_SIZE);
	atomic_set_bit8(h->adjs + k, c);
	pthread_mutex_unlock(lock);
}

void kmhash_put_adj(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	pthread_mutex_lock(lock);
	k = internal_kmhash_singleton_put(h, key);
	pthread_mutex_unlock(lock);

	while (k == KMHASH_MAX_SIZE) {
		/* ensure only one thread resize the hash table */
		if (atomic_bool_CAS32(&(h->status),
				KMHASH_IDLE, KMHASH_BUSY)) {
			kmhash_resize(h, 1);
			atomic_val_CAS32(&(h->status),
				KMHASH_BUSY, KMHASH_IDLE);
		}

		pthread_mutex_lock(lock);
		k = internal_kmhash_singleton_put(h, key);
		pthread_mutex_unlock(lock);
	}
}

void kmhash_put(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	pthread_mutex_lock(lock);
	k = internal_kmhash_singleton_put(h, key);
	pthread_mutex_unlock(lock);

	while (k == KMHASH_MAX_SIZE) {
		if (atomic_bool_CAS32(&(h->status),
				KMHASH_IDLE, KMHASH_BUSY)) {
			kmhash_resize(h, 0);
			atomic_val_CAS32(&(h->status),
				KMHASH_BUSY, KMHASH_IDLE);
		}

		pthread_mutex_lock(lock);
		k = internal_kmhash_singleton_put(h, key);
		pthread_mutex_unlock(lock);
	}
}

void kmhash_filter(struct kmhash_t *h, int adj_included)
{
	kmkey_t  *keys;
	uint32_t *sgts, *flag;
	uint8_t *adjs;

	kmint_t n_item, i, j, new_size, cur_size, cnt, n_probe, mask, step;
	kmkey_t x, xt, k;
	uint8_t a, at, current_flag;

	keys = h->keys;
	sgts = h->sgts;
	adjs = h->adjs;

	cur_size = h->size;

	n_item = 0;
	for (i = 0; i < cur_size; ++i) {
		if ((sgts[i >> 5] >> (i & 31)) & 1)
			++n_item;
		else
			keys[i] = TOMB_STONE;
	}

	new_size = n_item * 1.25;
	__round_up_kmint(new_size);
	/* if cur_size can hold, why not? */
	if (new_size > cur_size)
		new_size = cur_size;
	n_probe = estimate_probe_3(new_size);
	mask = new_size - 1;

	/* filtering singleton kmer */
	flag = calloc(cur_size >> 4, sizeof(uint32_t));

	for (i = 0; i < cur_size; ++i)
		if ((sgts[i >> 5] >> (i & 31)) & 1)
			rs_set_old(flag, i);

	for (i = 0; i < cur_size; ++i) {
		if (rs_is_old(flag, i)) {
			x = keys[i];
			keys[i] = TOMB_STONE;
			if (adj_included) {
				a = adjs[i];
				adjs[i] = 0;
			}
			rs_set_empty(flag, i);
			while (1) {
				k = __hash_int2(x);
				j = k & mask;
				step = 0;

				current_flag = KMMASK_NEW;
				while (step <= n_probe) {
					j = (j + step * (step + 1) / 2) & mask;
					current_flag = rs_get_flag(flag, j);
					if (current_flag == KMMASK_EMPTY) {
						rs_set_new(flag, j);
						keys[j] = x;
						if (adj_included)
							adjs[j] = a;
						break;
					} else if (current_flag == KMMASK_OLD) {
						rs_set_new(flag, j);
						xt = keys[j];
						keys[j] = x;
						x = xt;
						if (adj_included) {
							at = adjs[j];
							adjs[j] = a;
							a = at;
						}
						break;
					}
					++step;
				}
				if (current_flag == KMMASK_EMPTY)
					break;
				else if (current_flag == KMMASK_NEW)
					__ERROR("Shrink kmhash error");
			}
		}
	}
	h->size = new_size;
	h->n_probe = n_probe;
	h->n_item = n_item;
	if (new_size < cur_size)
		h->keys = realloc(h->keys, new_size * sizeof(kmkey_t));
	free(h->sgts);
	free(flag);
	h->sgts = NULL;
}

void kmhash_init(struct kmhash_t *h, kmint_t size, int n_threads, int adj_included)
{
	kmint_t i;
	int k;

	h->size = size;
	__round_up_kmint(h->size);
	h->keys = malloc(h->size * sizeof(kmkey_t));
	h->sgts = calloc(h->size >> 5, sizeof(uint32_t));
	if (adj_included)
		h->adjs = calloc(h->size, sizeof(uint8_t));
	else
		h->adjs = NULL;
	for (i = 0; i < h->size; ++i)
		h->keys[i] = TOMB_STONE;

	h->n_worker = n_threads;
	h->status = KMHASH_IDLE;
	h->locks = malloc(n_threads * sizeof(pthread_mutex_t));
	for (k = 0; k < n_threads; ++k)
		pthread_mutex_init(h->locks + k, NULL);
}

void kmhash_destroy(struct kmhash_t *h)
{
	if (!h) return;
	int i;
	for (i = 0; i < h->n_worker; ++i)
		pthread_mutex_destroy(h->locks + i);
	free(h->locks);
	free(h->keys);
	free(h->sgts);
	free(h->adjs);
}
