#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "attribute.h"
#include "io_utils.h"
#include "k31hash.h"
#include "utils.h"
#include "verbose.h"

#define rs_set_old(x, i) ((x)[(i) >> 4] = ((x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1)))) |       \
				((uint32_t)KMFLAG_OLD << (((i) & 15) << 1)))

#define rs_set_new(x, i) ((x)[(i) >> 4] = ((x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1)))) |       \
				((uint32_t)KMFLAG_NEW << (((i) & 15) << 1)))

#define rs_set_empty(x, i) ((x)[(i) >> 4] = (x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1))))

#define rs_get_flag(x, i) (((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)

#define rs_is_old(x, i) ((((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)  \
							== (uint32_t)KMFLAG_OLD)

#define rs_is_empty(x, i) ((((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)  \
							== (uint32_t)KMFLAG_EMPTY)

#define rs_is_new(x, i) ((((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)  \
							== (uint32_t)KMFLAG_NEW)

void save_k31hash(struct k31hash_t *h, const char *path)
{
	FILE *fp = xfopen(path, "wb");
	char *sig = "kh31";
	xfwrite(sig, 4, 1, fp);
	xfwrite(&h->size, sizeof(kmint_t), 1, fp);
	xfwrite(&h->n_item, sizeof(kmint_t), 1, fp);
	xfwrite(h->keys, sizeof(k31key_t), h->size, fp);
	xfwrite(h->adjs, sizeof(uint8_t), h->size, fp);
	fclose(fp);
}

int load_k31hash(struct k31hash_t *h, const char *path)
{
	FILE *fp = xfopen(path, "rb");
	if (!fp)
		return 0;
	char sig[5];
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, "kh31", 4))
		return 0;
	xfread(&h->size, sizeof(kmint_t), 1, fp);
	xfread(&h->n_item, sizeof(kmint_t), 1, fp);
	h->n_probe = estimate_probe_3(h->size);
	h->keys = malloc(h->size * sizeof(k31key_t));
	h->adjs = malloc(h->size * sizeof(uint8_t));
	xfread(h->keys, sizeof(k31key_t), h->size, fp);
	xfread(h->adjs, sizeof(uint8_t), h->size, fp);
	fclose(fp);
	return 1;
}

struct k31resize_bundle_t {
	struct k31hash_t *h;
	int n_threads;
	int thread_no;
	int adj_included;
	pthread_barrier_t *barrier;
};

static kmint_t internal_k31hash_put(struct k31hash_t *h, k31key_t key)
{
	kmint_t mask, step, i, n_probe;
	k31key_t cur_key;
	uint64_t k;

	mask = h->size - 1;
	k = __hash_k31(key);
	n_probe = h->n_probe;
	i = k & mask;
	step = 0;
	do {
		i = (i + step * (step + 1) / 2) & mask;
		cur_key = atomic_val_CAS_k31key(h->keys + i, K31_NULL, key);
		++step;
	} while (step <= n_probe && cur_key != key && cur_key != K31_NULL);
	if (cur_key == K31_NULL || cur_key == key) {
		if (cur_key == K31_NULL)
			atomic_add_and_fetch_kmint(&(h->n_item), 1);
		return i;
	}
	return KMHASH_END(h);
}

static kmint_t internal_k31hash_singleton_put(struct k31hash_t *h, k31key_t key)
{
	kmint_t mask, step, i, n_probe;
	k31key_t cur_key;
	uint64_t k;

	mask = h->size - 1;
	k = __hash_k31(key);
	n_probe = h->n_probe;
	i = k & mask;
	step = 0;
	do {
		i = (i + step * (step + 1) / 2) & mask;
		cur_key = atomic_val_CAS_k31key(h->keys + i, K31_NULL, key);
		++step;
	} while (step <= n_probe && cur_key != key && cur_key != K31_NULL);
	if (cur_key == K31_NULL || cur_key == key) {
		if (cur_key == K31_NULL)
			atomic_add_and_fetch_kmint(&(h->n_item), 1);
		else
			atomic_set_bit32(h->sgts + (i >> 5), i & 31);
		return i;
	}
	return KMHASH_END(h);
}

kmint_t k31hash_get(struct k31hash_t *h, k31key_t key)
{
	kmint_t step, i, n_probe, mask;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_k31(key);
	i = k & mask;
	n_probe = h->n_probe;
	step = 0;
	do {
		i = (i + step * (step + 1) / 2) & mask;
		if (h->keys[i] == key)
			return i;
		++step;
	} while (step <= n_probe && h->keys[i] != K31_NULL);
	return KMHASH_END(h);
}

kmint_t k31hash_get_deb(struct k31hash_t *h, k31key_t key)
{
	kmint_t step, i, n_probe, mask;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_k31(key);
	i = k & mask;
	n_probe = h->n_probe;
	step = 0;
	do {
		i = (i + step * (step + 1) / 2) & mask;
		fprintf(stderr, "i = %lu; keys[i] = %lu; key = %lu\n", i, h->keys[i], key);
		if (h->keys[i] == key)
			return i;
		++step;
	} while (step <= n_probe && h->keys[i] != K31_NULL);
	return KMHASH_END(h);
}

void *k31hash_resize_worker(void *data)
{
	struct k31resize_bundle_t *bundle = (struct k31resize_bundle_t *)data;
	struct k31hash_t *h;

	k31key_t *keys;
	uint32_t *sgts;
	uint8_t  *adjs, *flag;

	kmint_t i, j, step, l, r, cap, size, old_size, mask, n_probe;
	k31key_t x, xt;
	uint64_t k;
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
		keys[old_size + i] = K31_NULL;
		if (adj_included)
			adjs[old_size + i] = 0;
	}

	cap = old_size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), old_size);
	for (i = l; i < r; ++i) {
		if (keys[i] != K31_NULL)
			flag[i] = KMFLAG_OLD;
	}

	pthread_barrier_wait(bundle->barrier);

	/* Re-positioning items */
	for (i = l; i < r; ++i) {
		if (atomic_val_CAS8(flag + i, KMFLAG_OLD, KMFLAG_LOADING)
								== KMFLAG_OLD) {
			x = keys[i];
			keys[i] = K31_NULL;
			b = atomic_get_bit32(sgts + (i >> 5), i & 31);
			atomic_clear_bit32(sgts + (i >> 5), i & 31);
			if (adj_included) {
				a = adjs[i];
				adjs[i] = 0;
			}
			flag[i] = KMFLAG_EMPTY;
			while (1) { /* kick-out process */
				k = __hash_k31(x);
				j = k & mask;
				step = 0;
				current_flag = KMFLAG_NEW;
				while (step <= n_probe) {
					if ((current_flag =
						atomic_val_CAS8(flag + j,
							KMFLAG_OLD, KMFLAG_NEW))
								== KMFLAG_OLD) {
						xt = keys[j];
						keys[j] = x;
						x = xt;
						bt = atomic_get_bit32(sgts +
							(j >> 5), j & 31);
						atomic_set_bit_var32(sgts +
							(j >> 5), j & 31, b);
						b = bt;
						if (adj_included) {
							at = adjs[j];
							adjs[j] = a;
							a = at;
						}
						break;
					} else if ((current_flag =
						atomic_val_CAS8(flag + j,
							KMFLAG_EMPTY, KMFLAG_NEW))
							== KMFLAG_EMPTY) {
						keys[j] = x;
						atomic_set_bit_val32(sgts +
							(j >> 5), j & 31, b);
						if (adj_included)
							adjs[j] = a;
						break;
					} else if (current_flag == KMFLAG_LOADING) {
						/* not sure */
						continue;
					}
					++step;
					j = (j + step * (step + 1) / 2) & mask;
				}
				if (current_flag == KMFLAG_EMPTY)
					break;
				else if (current_flag == KMFLAG_NEW)
					__ERROR("[Multi-thread] Resize kmhash failed");
			}
		}
	}
	pthread_exit(NULL);
}

void k31hash_resize_multi(struct k31hash_t *h, int adj_included)
{
	int n_threads, i;
	n_threads = h->n_worker;

	h->old_size = h->size;
	h->size <<= 1;
	h->n_probe = estimate_probe_3(h->size);

	h->keys = realloc(h->keys, h->size * sizeof(k31key_t));
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

	struct k31resize_bundle_t *bundles;
	bundles = calloc(n_threads, sizeof(struct k31resize_bundle_t));

	for (i = 0; i < n_threads; ++i) {
		bundles[i].n_threads = n_threads;
		bundles[i].thread_no = i;
		bundles[i].h = h;
		bundles[i].barrier = &barrier;
		bundles[i].adj_included = adj_included;
		pthread_create(t + i, &attr, k31hash_resize_worker, bundles + i);
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

void k31hash_resize_single(struct k31hash_t *h, int adj_included)
{
	k31key_t *keys;
	uint32_t *sgts;
	uint8_t  *adjs, *flag;

	kmint_t i, j, mask, n_probe, step, size, old_size;
	k31key_t x, xt;
	uint64_t k;
	uint32_t b, bt;
	uint8_t current_flag, a, at;

	old_size = h->size;
	h->size <<= 1;
	size = h->size;
	mask = size - 1;
	h->n_probe = estimate_probe_3(h->size);
	n_probe = h->n_probe;

	h->keys = realloc(h->keys, h->size * sizeof(k31key_t));
	h->sgts = realloc(h->sgts, (h->size >> 5) * sizeof(uint32_t));
	memset(h->sgts + (old_size >> 5), 0,
			((size >> 5) - (old_size >> 5)) * sizeof(uint32_t));
	if (adj_included) {
		h->adjs = realloc(h->adjs, h->size * sizeof(uint8_t));
		adjs = h->adjs;
		memset(adjs + old_size, 0, (size - old_size) * sizeof(uint8_t));
	}
	flag = calloc(h->size, sizeof(uint8_t));
	keys = h->keys;
	sgts = h->sgts;

	for (i = old_size; i < size; ++i)
		keys[i] = K31_NULL;

	for (i = 0; i < old_size; ++i)
		if (keys[i] != K31_NULL)
			flag[i] = KMFLAG_OLD;

	for (i = 0; i < old_size; ++i) {
		if (flag[i] == KMFLAG_OLD) {
			x = keys[i];
			keys[i] = K31_NULL;
			b = get_bit32(sgts[i >> 5], i & 31);
			sgts[i >> 5] = clear_bit32(sgts[i >> 5], i & 31);
			if (adj_included) {
				a = adjs[i];
				adjs[i] = 0;
			}
			flag[i] = KMFLAG_EMPTY;
			while (1) { /* kick-out process */
				k = __hash_k31(x);
				j = k & mask;
				step = 0;
				current_flag = KMFLAG_NEW;
				while (step <= n_probe) {
					j = (j + step * (step + 1) / 2) & mask;
					current_flag = flag[j];
					if (current_flag == KMFLAG_EMPTY) {
						flag[j] = KMFLAG_NEW;
						keys[j] = x;
						sgts[j >> 5] = set_bit_val32(sgts[j >> 5],
									j & 31, b);
						if (adj_included)
							adjs[j] = a;
						break;
					} else if (current_flag == KMFLAG_OLD) {
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
						flag[j] = KMFLAG_NEW;
						break;
					}
					++step;
				}
				if (current_flag == KMFLAG_EMPTY)
					break;
				else if (current_flag == KMFLAG_NEW)
					__ERROR("[Single-thread] Resizing kmhash failed");
			}
		}
	}
	free(flag);
}

void k31hash_resize(struct k31hash_t *h, int adj_included)
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
		k31hash_resize_single(h, adj_included);
	else
		k31hash_resize_multi(h, adj_included);

	for (i = 0; i < h->n_worker; ++i)
		pthread_mutex_unlock(h->locks + i);
}

void k31hash_add_edge(struct k31hash_t *h, k31key_t key, int c, pthread_mutex_t *lock)
{
	kmint_t k;
	pthread_mutex_lock(lock);
	k = k31hash_get(h, key);
	assert(k != KMHASH_END(h));
	atomic_set_bit8(h->adjs + k, c);
	pthread_mutex_unlock(lock);
}

void k31hash_put_adj(struct k31hash_t *h, k31key_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	pthread_mutex_lock(lock);
	k = internal_k31hash_singleton_put(h, key);
	pthread_mutex_unlock(lock);

	while (k == KMHASH_END(h)) {
		/* ensure only one thread resize the hash table */
		if (atomic_bool_CAS32(&(h->status),
				KMHASH_IDLE, KMHASH_BUSY)) {
			k31hash_resize(h, 1);
			atomic_val_CAS32(&(h->status),
				KMHASH_BUSY, KMHASH_IDLE);
		}

		pthread_mutex_lock(lock);
		k = internal_k31hash_singleton_put(h, key);
		pthread_mutex_unlock(lock);
	}
}

void k31hash_put(struct k31hash_t *h, k31key_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	pthread_mutex_lock(lock);
	k = internal_k31hash_singleton_put(h, key);
	pthread_mutex_unlock(lock);

	while (k == KMHASH_END(h)) {
		if (atomic_bool_CAS32(&(h->status),
				KMHASH_IDLE, KMHASH_BUSY)) {
			k31hash_resize(h, 0);
			atomic_val_CAS32(&(h->status),
				KMHASH_BUSY, KMHASH_IDLE);
		}

		pthread_mutex_lock(lock);
		k = internal_k31hash_singleton_put(h, key);
		pthread_mutex_unlock(lock);
	}
}

void k31hash_filter(struct k31hash_t *h, int adj_included)
{
	k31key_t *keys = h->keys;
	uint32_t *sgts = h->sgts;
	uint8_t *adjs = h->adjs;

	kmint_t n_item, i, j, new_size, cur_size, n_probe, mask, step;
	k31key_t x, xt;
	uint64_t k;
	uint8_t a, at;

	cur_size = h->size;
	n_item = 0;
	for (i = 0; i < cur_size; ++i) {
		if ((sgts[i >> 5] >> (i & 31)) & 1)
			++n_item;
		else
			keys[i] = K31_NULL;
	}

	new_size = n_item;
	__round_up_kmint(new_size);
	/* if cur_size can hold, why not? */
	if (new_size > cur_size)
		new_size = cur_size;
	n_probe = estimate_probe_3(new_size);
	mask = new_size - 1;

	/* filtering singleton kmer */
	uint32_t *flag = calloc(cur_size >> 4, sizeof(uint32_t));

	for (i = 0; i < cur_size; ++i)
		if ((sgts[i >> 5] >> (i & 31)) & 1)
			rs_set_old(flag, i);

	int retry = 0;
loop_refill:
	if (retry) {
		fprintf(stderr, "Retry shrink #%d\n", retry);
		for (i = 0; i < cur_size; ++i) {
			if (rs_is_new(flag, i))
				rs_set_old(flag, i);
		}
		for (i = 0; i < cur_size; ++i) {
			if (rs_is_empty(flag, i)) {
				keys[i] = x;
				if (adj_included)
					adjs[i] = a;
				rs_set_old(flag, i);
				break;
			}
		}
		new_size <<= 1;
		n_probe = estimate_probe_3(new_size);
		mask = new_size - 1;
	}
	for (i = 0; i < cur_size; ++i) {
		if (rs_is_old(flag, i)) {
			assert(keys[i] != K31_NULL);
			x = keys[i];
			keys[i] = K31_NULL;
			if (adj_included) {
				a = adjs[i];
				adjs[i] = 0;
			}
			rs_set_empty(flag, i);
			while (1) {
				k = __hash_k31(x);
				j = k & mask;
				step = 0;
				uint32_t current_flag = KMFLAG_NEW;
				while (step <= n_probe) {
					j = (j + step * (step + 1) / 2) & mask;
					current_flag = rs_get_flag(flag, j);
					if (current_flag == KMFLAG_EMPTY) {
						rs_set_new(flag, j);
						keys[j] = x;
						if (adj_included)
							adjs[j] = a;
						break;
					} else if (current_flag == KMFLAG_OLD) {
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
				if (current_flag == KMFLAG_EMPTY)
					break;
				else if (current_flag == KMFLAG_NEW) {
					if (cur_size == new_size)
						__ERROR("Shrink kmhash error");
					++retry;
					goto loop_refill;
				}
			}
		}
	}
	h->size = new_size;
	h->n_probe = n_probe;
	h->n_item = n_item;
	if (new_size < cur_size) {
		h->keys = realloc(h->keys, new_size * sizeof(k31key_t));
		if (adj_included)
			h->adjs = realloc(h->adjs, new_size * sizeof(uint8_t));
	}
	free(h->sgts);
	free(flag);
	h->sgts = NULL;
}

void k31hash_init(struct k31hash_t *h, kmint_t size, int n_threads, int adj_included)
{
	kmint_t i;
	int k;

	h->size = size;
	__round_up_kmint(h->size);
	h->n_probe = estimate_probe_3(h->size);
	h->old_size = 0;
	h->n_item = 0;
	h->keys = malloc(h->size * sizeof(k31key_t));
	h->sgts = calloc(h->size >> 5, sizeof(uint32_t));
	if (adj_included)
		h->adjs = calloc(h->size, sizeof(uint8_t));
	else
		h->adjs = NULL;
	for (i = 0; i < h->size; ++i)
		h->keys[i] = K31_NULL;

	h->n_worker = n_threads;
	h->status = KMHASH_IDLE;
	h->locks = malloc(n_threads * sizeof(pthread_mutex_t));
	for (k = 0; k < n_threads; ++k)
		pthread_mutex_init(h->locks + k, NULL);
}

void k31hash_destroy(struct k31hash_t *h)
{
	if (!h) return;
	int i;
	for (i = 0; i < h->n_worker; ++i)
		pthread_mutex_destroy(h->locks + i);
	free(h->locks);
	free(h->keys);
	free(h->sgts);
	free(h->adjs);
	h->locks = NULL;
	h->keys = NULL;
	h->sgts = NULL;
	h->adjs = NULL;
}

/****************** Single threaded k31 hash table ****************************/

void k31_convert_table(struct k31hash_t *h, struct k31_idhash_t *dict)
{
	dict->size = h->size;
	dict->n_probe = h->n_probe;
	dict->n_item = h->n_item;

	dict->keys = h->keys;
	dict->adjs = h->adjs;
	dict->id   = calloc(dict->size, sizeof(gint_t));

	h->keys = NULL;
	h->adjs = NULL;
}

void k31_idhash_init(struct k31_idhash_t *h, kmint_t size, int adj_included)
{
	h->size = size;
	__round_up_kmint(h->size);
	h->n_probe = estimate_probe_3(h->size);
	h->n_item = 0;
	h->keys = malloc(h->size * sizeof(k31key_t));
	h->id = calloc(h->size, sizeof(gint_t));
	if (adj_included)
		h->adjs = calloc(h->size, sizeof(uint8_t));
	kmint_t i;
	for (i = 0; i < h->size; ++i)
		h->keys[i] = K31_NULL;
}

void k31_idhash_clean(struct k31_idhash_t *h)
{
	free(h->keys);
	free(h->adjs);
	free(h->id);
	h->keys = NULL;
	h->id = NULL;
	h->adjs = NULL;
	h->n_item = h->size = h->n_probe = 0;
}

kmint_t k31_idhash_get(struct k31_idhash_t *h, k31key_t key)
{
	kmint_t step, i, n_probe, mask;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_k31(key);
	i = k & mask;
	n_probe = h->n_probe;
	step = 0;
	do {
		i = (i + step * (step + 1) / 2) & mask;
		if (h->keys[i] == key)
			return i;
		++step;
	} while (step <= n_probe && h->keys[i] != K31_NULL);
	return IDHASH_END(h);
}

static inline kmint_t internal_k31_idhash_put(struct k31_idhash_t *h, k31key_t key)
{
	kmint_t mask, step, i, n_probe;
	uint64_t k;

	mask = h->size - 1;
	k = __hash_k31(key);
	n_probe = h->n_probe;
	i = k & mask;

	step = 0;
	do {
		i = (i + step * (step + 1) / 2) & mask;
		if (h->keys[i] == K31_NULL || h->keys[i] == key) {
			if (h->keys[i] == K31_NULL) {
				h->keys[i] = key;
				++h->n_item;
			}
			return i;
		}
		++step;
	} while (step <= n_probe);
	return IDHASH_END(h);
}

static void k31_idhash_resize(struct k31_idhash_t *h)
{
	kmint_t i, old_size, n_probe, mask;
	old_size = h->size;
	h->size <<= 1;
	h->n_probe = estimate_probe_3(h->size);
	h->keys = realloc(h->keys, h->size * sizeof(k31key_t));
	h->id = realloc(h->id, h->size * sizeof(gint_t));
	uint32_t *flag = calloc(h->size >> 4, sizeof(uint32_t));

	mask = h->size - 1;
	n_probe = h->n_probe;

	for (i = old_size; i < h->size; ++i) {
		h->keys[i] = K31_NULL;
		h->id[i] = -1;
	}

	for (i = 0; i < old_size; ++i) {
		if (h->keys[i] != K31_NULL)
			rs_set_old(flag, i);
	}

	for (i = 0; i < old_size; ++i) {
		if (rs_is_old(flag, i)) {
			k31key_t x = h->keys[i], xt;
			gint_t   y = h->id[i],   yt;
			rs_set_empty(flag, i);
			h->keys[i] = K31_NULL;
			h->id[i] = -1;
			while (1) {
				uint64_t k = __hash_k31(x);
				kmint_t j = k & mask;
				kmint_t step = 0;
				uint32_t current_flag = KMFLAG_NEW;

				while (step <= n_probe) {
					j = (j + step * (step + 1) / 2) & mask;
					current_flag = rs_get_flag(flag, j);
					if (current_flag == KMFLAG_EMPTY) {
						rs_set_new(flag, j);
						h->keys[j] = x;
						h->id[j] = y;
						break;
					} else if (current_flag == KMFLAG_OLD) {
						rs_set_new(flag, j);
						xt = h->keys[j];
						yt = h->id[j];
						h->keys[j] = x;
						h->id[j] = y;
						x = xt;
						y = yt;
						break;
					}
					++step;
				}
				if (current_flag == KMFLAG_EMPTY)
					break;
				else if (current_flag == KMFLAG_NEW)
					__ERROR("Resizing node table error");
			}
		}
	}
	free(flag);
}

kmint_t k31_idhash_put(struct k31_idhash_t *h, k31key_t key)
{
	kmint_t k;
	k = internal_k31_idhash_put(h, key);
	while (k == IDHASH_END(h)) {
		k31_idhash_resize(h);
		k = internal_k31_idhash_put(h, key);
	}
	return k;
}

/***************************** Barcode hash ***********************************/

void barcode_hash_init(struct barcode_hash_t *h, uint32_t size)
{
	h->size = size - 1;
	__round_up_32(h->size);
	h->n_item = 0;
	h->keys = malloc(h->size * sizeof(uint64_t));
	h->cnts = calloc(h->size, sizeof(uint32_t));
	uint32_t i;
	for (i = 0; i < h->size; ++i)
		h->keys[i] = (uint64_t)-1;
}

void barcode_hash_clean(struct barcode_hash_t *h)
{
	free(h->keys);
	free(h->cnts);
	h->keys = NULL;
	h->cnts = NULL;
	h->n_item = h->size = 0;
}

uint32_t barcode_hash_get(struct barcode_hash_t *h, uint64_t key)
{
	uint32_t i, last, mask, step = 0;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_int(key); i = k & mask;
	last = i;
	while (h->keys[i] != K31_NULL && h->keys[i] != key) {
		i = (i + (++step)) & mask;
		if (i == last) return BARCODE_HASH_END(h);
	}
	return h->keys[i] == key ? i : BARCODE_HASH_END(h);
}

static inline uint32_t internal_barcode_hash_put(struct barcode_hash_t *h, uint64_t key)
{
	if (h->n_item >= (uint32_t)(h->size * BARCODE_HASH_UPPER))
		return BARCODE_HASH_END(h);
	uint32_t i, last, mask, step = 0;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_int(key); i = k & mask;
	if (h->keys[i] == (uint64_t)-1) {
		h->keys[i] = key;
		++h->n_item;
		return i;
	}
	last = i;
	while (h->keys[i] != (uint64_t)-1 && h->keys[i] != key) {
		i = (i + (++step)) & mask;
		if (i == last)
			break;
	}
	if (h->keys[i] == (uint64_t)-1 || h->keys[i] == key) {
		if (h->keys[i] == (uint64_t)-1) {
			h->keys[i] = key;
			++h->n_item;
		}
		return i;
	} else {
		return BARCODE_HASH_END(h);
	}
}

static void barcode_hash_resize(struct barcode_hash_t *h)
{
	uint32_t old_size, mask, i;
	old_size = h->size;
	h->size <<= 1;
	mask = h->size - 1;
	h->keys = realloc(h->keys, h->size * sizeof(uint64_t));
	h->cnts = realloc(h->cnts, h->size * sizeof(uint32_t));
	uint8_t *flag = calloc(h->size, sizeof(uint8_t));

	for (i = old_size; i < h->size; ++i) {
		h->keys[i] = K31_NULL;
		h->cnts[i] = 0;
	}

	for (i = 0; i < old_size; ++i) {
		if (h->keys[i] != K31_NULL)
			flag[i] = KMFLAG_OLD;
	}

	for (i = 0; i < old_size; ++i) {
		if (flag[i] == KMFLAG_OLD) {
			uint64_t x = h->keys[i], xt;
			uint32_t y = h->cnts[i], yt;
			flag[i] = KMFLAG_EMPTY;
			h->keys[i] = K31_NULL;
			h->cnts[i] = 0;
			while (1) {
				uint64_t k = __hash_k31(x);
				uint32_t j = k & mask, step = 0, last;
				last = j;
				while (flag[j] != KMFLAG_EMPTY && flag[j] != KMFLAG_OLD) {
					j = (j + (++step)) & mask;
					if (j == last)
						break;
				}
				if (flag[j] == KMFLAG_EMPTY) {
					flag[j] = KMFLAG_NEW;
					h->keys[j] = x;
					h->cnts[j] = y;
					break;
				} else if (flag[j] == KMFLAG_OLD) {
					flag[j] = KMFLAG_NEW;
					xt = h->keys[j];
					yt = h->cnts[j];
					h->keys[j] = x;
					h->cnts[j] = y;
					x = xt; y = yt;
				} else {
					__ERROR("Resize barcode hash failed");
				}
			}
		}
	}
	free(flag);
}

uint32_t barcode_hash_put(struct barcode_hash_t *h, uint64_t key)
{
	uint32_t k;
	// pthread_mutex_lock(h->lock);
	k = internal_barcode_hash_put(h, key);
	while (k == BARCODE_HASH_END(h)) {
		barcode_hash_resize(h);
		k = internal_barcode_hash_put(h, key);
	}
	// pthread_mutex_unlock(h->lock);
	return k;
}

uint32_t barcode_hash_inc_count(struct barcode_hash_t *h, uint64_t key)
{
	uint32_t k;
	// pthread_mutex_lock(h->lock);
	k = internal_barcode_hash_put(h, key);
	while (k == BARCODE_HASH_END(h)) {
		barcode_hash_resize(h);
		k = internal_barcode_hash_put(h, key);
	}
	++h->cnts[k];
	// pthread_mutex_unlock(h->lock);
	return k;
}

void barcode_hash_merge(struct barcode_hash_t *dst, struct barcode_hash_t *src)
{
	uint32_t k, j;
	for (k = 0; k < src->size; ++k) {
		if (src->keys[k] == K31_NULL)
			continue;
		j = barcode_hash_put(dst, src->keys[k]);
		dst->cnts[j] += src->cnts[k];
	}
}

void barcode_hash_clone(struct barcode_hash_t *dst, struct barcode_hash_t *src)
{
	dst->size = src->size;
	dst->n_item = src->n_item;
	dst->keys = malloc(dst->size * sizeof(uint64_t));
	dst->cnts = malloc(dst->size * sizeof(uint32_t));
	memcpy(dst->keys, src->keys, dst->size * sizeof(uint64_t));
	memcpy(dst->cnts, src->cnts, dst->size * sizeof(uint32_t));
}
