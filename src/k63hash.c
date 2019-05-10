#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "attribute.h"
#include "io_utils.h"
#include "k63hash.h"
#include "utils.h"
#include "verbose.h"

#define KMHASH_IDLE			0
#define KMHASH_BUSY			1

#define __atomic_k63_equal(x, y)					       \
	(*(volatile uint64_t *)(&((x).bin[0])) == *(volatile uint64_t *)(&((y).bin[0])) && \
	 *(volatile uint64_t *)(&((x).bin[1])) == *(volatile uint64_t *)(&((y).bin[1])))

void save_k63hash(struct k63hash_t *h, const char *path)
{
	FILE *fp = xfopen(path, "wb");
	char *sig = "kh63";
	xfwrite(sig, 4, 1, fp);
	xfwrite(&h->size, sizeof(kmint_t), 1, fp);
	xfwrite(&h->n_item, sizeof(kmint_t), 1, fp);
	xfwrite(h->keys, sizeof(k63key_t), h->size, fp);
	xfwrite(h->adjs, sizeof(uint8_t), h->size, fp);
	xfwrite(h->flag, sizeof(uint8_t), h->size, fp);
	fclose(fp);
}

int load_k63hash(struct k63hash_t *h, const char *path)
{
	FILE *fp = xfopen(path, "rb");
	if (!fp)
		return 0;
	char sig[5];
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, "kh63", 4))
		return 0;
	xfread(&h->size, sizeof(kmint_t), 1, fp);
	xfread(&h->n_item, sizeof(kmint_t), 1, fp);
	h->n_probe = estimate_probe_3(h->size);
	h->keys = malloc(h->size * sizeof(k63key_t));
	h->adjs = malloc(h->size * sizeof(uint8_t));
	h->flag = malloc(h->size * sizeof(uint8_t));
	xfread(h->keys, sizeof(k63key_t), h->size, fp);
	xfread(h->adjs, sizeof(uint8_t), h->size, fp);
	xfread(h->flag, sizeof(uint8_t), h->size, fp);
	fclose(fp);
	return 1;
}

struct k63resize_bundle_t {
	struct k63hash_t *h;
	int n_threads;
	int thread_no;
	int adj_included;
	pthread_barrier_t *barrier;
};

static kmint_t internal_k63hash_put(struct k63hash_t *h, k63key_t key)
{
	kmint_t mask, step, i, n_probe;
	uint64_t k;
	uint8_t cur_flag;

	mask = h->size - 1;
	k = __hash_k63(key);
	n_probe = h->n_probe;
	i = k & mask;
	step = 0;
	do {
		cur_flag = atomic_val_CAS8(&(h->flag[i]), KMFLAG_EMPTY,
								KMFLAG_LOADING);
		if (cur_flag == KMFLAG_EMPTY) {
			*((volatile k63key_t *)(h->keys + i)) = key;
			atomic_add_and_fetch_kmint(&(h->n_item), 1);
			*((volatile uint8_t *)(h->flag + i)) = KMFLAG_NEW;
			return i;
		} else if (cur_flag == KMFLAG_LOADING) {
			continue;
		} else { /* KMFLAG_NEW */
			if (__atomic_k63_equal(h->keys[i], key))
				return i;
			++step;
			i = (i + step * (step + 1) / 2) & mask;
		}
	} while (step <= n_probe);
	return KMHASH_END(h);
}

static kmint_t internal_k63hash_singleton_put(struct k63hash_t *h, k63key_t key)
{
	kmint_t mask, step, i, n_probe;
	uint64_t k;
	uint8_t cur_flag;

	mask = h->size - 1;
	k = __hash_k63(key);
	n_probe = h->n_probe;
	i = k & mask;
	step = 0;
	do {
		cur_flag = atomic_val_CAS8(&(h->flag[i]), KMFLAG_EMPTY,
								KMFLAG_LOADING);
		if (cur_flag == KMFLAG_EMPTY) {
			*((volatile k63key_t *)(h->keys + i)) = key;
			atomic_add_and_fetch_kmint(&(h->n_item), 1);
			*((volatile uint8_t *)(&(h->flag[i]))) = KMFLAG_NEW;
			return i;
		} else if (cur_flag == KMFLAG_LOADING) {
			continue;
		} else { /* KMFLAG_NEW */
			assert(cur_flag == KMFLAG_NEW);
			if (__atomic_k63_equal(h->keys[i], key)) {
				atomic_set_bit32(h->sgts + (i >> 5), i & 31);
				return i;
			}
			++step;
			i = (i + step * (step + 1) / 2) & mask;
		}
	} while (step <= n_probe);
	return KMHASH_END(h);
}

kmint_t k63hash_get(struct k63hash_t *h, k63key_t key)
{
	kmint_t step, i, n_probe, mask;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_k63(key);
	i = k & mask;
	n_probe = h->n_probe;
	step = 0;
	do {
		if (*((volatile uint8_t *)(h->flag + i)) == KMFLAG_EMPTY)
			return KMHASH_END(h);
		else if (*((volatile uint8_t *)(h->flag + i)) == KMFLAG_LOADING)
			continue;
		if (__atomic_k63_equal(h->keys[i], key))
			return i;
		++step;
		i = (i + step * (step + 1) / 2) & mask;
	} while (step <= n_probe);
	return KMHASH_END(h);
}

void *k63hash_resize_worker(void *data)
{
	struct k63resize_bundle_t *bundle = (struct k63resize_bundle_t *)data;
	struct k63hash_t *h = bundle->h;
	int adj_included = bundle->adj_included;

	kmint_t i, j, step, l, r, cap, size, old_size, mask, n_probe;
	k63key_t x, xt;
	uint64_t k;
	uint32_t b, bt;
	uint8_t  a, at, current_flag;

	size = h->size;
	old_size = h->old_size;
	mask = size - 1;
	n_probe = h->n_probe;

	/* Init all keys */
	cap = (size - old_size) / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), size - old_size);
	for (i = l; i < r; ++i) {
		if (adj_included)
			h->adjs[old_size + i] = 0;
		h->flag[old_size + i] = KMFLAG_EMPTY;
	}

	cap = old_size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), old_size);
	for (i = l; i < r; ++i) {
		if (h->flag[i] == KMFLAG_NEW)
			h->flag[i] = KMFLAG_OLD;
	}

	pthread_barrier_wait(bundle->barrier);

	/* Re-positioning items */
	for (i = l; i < r; ++i) {
		if (atomic_val_CAS8(h->flag + i, KMFLAG_OLD,
				KMFLAG_LOADING) == KMFLAG_OLD) {
			x = h->keys[i];
			b = atomic_get_bit32(h->sgts + (i >> 5), i & 31);
			atomic_clear_bit32(h->sgts + (i >> 5), i & 31);
			if (adj_included) {
				a = h->adjs[i];
				h->adjs[i] = 0;
			}
			*((volatile uint8_t *)(h->flag + i)) = KMFLAG_EMPTY;
			while (1) { /* kick-out process */
				k = __hash_k63(x);
				j = k & mask;
				step = 0;
				current_flag = KMFLAG_NEW;
				while (step <= n_probe) {
					if ((current_flag =
						atomic_val_CAS8(h->flag + j,
							KMFLAG_OLD, KMFLAG_NEW))
								== KMFLAG_OLD) {
						xt = h->keys[j];
						h->keys[j] = x;
						x = xt;
						bt = atomic_get_bit32(h->sgts +
							(j >> 5), j & 31);
						atomic_set_bit_var32(h->sgts +
							(j >> 5), j & 31, b);
						b = bt;
						if (adj_included) {
							at = h->adjs[j];
							h->adjs[j] = a;
							a = at;
						}
						break;
					} else if ((current_flag =
						atomic_val_CAS8(h->flag + j,
							KMFLAG_EMPTY, KMFLAG_NEW))
								== KMFLAG_EMPTY) {
						h->keys[j] = x;
						atomic_set_bit_val32(h->sgts +
							(j >> 5), j & 31, b);
						if (adj_included)
							h->adjs[j] = a;
						break;
					}  else if (current_flag == KMFLAG_LOADING) {
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

void k63hash_resize_multi(struct k63hash_t *h, int adj_included)
{
	int n_threads, i;
	n_threads = h->n_worker;

	h->old_size = h->size;
	h->size <<= 1;
	h->n_probe = estimate_probe_3(h->size);

	h->keys = realloc(h->keys, h->size * sizeof(k63key_t));
	h->sgts = realloc(h->sgts, (h->size >> 5) * sizeof(uint32_t));
	if (adj_included)
		h->adjs = realloc(h->adjs, h->size * sizeof(uint8_t));
	h->flag = realloc(h->flag, h->size * sizeof(uint8_t));
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

	struct k63resize_bundle_t *bundles;
	bundles = calloc(n_threads, sizeof(struct k63resize_bundle_t));

	for (i = 0; i < n_threads; ++i) {
		bundles[i].n_threads = n_threads;
		bundles[i].thread_no = i;
		bundles[i].h = h;
		bundles[i].barrier = &barrier;
		bundles[i].adj_included = adj_included;
		pthread_create(t + i, &attr, k63hash_resize_worker, bundles + i);
	}

	for (i = 0; i < n_threads; ++i)
		pthread_join(t[i], NULL);

	pthread_attr_destroy(&attr);
	pthread_barrier_destroy(&barrier);

	free(t);
	free(bundles);
}

void k63hash_resize_single(struct k63hash_t *h, int adj_included)
{
	k63key_t *keys;
	uint32_t *sgts;
	uint8_t  *adjs, *flag;

	kmint_t i, j, mask, n_probe, step, size, old_size;
	k63key_t x, xt;
	uint64_t k;
	uint32_t b, bt;
	uint8_t current_flag, a, at;

	old_size = h->size;
	h->size <<= 1;
	size = h->size;
	mask = size - 1;
	h->n_probe = estimate_probe_3(h->size);
	n_probe = h->n_probe;

	h->keys = realloc(h->keys, h->size * sizeof(k63key_t));
	h->sgts = realloc(h->sgts, (h->size >> 5) * sizeof(uint32_t));
	memset(h->sgts + (old_size >> 5), 0,
			((size >> 5) - (old_size >> 5)) * sizeof(uint32_t));
	if (adj_included) {
		h->adjs = realloc(h->adjs, h->size * sizeof(uint8_t));
		adjs = h->adjs;
		memset(adjs + old_size, 0, (size - old_size) * sizeof(uint8_t));
	}
	h->flag = realloc(h->flag, h->size * sizeof(uint8_t));
	memset(h->flag + old_size, KMFLAG_EMPTY,
		(size - old_size) * sizeof(uint8_t));

	keys = h->keys;
	sgts = h->sgts;
	flag = h->flag;

	for (i = 0; i < old_size; ++i)
		if (flag[i] == KMFLAG_NEW)
			flag[i] = KMFLAG_OLD;

	for (i = 0; i < old_size; ++i) {
		if (flag[i] == KMFLAG_OLD) {
			x = keys[i];
			b = get_bit32(sgts[i >> 5], i & 31);
			sgts[i >> 5] = clear_bit32(sgts[i >> 5], i & 31);
			if (adj_included) {
				a = adjs[i];
				adjs[i] = 0;
			}
			flag[i] = KMFLAG_EMPTY;
			while (1) { /* kick-out process */
				k = __hash_k63(x);
				j = k & mask;
				step = 0;
				current_flag = KMFLAG_NEW;
				while (step <= n_probe) {
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
					j = (j + step * (step + 1) / 2) & mask;
				}
				if (current_flag == KMFLAG_EMPTY)
					break;
				else if (current_flag == KMFLAG_NEW)
					__ERROR("[Single-thread] Resizing kmhash failed");
			}
		}
	}
}

void k63hash_resize(struct k63hash_t *h, int adj_included)
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
		k63hash_resize_single(h, adj_included);
	else
		k63hash_resize_multi(h, adj_included);

	for (i = 0; i < h->n_worker; ++i)
		pthread_mutex_unlock(h->locks + i);
}

void k63hash_add_edge(struct k63hash_t *h, k63key_t key, int c, pthread_mutex_t *lock)
{
	kmint_t k;
	pthread_mutex_lock(lock);
	k = k63hash_get(h, key);
	assert(k != KMHASH_END(h));
	atomic_set_bit8(h->adjs + k, c);
	pthread_mutex_unlock(lock);
}

void k63hash_put_adj(struct k63hash_t *h, k63key_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	pthread_mutex_lock(lock);
	k = internal_k63hash_singleton_put(h, key);
	pthread_mutex_unlock(lock);

	while (k == KMHASH_END(h)) {
		/* ensure only one thread resize the hash table */
		if (atomic_bool_CAS32(&(h->status),
				KMHASH_IDLE, KMHASH_BUSY)) {
			k63hash_resize(h, 1);
			atomic_val_CAS32(&(h->status),
				KMHASH_BUSY, KMHASH_IDLE);
		}

		pthread_mutex_lock(lock);
		k = internal_k63hash_singleton_put(h, key);
		pthread_mutex_unlock(lock);
	}
}

void k63hash_put(struct k63hash_t *h, k63key_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	pthread_mutex_lock(lock);
	k = internal_k63hash_singleton_put(h, key);
	pthread_mutex_unlock(lock);

	while (k == KMHASH_END(h)) {
		if (atomic_bool_CAS32(&(h->status),
				KMHASH_IDLE, KMHASH_BUSY)) {
			k63hash_resize(h, 0);
			atomic_val_CAS32(&(h->status),
				KMHASH_BUSY, KMHASH_IDLE);
		}

		pthread_mutex_lock(lock);
		k = internal_k63hash_singleton_put(h, key);
		pthread_mutex_unlock(lock);
	}
}

void k63hash_filter(struct k63hash_t *h, int adj_included)
{
	kmint_t n_item, i, j, new_size, cur_size, n_probe, mask, step;
	k63key_t x, xt;
	uint64_t k;
	uint8_t a, at, current_flag;

	cur_size = h->size;
	n_item = 0;
	for (i = 0; i < cur_size; ++i) {
		if ((h->sgts[i >> 5] >> (i & 31)) & 1)
			++n_item;
	}

	new_size = n_item;
	__round_up_kmint(new_size);
	/* if cur_size can hold, why not? */
	if (new_size > cur_size)
		new_size = cur_size;
	n_probe = estimate_probe_3(new_size);
	mask = new_size - 1;

	/* filtering singleton kmer */
	for (i = 0; i < cur_size; ++i) {
		assert(h->flag[i] == KMFLAG_NEW || h->flag[i] == KMFLAG_EMPTY);
		if (h->flag[i] == KMFLAG_NEW &&
			((h->sgts[i >> 5] >> (i & 31)) & 1))
			h->flag[i] = KMFLAG_OLD;
		else
			h->flag[i] = KMFLAG_EMPTY;
	}

	int retry = 0;
loop_refill:
	if (retry) {
		for (i = 0; i < cur_size; ++i) {
			if (h->flag[i] == KMFLAG_NEW)
				h->flag[i] = KMFLAG_OLD;
		}
		for (i = 0; i < cur_size; ++i) {
			if (h->flag[i] == KMFLAG_EMPTY) {
				h->keys[i] = x;
				if (adj_included)
					h->adjs[i] = a;
				h->flag[i] = KMFLAG_OLD;
				break;
			}
		}
		new_size <<= 1;
		n_probe = estimate_probe_3(new_size);
		mask = new_size - 1;
		retry = 0;
	}
	for (i = 0; i < cur_size; ++i) {
		if (h->flag[i] == KMFLAG_OLD) {
			x = h->keys[i];
			if (adj_included) {
				a = h->adjs[i];
				h->adjs[i] = 0;
			}
			h->flag[i] = KMFLAG_EMPTY;
			while (1) {
				k = __hash_k63(x);
				j = k & mask;
				step = 0;

				current_flag = KMFLAG_NEW;
				while (step <= n_probe) {
					current_flag = h->flag[j];
					if (current_flag == KMFLAG_EMPTY) {
						h->flag[j] = KMFLAG_NEW;
						h->keys[j] = x;
						if (adj_included)
							h->adjs[j] = a;
						break;
					} else if (current_flag == KMFLAG_OLD) {
						h->flag[j] = KMFLAG_NEW;
						xt = h->keys[j];
						h->keys[j] = x;
						x = xt;
						if (adj_included) {
							at = h->adjs[j];
							h->adjs[j] = a;
							a = at;
						}
						break;
					}
					++step;
					j = (j + step * (step + 1) / 2) & mask;
				}
				if (current_flag == KMFLAG_EMPTY) {
					break;
				} else if (current_flag == KMFLAG_NEW) {
					if (new_size == cur_size)
						__ERROR("Shrink kmhash error");
					retry = 1;
					goto loop_refill;
				}
			}
		}
	}
	h->size = new_size;
	h->n_probe = n_probe;
	h->n_item = n_item;
	if (new_size < cur_size) {
		h->keys = realloc(h->keys, new_size * sizeof(k63key_t));
		h->flag = realloc(h->flag, new_size * sizeof(uint8_t));
		if (adj_included)
			h->adjs = realloc(h->adjs, new_size * sizeof(uint8_t));
	}
	free(h->sgts);
	h->sgts = NULL;
}

void k63hash_test(struct k63hash_t *h)
{
	kmint_t i, k;
	for (i = 0; i < h->size; ++i) {
		if (h->flag[i] == KMFLAG_EMPTY)
			continue;
		k = k63hash_get(h, h->keys[i]);
		if (k != i) {
			__ERROR("Hash fail at postion: get_hash_pos(h[%lu] = (%lu, %lu)) = %lu\n",
				i, h->keys[i].bin[0], h->keys[i].bin[1], k);
		}
	}
}

void k63hash_init(struct k63hash_t *h, kmint_t size, int n_threads, int adj_included)
{
	int k;

	h->size = size;
	__round_up_kmint(h->size);
	h->n_probe = estimate_probe_3(h->size);
	h->keys = malloc(h->size * sizeof(k63key_t));
	h->sgts = calloc(h->size >> 5, sizeof(uint32_t));
	h->flag = calloc(h->size, sizeof(uint8_t));
	h->n_item = 0;
	if (adj_included)
		h->adjs = calloc(h->size, sizeof(uint8_t));
	else
		h->adjs = NULL;

	h->n_worker = n_threads;
	h->status = KMHASH_IDLE;
	h->locks = malloc(n_threads * sizeof(pthread_mutex_t));
	for (k = 0; k < n_threads; ++k)
		pthread_mutex_init(h->locks + k, NULL);
}

void k63hash_destroy(struct k63hash_t *h)
{
	if (!h) return;
	int i;
	for (i = 0; i < h->n_worker; ++i)
		pthread_mutex_destroy(h->locks + i);
	free(h->locks);
	free(h->keys);
	free(h->sgts);
	free(h->adjs);
	free(h->flag);
	h->locks = NULL;
	h->keys = NULL;
	h->sgts = NULL;
	h->adjs = NULL;
	h->flag = NULL;
}

/****************** Single threaded k63 hash table ****************************/

void k63_convert_table(struct k63hash_t *h, struct k63_idhash_t *dict)
{
	dict->size = h->size;
	dict->n_probe = h->n_probe;
	dict->n_item = h->n_item;

	dict->keys = h->keys;
	dict->adjs = h->adjs;
	dict->flag = h->flag;
	dict->id   = calloc(dict->size, sizeof(gint_t));

	h->keys = NULL;
	h->adjs = NULL;
	h->flag = NULL;
}

void k63_idhash_init(struct k63_idhash_t *h, kmint_t size, int adj_included)
{
	h->size = size;
	__round_up_kmint(h->size);
	h->n_probe = estimate_probe_3(h->size);
	h->n_item = 0;
	h->keys = malloc(h->size * sizeof(k63key_t));
	h->id = calloc(h->size, sizeof(gint_t));
	h->flag = calloc(h->size, sizeof(uint8_t));
	if (adj_included)
		h->adjs = calloc(h->size, sizeof(uint8_t));
	else
		h->adjs = NULL;
}

void k63_idhash_clean(struct k63_idhash_t *h)
{
	free(h->keys);
	free(h->adjs);
	free(h->id);
	free(h->flag);
	h->keys = NULL;
	h->id = NULL;
	h->adjs = NULL;
	h->flag = NULL;
	h->n_item = h->size = h->n_probe = 0;
}

kmint_t k63_idhash_get(struct k63_idhash_t *h, k63key_t key)
{
	kmint_t step, i, n_probe, mask;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_k63(key);
	i = k & mask;
	n_probe = h->n_probe;
	step = 0;
	do {
		if (h->flag[i] == KMFLAG_EMPTY)
			return IDHASH_END(h);
		if (__k63_equal(h->keys[i], key))
			return i;
		++step;
		i = (i + step * (step + 1) / 2) & mask;
	} while (step <= n_probe);
	return IDHASH_END(h);
}

static inline kmint_t internal_k63_idhash_put(struct k63_idhash_t *h, k63key_t key)
{
	kmint_t mask, step, i, n_probe;
	uint64_t k;

	mask = h->size - 1;
	k = __hash_k63(key);
	n_probe = h->n_probe;
	i = k & mask;
	step = 0;
	do {
		if (h->flag[i] == KMFLAG_EMPTY || (h->flag[i] == KMFLAG_NEW &&
						   __k63_equal(h->keys[i], key))) {
			if (h->flag[i] == KMFLAG_EMPTY) {
				h->keys[i] = key;
				h->flag[i] = KMFLAG_NEW;
				++h->n_item;
			}
			return i;
		}
		++step;
		i = (i + step * (step + 1) / 2) & mask;
	} while (step <= n_probe);
	return IDHASH_END(h);
}

static void k63_idhash_resize(struct k63_idhash_t *h)
{
	kmint_t i, old_size, n_probe, mask;
	old_size = h->size;
	h->size <<= 1;
	h->n_probe = estimate_probe_3(h->size);
	h->keys = realloc(h->keys, h->size * sizeof(k63key_t));
	h->id = realloc(h->id, h->size * sizeof(gint_t));
	h->flag = realloc(h->flag, h->size * sizeof(uint8_t));

	mask = h->size - 1;
	n_probe = h->n_probe;

	for (i = old_size; i < h->size; ++i) {
		h->flag[i] = KMFLAG_EMPTY;
		h->id[i] = -1;
	}

	for (i = 0; i < old_size; ++i) {
		if (h->flag[i] == KMFLAG_NEW)
			h->flag[i] = KMFLAG_OLD;
	}

	for (i = 0; i < old_size; ++i) {
		if (h->flag[i] == KMFLAG_OLD) {
			k63key_t x = h->keys[i], xt;
			gint_t   y = h->id[i],   yt;
			h->flag[i] = KMFLAG_EMPTY;
			h->id[i] = -1;
			while (1) {
				uint64_t k = __hash_k63(x);
				kmint_t j = k & mask;
				kmint_t step = 0;
				uint8_t current_flag = KMFLAG_NEW;
				while (step <= n_probe) {
					current_flag = h->flag[j];
					if (current_flag == KMFLAG_EMPTY) {
						h->flag[j] = KMFLAG_NEW;
						h->keys[j] = x;
						h->id[j] = y;
						break;
					} else if (current_flag == KMFLAG_OLD) {
						h->flag[j] = KMFLAG_NEW;
						xt = h->keys[j];
						yt = h->id[j];
						h->keys[j] = x;
						h->id[j] = y;
						x = xt;
						y = yt;
						break;
					}
					++step;
					j = (j + step * (step + 1) / 2) & mask;
				}
				if (current_flag == KMFLAG_EMPTY)
					break;
				else if (current_flag == KMFLAG_NEW)
					__ERROR("Resizing node table error");
			}
		}
	}
}

kmint_t k63_idhash_put(struct k63_idhash_t *h, k63key_t key)
{
	kmint_t k;
	k = internal_k63_idhash_put(h, key);
	while (k == IDHASH_END(h)) {
		k63_idhash_resize(h);
		k = internal_k63_idhash_put(h, key);
	}
	return k;
}

