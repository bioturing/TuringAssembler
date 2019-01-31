#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "attribute.h"
#include "k63hash.h"
#include "utils.h"
#include "verbose.h"

#define KMHASH_IDLE			0
#define KMHASH_BUSY			1

#define atomic_add_and_fetch_kmint	atomic_add_and_fetch64

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
		i = (i + step * (step + 1) / 2) & mask;
		cur_flag = atomic_val_CAS8(&(h->flag[i]), KMFLAG_EMPTY,
								KMFLAG_LOADING);
		if (cur_flag == KMFLAG_EMPTY) {
			h->keys[i] = key;
			h->flag[i] = KMFLAG_NEW;
			atomic_add_and_fetch_kmint(&(h->n_item), 1);
			return i;
		} else if (cur_flag == KMFLAG_LOADING) {
			continue;
		} else { /* KMFLAG_NEW */
			if (k63_equal(h->keys[i], key))
				return i;
			++step;
		}
	} while (step <= n_probe);
	return KMHASH_MAX_SIZE;
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
		i = (i + step * (step + 1) / 2) & mask;
		cur_flag = atomic_val_CAS8(&(h->flag[i]), KMFLAG_EMPTY,
								KMFLAG_LOADING);
		if (cur_flag == KMFLAG_EMPTY) {
			h->keys[i] = key;
			h->flag[i] = KMFLAG_NEW;
			atomic_add_and_fetch_kmint(&(h->n_item), 1);
			return i;
		} else if (cur_flag == KMFLAG_LOADING) {
			continue;
		} else { /* KMFLAG_NEW */
			if (k63_equal(h->keys[i], key)) {
				atomic_set_bit32(h->sgts + (i >> 5), i & 31);
				return i;
			}
			++step;
		}
	} while (step <= n_probe);
	return KMHASH_MAX_SIZE;
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
		i = (i + step * (step + 1) / 2) & mask;
		if (h->flag[i] == KMFLAG_EMPTY)
			return KMHASH_MAX_SIZE;
		if (k63_equal(h->keys[i], key))
			return i;
		++step;
	} while (step <= n_probe);
	return KMHASH_MAX_SIZE;
}

void *k63hash_resize_worker(void *data)
{
	struct k63resize_bundle_t *bundle = (struct k63resize_bunlde_t *)data;
	struct k63hash_t *h;

	k63key_t *keys;
	uint32_t *sgts;
	uint8_t  *adjs, *flag;

	kmint_t i, j, step, l, r, cap, size, old_size, mask, n_probe;
	k63key_t x, xt;
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
		if (adj_included)
			adjs[old_size + i] = 0;
	}

	cap = old_size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), old_size);
	for (i = l; i < r; ++i) {
		if (flag[i] == KMFLAG_NEW)
			flag[i] = KMFLAG_OLD;
	}

	pthread_barrier_wait(bundle->barrier);

	/* Re-positioning items */
	for (i = l; i < r; ++i) {
		if (atomic_val_CAS8(flag + i, KMFLAG_OLD, KMFLAG_LOADING)
								== KMFLAG_OLD) {
			x = keys[i];
			b = atomic_get_bit32(sgts + (i >> 5), i & 31);
			atomic_clear_bit32(sgts + (i >> 5), i & 31);
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
					j = (j + step * (step + 1) / 2) & mask;
					if ((current_flag =
						atomic_val_CAS8(flag + j,
								KMFLAG_EMPTY,
								KMFLAG_NEW))
							== KMFLAG_EMPTY) {
						keys[j] = x;
						atomic_set_bit_val32(sgts +
								(j >> 5),
								j & 31, b);
						if (adj_included)
							adjs[j] = a;
						break;
					} else if ((current_flag =
						atomic_val_CAS8(flag + j,
								KMFLAG_OLD,
								KMFLAG_NEW))
							== KMFLAG_OLD) {
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
	memset(h->flag + h->old_size, KMFLAG_EMPTY,
		(h->size - h->old_size) * sizeof(uint8_t));
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
	memset(sgts + (old_size >> 5), 0,
			((size >> 5) - (old_size >> 5)) * sizeof(uint32_t));
	if (adj_included) {
		adjs = h->adjs;
		h->adjs = realloc(h->adjs, h->size * sizeof(uint8_t));
		memset(adjs + old_size, 0, (size - old_size) * sizeof(uint8_t));
	}
	h->flag = realloc(h->flag, h->size * sizeof(uint8_t));
	memset(h->flag + old_size, KMFLAG_EMPTY,
		(h->size - h->old_size) * sizeof(uint8_t));

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
	assert(k != KMHASH_MAX_SIZE);
	atomic_set_bit8(h->adjs + k, c);
	pthread_mutex_unlock(lock);
}

void k63hash_put_adj(struct k63hash_t *h, k63key_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	pthread_mutex_lock(lock);
	k = internal_k63hash_singleton_put(h, key);
	pthread_mutex_unlock(lock);

	while (k == KMHASH_MAX_SIZE) {
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

	while (k == KMHASH_MAX_SIZE) {
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
	k63key_t *keys;
	uint32_t *sgts;
	uint8_t *adjs, *flag;

	kmint_t n_item, i, j, new_size, cur_size, cnt, n_probe, mask, step;
	k63key_t x, xt;
	uint64_t k;
	uint8_t a, at, current_flag;

	keys = h->keys;
	sgts = h->sgts;
	adjs = h->adjs;
	flag = h->flag;

	cur_size = h->size;

	n_item = 0;
	for (i = 0; i < cur_size; ++i) {
		if ((sgts[i >> 5] >> (i & 31)) & 1)
			++n_item;
	}

	new_size = n_item * 1.25;
	__round_up_kmint(new_size);
	/* if cur_size can hold, why not? */
	if (new_size > cur_size)
		new_size = cur_size;
	n_probe = estimate_probe_3(new_size);
	mask = new_size - 1;

	/* filtering singleton kmer */
	// flag = calloc(cur_size >> 4, sizeof(uint32_t));

	// for (i = 0; i < cur_size; ++i)
	// 	if ((sgts[i >> 5] >> (i & 31)) & 1)
	// 		rs_set_old(flag, i);
	for (i = 0; i < cur_size; ++i) {
		if (flag[i] == KMFLAG_NEW &&
			((sgts[i >> 5] >> (i & 31)) & 1))
			flag[i] = KMFLAG_OLD;
		else
			flag[i] = KMFLAG_EMPTY;
	}

	for (i = 0; i < cur_size; ++i) {
		if (flag[i] == KMFLAG_OLD) {
			x = keys[i];
			if (adj_included) {
				a = adjs[i];
				adjs[i] = 0;
			}
			flag[i] = KMFLAG_EMPTY;
			while (1) {
				k = __hash_k63(x);
				j = k & mask;
				step = 0;

				current_flag = KMFLAG_NEW;
				while (step <= n_probe) {
					j = (j + step * (step + 1) / 2) & mask;
					current_flag = flag[j];
					if (current_flag == KMFLAG_EMPTY) {
						flag[j] = KMFLAG_NEW;
						keys[j] = x;
						if (adj_included)
							adjs[j] = a;
						break;
					} else if (current_flag == KMFLAG_OLD) {
						flag[j] = KMFLAG_NEW;
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
				else if (current_flag == KMFLAG_NEW)
					__ERROR("Shrink kmhash error");
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

void k63hash_init(struct k63hash_t *h, kmint_t size, int n_threads, int adj_included)
{
	kmint_t i;
	int k;

	h->size = size;
	__round_up_kmint(h->size);
	h->keys = malloc(h->size * sizeof(k63key_t));
	h->sgts = calloc(h->size >> 5, sizeof(uint32_t));
	h->flag = calloc(h->size, sizeof(uint8_t));
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
	h->locks = NULL;
	h->keys = NULL;
	h->sgts = NULL;
	h->adjs = NULL;
}
