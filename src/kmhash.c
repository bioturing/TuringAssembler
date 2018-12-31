#include <stdlib.h>

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

struct kmresize_bundle_t {
	struct kmhash_t *h;
	int n_threads;
	int thread_no;
	pthread_barrier_t *barrier;
};

static kmint_t estimate_probe_3(kmint_t size)
{
	kmint_t s, i;
	i = s = 0;
	while (s < size) {
		++i;
		s += i * i * i * 64;
	}
	return i;
}

static kmint_t estimate_probe_prime(kmint_t size)
{
	kmint_t s, i;
	i = s = 0;
	while (s < size) {
		++i;
		s += i * i * 256;
	}
	return i;
}

// static kmint_t estimate_probe(kmint_t size)
// {
// 	kmint_t s, i;
// 	i = s = 0;
// 	while (s < size) {
// 		++i;
// 		s += i * i * 2048;
// 	}
// 	return i;
// }

static void atomic_set_bit_kmval(kmval_t *ptr, kmint_t i, kmval_t b)
{
	kmval_t old_bin, new_bin, cur_bin;
	cur_bin = *(volatile kmval_t *)ptr;
	do {
		old_bin = cur_bin;
		new_bin = old_bin | ((kmval_t)b << i);
		cur_bin = __sync_val_compare_and_swap(ptr, old_bin, new_bin);
	} while (cur_bin != old_bin);
}

static kmint_t kmhash_set(struct kmhash_t *h, kmkey_t key)
{
	kmint_t mask, step, i, n_probe;
	kmkey_t cur_key, k;

	mask = h->size - 1;
	k = __hash_int2(key);
	n_probe = h->n_probe;
	i = k & mask;
	{
		cur_key = __sync_val_compare_and_swap(&(h->keys[i]), TOMB_STONE, key);
		if (cur_key == TOMB_STONE || cur_key == key) {
			if (cur_key == TOMB_STONE)
				__sync_fetch_and_add(&(h->n_items), 1);
			else
				atomic_set_bit_kmval(h->vals + (i >> KMVAL_LOG), i & KMVAL_MASK, 1);
			return i;
		}
	}
	step = 0;
	do {
		++step;
		i = (i + (step * (step + 1)) / 2) & mask;
		cur_key = __sync_val_compare_and_swap(&(h->keys[i]), TOMB_STONE, key);
	} while (step < n_probe && cur_key != key && cur_key != TOMB_STONE);
	if (cur_key == TOMB_STONE || cur_key == key) {
		if (cur_key == TOMB_STONE)
			__sync_fetch_and_add(&(h->n_items), 1);
		else
			atomic_set_bit_kmval(h->vals + (i >> KMVAL_LOG), i & KMVAL_MASK, 1);
		return i;
	}
	return KMHASH_MAX_SIZE;
}

// static kmint_t kmhash_put(struct kmhash_t *h, kmkey_t key)
// {
// 	kmint_t mask, step, i, n_probe;
// 	kmkey_t cur_key, k;

// 	mask = h->size - 1;
// 	k = __hash_int2(key);
// 	n_probe = h->n_probe;
// 	i = k & mask;
// 	{
// 		cur_key = __sync_val_compare_and_swap(&(h->keys[i]), TOMB_STONE, key);
// 		if (cur_key == TOMB_STONE || cur_key == key) {
// 			if (cur_key == TOMB_STONE)
// 				__sync_fetch_and_add(&(h->n_items), 1);
// 			return i;
// 		}
// 	}
// 	step = 0;
// 	do {
// 		++step;
// 		i = (i + (step * (step + 1)) / 2) & mask;
// 		cur_key = __sync_val_compare_and_swap(&(h->keys[i]), TOMB_STONE, key);
// 	} while (step < n_probe && cur_key != key && cur_key != TOMB_STONE);
// 	if (cur_key == TOMB_STONE || cur_key == key) {
// 		if (cur_key == TOMB_STONE)
// 			__sync_fetch_and_add(&(h->n_items), 1);
// 		return i;
// 	}
// 	return KMHASH_MAX_SIZE;
// }

void *kmresize_worker(void *data)
{
	struct kmresize_bundle_t *bundle = (struct kmresize_bundle_t *)data;
	struct kmhash_t *h;
	kmint_t i, j, step, l, r, cap, size, old_size, mask, n_probe;
	kmkey_t *keys;
	kmval_t *vals, *old_vals;
	uint8_t *rs_flag;
	kmkey_t x, y, k;
	kmval_t b;
	uint8_t cur_flag;

	h = bundle->h;

	size = h->size;
	old_size = h->old_size;
	mask = size - 1;
	n_probe = h->n_probe;

	keys = h->keys;
	vals = h->vals;
	old_vals = h->old_vals;
	rs_flag = h->rs_flag;

	// Init all keys
	cap = (size - old_size) / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), size - old_size);
	for (i = l; i < r; ++i)
		keys[old_size + i] = TOMB_STONE;

	cap = old_size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), old_size);
	for (i = l; i < r; ++i) {
		if (keys[i] != TOMB_STONE)
			rs_flag[i] = KMMASK_OLD;
	}

	pthread_barrier_wait(bundle->barrier);

	// Fill buckets
	for (i = l; i < r; ++i) {
		if (__sync_val_compare_and_swap(rs_flag + i, KMMASK_OLD, KMMASK_EMPTY) == KMMASK_OLD) {
			x = keys[i];
			b = (old_vals[i >> 5] >> (i & 31)) & 1;
			__sync_val_compare_and_swap(keys + i, x, TOMB_STONE);
			// __sync_val_compare_and_swap(rs_flag + i, KMMASK_LOADING, KMMASK_EMPTY);
			// rs_flag[i] = KMMASK_EMPTY;
			while (1) {
				k = __hash_int2(x);
				j = k & mask;
				step = 0;

				cur_flag = KMMASK_NEW;
				while (step <= n_probe) {
					j = (j + (step * (step + 1)) / 2) & mask;
					if ((cur_flag = __sync_val_compare_and_swap(rs_flag + j, KMMASK_EMPTY, KMMASK_NEW)) == KMMASK_EMPTY) {
					// cur_flag = rs_flag[j];
					// if (cur_flag == KMMASK_EMPTY) {
						// rs_flag[j] = KMMASK_NEW;
						keys[j] = x;
						atomic_set_bit_kmval(vals + (j >> KMVAL_LOG), j & KMVAL_MASK, b);
						// vals[j >> 5] |= (b << (j & 31));
						break;
					// } else if (cur_flag == KMMASK_OLD) {
					} else if ((cur_flag = __sync_val_compare_and_swap(rs_flag + j, KMMASK_OLD, KMMASK_NEW) == KMMASK_OLD)) {
						y = keys[j];
						// rs_flag[j] = KMMASK_NEW;
						keys[j] = x;
						atomic_set_bit_kmval(vals + (j >> KMVAL_LOG), j & KMVAL_MASK, b);
						b = (old_vals[j >> 5] >> (j & 31)) & 1;
						// vals[j >> 5] |= (b << (j & 31));
						x = y;
						break;
					}
					++step;
				}
				if (cur_flag == KMMASK_EMPTY)
					break;
				else if (cur_flag == KMMASK_NEW)
					__ERROR("[Multi thread] Resize kmer count hash table fail");
			}
		}
	}

	pthread_exit(NULL);
}

void kmhash_resize_multi(struct kmhash_t *h)
{
	int n_threads, i;
	n_threads = h->n_workers;

	h->old_size = h->size;
	// h->old_keys = h->keys;
	h->old_vals = h->vals;

	h->size <<= 1;
	h->n_probe = estimate_probe_3(h->size);
	// h->keys = malloc(h->size * sizeof(kmkey_t));
	// h->vals = malloc(h->size * sizeof(kmval_t));
	h->keys = realloc(h->keys, h->size * sizeof(kmkey_t));
	h->vals = calloc(h->size >> 5, sizeof(uint32_t));
	h->rs_flag = calloc(h->size, sizeof(uint8_t));

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_t *t;
	t = calloc(h->n_workers, sizeof(pthread_t));

	pthread_barrier_t barrier;
	pthread_barrier_init(&barrier, NULL, n_threads);

	struct kmresize_bundle_t *bundles;
	bundles = calloc(n_threads, sizeof(struct kmresize_bundle_t));

	for (i = 0; i < n_threads; ++i) {
		bundles[i].n_threads = n_threads;
		bundles[i].thread_no = i;
		bundles[i].h = h;
		bundles[i].barrier = &barrier;
		pthread_create(t + i, &attr, kmresize_worker, bundles + i);
	}

	for (i = 0; i < n_threads; ++i)
		pthread_join(t[i], NULL);

	pthread_attr_destroy(&attr);
	pthread_barrier_destroy(&barrier);

	// kmint_t cnt;
	// cnt = 0;
	// for (i = 0; i < h->size; ++i) {
	// 	if (!((h->rs_flag[i] == KMMASK_NEW) == (h->keys[i] != TOMB_STONE))) {
	// 		fprintf(stderr, "\ni = %llu;\nkey = %llu\n;flag = %u;\n",
	// 			(long long unsigned)i, (long long unsigned)h->keys[i],
	// 			(unsigned)h->rs_flag[i]);
	// 	}
	// 	assert((h->rs_flag[i] == KMMASK_NEW) == (h->keys[i] != TOMB_STONE));
	// 	cnt += (h->rs_flag[i] == KMMASK_NEW);
	// }
	// fprintf(stderr, "[Multi] Number of replaced items: %llu\n", (long long unsigned)cnt);

	free(t);
	free(bundles);
	free(h->rs_flag);
	free(h->old_vals);
}

void kmhash_resize_single(struct kmhash_t *h)
{
	kmint_t i, j, mask, n_probe, step, size, old_size;
	kmkey_t k, x, y;
	kmval_t b;
	uint8_t cur_flag;

	kmval_t *vals, *old_vals;
	kmkey_t *keys;
	uint8_t *rs_flag;

	old_size = h->size;
	old_vals = h->vals;

	h->size <<= 1;
	h->n_probe = estimate_probe_3(h->size);
	h->keys = realloc(h->keys, h->size * sizeof(kmkey_t));
	h->vals = calloc(h->size >> 5, sizeof(kmval_t));
	rs_flag = calloc(h->size, sizeof(uint8_t));
	size = h->size;
	mask = size - 1;
	keys = h->keys;
	vals = h->vals;
	n_probe = h->n_probe;

	for (i = old_size; i < size; ++i)
		keys[i] = TOMB_STONE;

	for (i = 0; i < old_size; ++i) {
		if (keys[i] != TOMB_STONE)
			rs_flag[i] = KMMASK_OLD;
	}

	for (i = 0; i < old_size; ++i) {
		if (rs_flag[i] == KMMASK_OLD) {
			x = keys[i];
			b = (old_vals[i >> 5] >> (i & 31)) & 1;
			rs_flag[i] = KMMASK_EMPTY;
			keys[i] = TOMB_STONE;
			while (1) {
				k = __hash_int2(x);
				j = k & mask;
				step = 0;

				cur_flag = KMMASK_NEW;
				while (step <= n_probe) {
					j = (j + (step * (step + 1)) / 2) & mask;
					cur_flag = rs_flag[j];
					if (cur_flag == KMMASK_EMPTY) {
						rs_flag[j] = KMMASK_NEW;
						keys[j] = x;
						vals[j >> 5] |= (b << (j & 31));
						break;
					} else if (cur_flag == KMMASK_OLD) {
						y = keys[j];
						rs_flag[j] = KMMASK_NEW;
						keys[j] = x;
						vals[j >> 5] |= (b << (j & 31));
						x = y;
						b = (old_vals[j >> 5] >> (i & 31)) & 1;
						break;
					}
					++step;
				}
				if (cur_flag == KMMASK_EMPTY)
					break;
				else if (cur_flag == KMMASK_NEW)
					__ERROR("[Single thread] Resize kmer count hash table fail");
			}
		}
	}
	// kmint_t cnt;
	// cnt = 0;
	// for (i = 0; i < size; ++i)
	// 	cnt += (rs_flag[i] == KMMASK_NEW);
	// fprintf(stderr, "[Single] Number of replaced items: %llu\n", (long long unsigned)cnt);
	free(old_vals);
	free(rs_flag);
}

void kmhash_resize(struct kmhash_t *h)
{
	int i;
	for (i = 0; i < h->n_workers; ++i)
		// sem_wrap_wait(&(h->gsem));
		pthread_mutex_lock(h->locks + i);

	// fprintf(stderr, "[DEBUG] resizing, current number of items = %llu\n", (long long unsigned)h->n_items);

	// kmint_t k, cnt;
	// cnt = 0;
	// for (k = 0; k < h->size; ++k)
	// 	cnt += (h->vals[k >> KMVAL_LOG] >> (k & KMVAL_MASK)) & 1;
	// fprintf(stderr, "\n[DEBUG] Before - Number of non-singleton kmer: %llu\n", cnt);
	if (h->size == KMHASH_MAX_SIZE)
		__ERROR("Unable to expand the hash table (max size = %llu)",
			(unsigned long long)KMHASH_MAX_SIZE);

	if (h->size <= KMHASH_SINGLE_RESIZE)
		kmhash_resize_single(h);
	else
		kmhash_resize_multi(h);

	// cnt = 0;
	// for (k = 0; k < h->size; ++k)
	// 	cnt += (h->vals[k >> KMVAL_LOG] >> (k & KMVAL_MASK)) & 1;
	// fprintf(stderr, "\n[DEBUG] After -  Number of non-singleton kmer: %llu\n", cnt);

	for (i = 0; i < h->n_workers; ++i)
		// sem_wrap_post(&(h->gsem));
		pthread_mutex_unlock(h->locks + i);
}

void kmhash_put(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	pthread_mutex_lock(lock);
	k = kmhash_set(h, key);
	pthread_mutex_unlock(lock);

	if (k == KMHASH_MAX_SIZE) {
		do {
			if (__sync_bool_compare_and_swap(&(h->status), KMHASH_IDLE, KMHASH_BUSY)) {
				kmhash_resize(h);
				__sync_val_compare_and_swap(&(h->status), KMHASH_BUSY, KMHASH_IDLE);
			}

			pthread_mutex_lock(lock);
			k = kmhash_set(h, key);
			pthread_mutex_unlock(lock);
		} while (k == KMHASH_MAX_SIZE);
	}
}

kmint_t kmphash_get(struct kmhash_t *h, kmkey_t key)
{
	kmint_t size, step, i, n_probe;
	kmkey_t k;
	size = h->size;
	k = __hash_int2(key);
	i = k % size;
	// i = key % size;
	n_probe = h->n_probe;
	if (h->keys[i] == key)
		return i;
	step = 0;
	do {
		++step;
		i = (i + step * (step + 1) / 2) % size;
		// if (i >= size)
			// i %= size;
		if (h->keys[i] == key)
			return i;
	} while (step < n_probe);
	return KMHASH_MAX_SIZE;
}

// kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key)
// {
// 	kmint_t mask, step, i, n_probe;
// 	kmkey_t k;
// 	mask = h->size - 1;
// 	k = __hash_int2(key);
// 	i = k & mask;
// 	n_probe = h->n_probe;
// 	if (h->keys[i] == key)
// 		return i;
// 	step = 0;
// 	do {
// 		++step;
// 		i = (i + (step * (step + 1)) / 2) & mask;
// 		if (h->keys[i] == key)
// 			return i;
// 	} while (step < n_probe);
// 	return KMHASH_MAX_SIZE;
// }

struct kmhash_t *init_filter_kmhash(kmint_t size, int n_threads)
{
	struct kmhash_t *h;
	kmint_t i;
	int k;

	h = calloc(1, sizeof(struct kmhash_t));
	h->size = size;
	__round_up_kmint(h->size);
	h->keys = malloc(h->size * sizeof(kmkey_t));
	h->vals = calloc((h->size >> 5), sizeof(kmval_t));
	for (i = 0; i < h->size; ++i)
		h->keys[i] = TOMB_STONE;

	h->n_workers = n_threads;
	// sem_wrap_init(&(h->gsem), n_threads);
	h->status = KMHASH_IDLE;
	h->locks = malloc(n_threads * sizeof(pthread_mutex_t));
	for (k = 0; k < n_threads; ++k)
		pthread_mutex_init(h->locks + k, NULL);

	return h;
}

struct kmhash_t *init_count_kmhash(kmint_t size, int n_threads)
{
	struct kmhash_t *h;
	kmint_t i;
	int k;

	h = calloc(1, sizeof(struct kmhash_t));
	h->size = size;
	h->keys = malloc(h->size * sizeof(kmkey_t));
	h->vals = malloc(h->size * sizeof(kmval_t));
	for (i = 0; i < h->size; ++i) {
		h->keys[i] = TOMB_STONE;
		h->vals[i] = 0;
	}

	h->n_workers = n_threads;
	// sem_wrap_init(&(h->gsem), n_threads);
	h->status = KMHASH_IDLE;
	h->locks = malloc(n_threads * sizeof(pthread_mutex_t));
	for (k = 0; k < n_threads; ++k)
		pthread_mutex_init(h->locks + k, NULL);

	return h;
}

static int is_prime(kmint_t x)
{
	kmint_t i;
	for (i = 2; i * i <= x; ++i)
		if (x % i == 0)
			return 0;
	return 1;
}

kmint_t find_next_prime(kmint_t x)
{
	while (!is_prime(x))
		++x;
	return x;
}

#define rs_set_old(x, i) ((x)[(i) >> 4] = ((x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1)))) |       \
				((uint32_t)KMMASK_OLD << (((i) & 15) << 1)))

#define rs_set_new(x, i) ((x)[(i) >> 4] = ((x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1)))) |       \
				((uint32_t)KMMASK_NEW << (((i) & 15) << 1)))

#define rs_set_empty(x, i) ((x)[(i) >> 4] = (x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1))))

#define rs_get_flag(x, i) (((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)

#define rs_is_old(x, i) ((((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3) == (uint32_t)KMMASK_OLD)

void kmhash_filter_singleton(struct kmhash_t *h, int reinit_val)
{
	kmkey_t *keys;
	kmval_t *vals;

	kmint_t cnt, i, j, new_size, cur_size, step, n_probe;
	kmkey_t x, y, k;
	uint32_t *rs_flag;
	uint32_t cur_flag;

	keys = h->keys;
	vals = h->vals;

	cur_size = h->size;
	// fprintf(stderr, "old_size = %llu\n", h->size);
	// n_probe = h->n_probe;
	// fprintf(stderr, "old_probe = %llu\n", h->n_probe);

	cnt = 0;
	for (i = 0; i < cur_size; ++i) {
		if ((vals[i >> KMVAL_LOG] >> (i & KMVAL_MASK)) & 1) {
			++cnt;
			assert(keys[i] != TOMB_STONE);
		} else
			keys[i] = TOMB_STONE;
	}
	// fprintf(stderr, "cnt = %llu\n", (long long unsigned)cnt);
	// new_size = find_next_prime(cnt * 1.136363);
	new_size = find_next_prime(cnt * 1.432);
	n_probe = estimate_probe_prime(new_size);
	// fprintf(stderr, "probe = %llu\n", n_probe);
	// fprintf(stderr, "new_size = %llu\n", (long long unsigned)new_size);

	rs_flag = calloc((__max(cur_size, new_size) + 15) >> 4, sizeof(uint32_t));
	if (new_size > cur_size) {
		h->keys = keys = realloc(keys, new_size * sizeof(kmkey_t));
		for (i = cur_size; i < new_size; ++i)
			keys[i] = TOMB_STONE;
	}
	for (i = 0; i < cur_size; ++i) {
		if ((vals[i >> KMVAL_LOG] >> (i & KMVAL_MASK)) & 1)
			rs_set_old(rs_flag, i);
	}
	for (i = 0; i < cur_size; ++i) {
		if (rs_is_old(rs_flag, i)) {
			x = keys[i];
			rs_set_empty(rs_flag, i);
			keys[i] = TOMB_STONE;
			while (1) {
				k = __hash_int2(x);
				j = k % new_size;
				// j = x % new_size;
				step = 0;

				cur_flag = KMMASK_NEW;
				while (step <= n_probe) {
					j = (j + step * (step + 1) / 2) % new_size;
					// if (j >= new_size)
					// 	j %= new_size;
					cur_flag = rs_get_flag(rs_flag, j);
					if (cur_flag == KMMASK_EMPTY) {
						rs_set_new(rs_flag, j);
						// rs_flag[j] = KMMASK_NEW;
						keys[j] = x;
						break;
					} else if (cur_flag == KMMASK_OLD) {
						y = keys[j];
						// rs_flag[j] = KMMASK_NEW;
						rs_set_new(rs_flag, j);
						keys[j] = x;
						x = y;
						break;
					}
					++step;
				}
				if (cur_flag == KMMASK_EMPTY)
					break;
				else if (cur_flag == KMMASK_NEW)
					__ERROR("Shrink kmhash size error");
			}
		}
	}
	if (new_size < cur_size)
		h->keys = realloc(h->keys, new_size * sizeof(kmkey_t));
	free(h->vals);
	free(rs_flag);
	h->vals = 0;
	h->size = new_size;
	h->n_items = cnt;
	// cnt = 0;
	// for (i = 0; i < h->size; ++i)
	// 	cnt += (h->keys[i] != TOMB_STONE);
	// if (cnt != h->n_items) {
	// 	fprintf(stderr, "cnt = %llu\nn_items = %llu\n", (long long unsigned)cnt, (long long unsigned)h->n_items);
	// 	assert(cnt == h->n_items);
	// }
	if (reinit_val)
		h->vals = calloc(h->size, sizeof(kmval_t));
}

// void filter_by_count(struct kmhash_t *h, kmval_t min_count)
// {
// 	uint32_t *mask;
// 	kmint_t n_items;
// 	for (i = 0; i < h->size; ++i) {
// 		if (h->keys[i] != TOMB_STONE && h->vals[i] > min_count) {
// 			__set_present(mask, i);
// 			++n_items;
// 		}
// 	}
// 	new_size = n_items;
// 	__round_up_kmint(new_size);
// 	if (n_items > (kmint_t)(0.88 * new_size))
// 		new_size <<= 1;
// 	if (new_size > h->size)
// 		new_size = h->size;
// 	mask_size = (new_size + 15) >> 4;
// 	mask = realloc(mask, mask_size * sizeof(uint32_t));
// 	memset(mask + old_mask_size, (mask_size - old_mask_size) * sizeof(uint32_t));

// 	while (fail) {
// 		for (i = 0; i < cur_size; ++i) {
// 			if (__is_present(mask, i)) {
// 				x = h->keys[i];
// 				__set_empty(
// 			}
// 		}
// 	}
// }

void kmhash_destroy(struct kmhash_t *h)
{
	if (!h) return;
	int i;
	for (i = 0; i < h->n_workers; ++i)
		pthread_mutex_destroy(h->locks + i);
	free(h->locks);
	free(h->keys);
	free(h->vals);
	// sem_wrap_destroy(&(h->gsem));
	free(h);
}
