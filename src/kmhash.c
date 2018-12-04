#include <stdlib.h>

#include "kmhash.h"
#include "semaphore_wrapper.h"
#include "utils.h"
#include "verbose.h"

#define __round_up_kmint(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16,	       \
				 ++(x))

#define HM_MAGIC_1			UINT64_C(0xbf58476d1ce4e5b9)
#define HM_MAGIC_2			UINT64_C(0x94d049bb133111eb)

#define KMHASH_IDLE			0
#define KMHASH_BUSY			1

struct kmresize_bundle_t {
	struct kmhash_t *h;
	int n_threads;
	int thread_no;
	pthread_barrier_t *barrier;
};

static kmkey_t __hash_int2(kmkey_t k)
{
	kmkey_t x = k;
	x = (x ^ (x >> 30)) * HM_MAGIC_1;
	x = (x ^ (x >> 27)) * HM_MAGIC_2;
	x ^= (x > 31);
	return x;
}

static kmint_t estimate_probe(kmint_t size)
{
	kmint_t s, i;
	i = s = 0;
	while (s < size) {
		++i;
		s += i * i * 2048;
	}
	return i;
}

static kmint_t kmhash_put(struct kmhash_t *h, kmkey_t key)
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
		return i;
	}
	return KMHASH_MAX_SIZE;
}

void *kmresize_worker(void *data)
{
	struct kmresize_bundle_t *bundle = (struct kmresize_bundle_t *)data;
	struct kmhash_t *h;
	kmint_t i, k, l, r, cap;

	h = bundle->h;
	// Init all bucket
	cap = h->size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), h->size);
	for (i = l; i < r; ++i) {
		h->keys[i] = TOMB_STONE;
		h->vals[i] = 0;
	}

	pthread_barrier_wait(bundle->barrier);

	// Fill buckets
	cap = h->old_size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), h->old_size);
	for (i = l; i < r; ++i) {
		if (h->old_keys[i] == TOMB_STONE)
			continue;
		k = kmhash_put(h, h->old_keys[i]);
		if (k == KMHASH_MAX_SIZE)
			__ERROR("Resizing hash table fail");
		__sync_fetch_and_add(&(h->vals[k]), h->old_vals[i]);
	}

	pthread_exit(NULL);
}

void kmhash_resize_multi(struct kmhash_t *h)
{
	int n_threads, i;
	n_threads = h->n_workers;

	h->old_size = h->size;
	h->old_keys = h->keys;
	h->old_vals = h->vals;

	h->size <<= 1;
	h->n_probe = estimate_probe(h->size);
	h->keys = malloc(h->size * sizeof(kmkey_t));
	h->vals = malloc(h->size * sizeof(kmval_t));

	h->n_items = 0;

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

	free(t);
	free(bundles);
	free(h->old_keys);
	free(h->old_vals);
}

void kmhash_resize_single(struct kmhash_t *h)
{
	kmint_t i, k;

	h->old_size = h->size;
	h->old_keys = h->keys;
	h->old_vals = h->vals;

	h->size <<= 1;
	h->n_probe = estimate_probe(h->size);
	h->keys = malloc(h->size * sizeof(kmkey_t));
	h->vals = malloc(h->size * sizeof(kmval_t));

	h->n_items = 0;
	for (i = 0; i < h->size; ++i) {
		h->keys[i] = TOMB_STONE;
		h->vals[i] = 0;
	}

	for (i = 0; i < h->old_size; ++i) {
		if (h->old_keys[i] == TOMB_STONE)
			continue;
		k = kmhash_put(h, h->old_keys[i]);
		if (k == KMHASH_MAX_SIZE)
			__ERROR("Resizing hash table fail");
		__sync_add_and_fetch(&(h->vals[k]), h->old_vals[i]);
	}
	free(h->old_keys);
	free(h->old_vals);
}

void kmhash_resize(struct kmhash_t *h)
{
	int i;
	for (i = 0; i < h->n_workers; ++i)
		// sem_wrap_wait(&(h->gsem));
		pthread_mutex_lock(h->locks + i);

	if (h->size == KMHASH_MAX_SIZE)
		__ERROR("Unable to expand the hash table (max size = %llu)",
			(unsigned long long)KMHASH_MAX_SIZE);

	if (h->size <= KMHASH_SINGLE_RESIZE)
		kmhash_resize_single(h);
	else
		kmhash_resize_multi(h);

	for (i = 0; i < h->n_workers; ++i)
		// sem_wrap_post(&(h->gsem));
		pthread_mutex_unlock(h->locks + i);
}


void kmhash_put_wrap(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	// sem_wrap_wait(&(h->gsem));
	pthread_mutex_lock(lock);
	k = kmhash_put(h, key);
	// sem_wrap_post(&(h->gsem));
	pthread_mutex_unlock(lock);

	if (k == KMHASH_MAX_SIZE) {
		do {
			if (__sync_bool_compare_and_swap(&(h->status), KMHASH_IDLE, KMHASH_BUSY)) {
				kmhash_resize(h);
				__sync_val_compare_and_swap(&(h->status), KMHASH_BUSY, KMHASH_IDLE);
			}

			// sem_wrap_wait(&(h->gsem));
			pthread_mutex_lock(lock);
			k = kmhash_put(h, key);
			// sem_wrap_post(&(h->gsem));
			pthread_mutex_unlock(lock);
		} while (k == KMHASH_MAX_SIZE);
	}
}

void kmhash_inc_val(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	// sem_wrap_wait(&(h->gsem));
	pthread_mutex_lock(lock);
	k = kmhash_put(h, key);
	if (k < KMHASH_MAX_SIZE)
		__sync_add_and_fetch(&(h->vals[k]), 1);
	// sem_wrap_post(&(h->gsem));
	pthread_mutex_unlock(lock);

	if (k == KMHASH_MAX_SIZE) {
		do {
			if (__sync_bool_compare_and_swap(&(h->status), KMHASH_IDLE, KMHASH_BUSY)) {
				kmhash_resize(h);
				__sync_val_compare_and_swap(&(h->status), KMHASH_BUSY, KMHASH_IDLE);
			}

			// sem_wrap_wait(&(h->gsem));
			pthread_mutex_lock(lock);
			k = kmhash_put(h, key);
			if (k < KMHASH_MAX_SIZE)
				__sync_add_and_fetch(&(h->vals[k]), 1);
			// sem_wrap_post(&(h->gsem));
			pthread_mutex_unlock(lock);
		} while (k == KMHASH_MAX_SIZE);
	}
}

kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key)
{
	kmint_t mask, step, i, n_probe;
	kmkey_t k;
	mask = h->size - 1;
	k = __hash_int2(key);
	i = k & mask;
	n_probe = h->n_probe;
	if (h->keys[i] == key)
		return i;
	step = 0;
	do {
		++step;
		i = (i + (step * (step + 1)) / 2) & mask;
		if (h->keys[i] == key)
			return i;
	} while (step < n_probe);
	return KMHASH_MAX_SIZE;
}

struct kmhash_t *init_kmhash(kmint_t size, int n_threads)
{
	struct kmhash_t *h;
	kmint_t i;
	int k;

	h = calloc(1, sizeof(struct kmhash_t));
	h->size = size;
	__round_up_kmint(h->size);
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
