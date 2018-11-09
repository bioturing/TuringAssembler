#include <stdlib.h>

#include "kmhash.h"
#include "utils.h"
#include "verbose.h"

#define __hash_int(x) (int32_t)((x) >> 33 ^ (x) ^ (x) << 11)

#define __round_up_kmint(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16,	       \
				 ++(x))

#define HM_MAGIC_1                      UINT64_C(0xbf58476d1ce4e5b9)
#define HM_MAGIC_2                      UINT64_C(0x94d049bb133111eb)

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
		s += i * i * 128;
	}
	return i;
}

static kmint_t kmhash_put(struct kmhash_t *h, kmkey_t key)
{
	kmint_t mask, step, i, n_probe;
	kmkey_t cur_key, k, tombstone;

	mask = h->size - 1;
	k = __hash_int2(key);
	n_probe = h->n_probe;
	tombstone = (kmkey_t)-1;
	i = k & mask;
	{
		cur_key = __sync_val_compare_and_swap(&(h->bucks[i].idx), tombstone, key);
		if (cur_key == tombstone || cur_key == key) {
			if (cur_key == tombstone)
				__sync_fetch_and_add(&(h->n_items), 1);
			return i;
		}
	}
	step = 0;
	do {
		i = (i + (step * (step + 1)) / 2) & mask;
		++step;
		cur_key = __sync_val_compare_and_swap(&(h->bucks[i].idx), tombstone, key);
	} while (step < n_probe && cur_key != key && cur_key != tombstone);
	if (cur_key == tombstone || cur_key == key) {
		if (cur_key == tombstone)
			__sync_fetch_and_add(&(h->n_items), 1);
		return i;
	} else {
		return h->size;
	}
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
	kmkey_t tombstone;
	tombstone = (kmkey_t)-1;
	for (i = l; i < r; ++i) {
		h->bucks[i].idx = tombstone;
		h->bucks[i].cnt = 0;
	}

	pthread_barrier_wait(bundle->barrier);

	cap = h->old_size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), h->old_size);
	for (i = l; i < r; ++i) {
		if (h->old_bucks[i].idx == tombstone)
			continue;
		k = kmhash_put(h, h->old_bucks[i].idx);
		if (k == h->size)
			__ERROR("Resizing kmer hash table fail");
		__sync_add_and_fetch(&(h->bucks[k].cnt), h->old_bucks[i].cnt);
	}

	pthread_exit(NULL);
}

void kmhash_enlarge(struct kmhash_t *h)
{
        __VERBOSE("I need to double hash table size\n");
	int i;
	for (i = 0; i < h->n_workers; ++i)
		pthread_mutex_lock(h->locks + i);

	if (h->size == KMHASH_MAX_SIZE)
		__ERROR("The kmer hash table is too big (exceeded %llu)", (unsigned long long)KMHASH_MAX_SIZE);

	h->old_size = h->size;
	h->old_bucks = h->bucks;

	h->size <<= 1;
	h->bucks = malloc(h->size * sizeof(struct kmbucket_t));
	// FIXME: calculate new prob value
	h->n_probe = estimate_probe(h->size);
        __VERBOSE("New probe: %d; New size: %llu\n", (int)h->n_probe, (unsigned long long)h->size);
        __VERBOSE("Current number of items: %llu\n", (unsigned long long)h->n_items);
        h->n_items = 0;

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_t *th;
	th = calloc(h->n_workers, sizeof(pthread_t));

	struct kmresize_bundle_t *bundles;
	bundles = calloc(h->n_workers, sizeof(struct kmresize_bundle_t));

	pthread_barrier_t barrier;
	pthread_barrier_init(&barrier, NULL, h->n_workers);

	for (i = 0; i < h->n_workers; ++i) {
		bundles[i].n_threads = h->n_workers;
		bundles[i].thread_no = i;
		bundles[i].barrier = &barrier;
		bundles[i].h = h;
		pthread_create(th + i, &attr, kmresize_worker, bundles + i);
	}

	for (i = 0; i < h->n_workers; ++i)
		pthread_join(th[i], NULL);

        for (i = 0; i < h->n_workers; ++i)
                pthread_mutex_unlock(h->locks + i);

	pthread_barrier_destroy(&barrier);
	pthread_attr_destroy(&attr);

	free(th);
	free(bundles);
	free(h->old_bucks);
        __sync_bool_compare_and_swap(&(h->status), KMHASH_RESIZE, KMHASH_IDLE);
        __VERBOSE("Done resize\n");
}

void kmhash_put_wrap(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	pthread_mutex_lock(lock);
	k = kmhash_put(h, key);
	pthread_mutex_unlock(lock);

	if (k == h->size) {
		if (__sync_bool_compare_and_swap(&(h->status), KMHASH_IDLE, KMHASH_RESIZE))
			kmhash_enlarge(h);

		pthread_mutex_lock(lock);
		k = kmhash_put(h, key);
		pthread_mutex_unlock(lock);
	}
}

void kmhash_inc_val_wrap(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock)
{
	kmint_t k;

	pthread_mutex_lock(lock);
	k = kmhash_put(h, key);
	__sync_add_and_fetch(&(h->bucks[k].cnt), 1);
	pthread_mutex_unlock(lock);

	if (k == h->size) {
		if (__sync_bool_compare_and_swap(&(h->status), KMHASH_IDLE, KMHASH_RESIZE))
			kmhash_enlarge(h);

		pthread_mutex_lock(lock);
		k = kmhash_put(h, key);
		__sync_add_and_fetch(&(h->bucks[k].cnt), 1);
		pthread_mutex_unlock(lock);
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
	if (h->bucks[i].idx == key)
		return i;
	step = 1;
	do {
		i = (i + (step * (step + 1)) / 2) & mask;
		if (h->bucks[i].idx == key)
			return i;
		++step;
	} while (step < n_probe);
	return h->size;
}

struct kmhash_t *init_kmhash(kmint_t size, int n_threads)
{
        __VERBOSE("Initilizing hash table\n");
	struct kmhash_t *h;
	kmint_t i;
	kmkey_t tombstone;
	h = calloc(1, sizeof(struct kmhash_t));
	h->size = size;
	__round_up_kmint(h->size);
	h->bucks = malloc(h->size * sizeof(struct kmbucket_t));
	h->n_probe = estimate_probe(h->size);
        __VERBOSE("Probe number: %d\n", (int)h->n_probe);

	tombstone = (kmkey_t)-1;
	for (i = 0; i < h->size; ++i)
		h->bucks[i].idx = tombstone;

	h->n_workers = n_threads;
	h->locks = calloc(n_threads, sizeof(pthread_mutex_t));
	int k;
	for (k = 0; k < n_threads; ++k)
		pthread_mutex_init(h->locks + k, NULL);
	h->status = KMHASH_IDLE;

	return h;
}

void kmhash_destroy(struct kmhash_t *h)
{
	if (!h) return;
	free(h->bucks);
	int i;
	for (i = 0; i < h->n_workers; ++i)
		pthread_mutex_destroy(h->locks + i);
	free(h->locks);
	free(h);
}
