#ifndef _KMHASH_H_
#define _KMHASH_H_

#include <pthread.h>
#include <stdint.h>

#define HM_MAGIC_1			UINT64_C(0xbf58476d1ce4e5b9)
#define HM_MAGIC_2			UINT64_C(0x94d049bb133111eb)

#define KMHASH_MAX_SIZE			UINT64_C(0x80000000)
#define KMHASH_SINGLE_RESIZE		UINT64_C(0x100000)

#define TOMB_STONE			((kmkey_t)-1)

#define __round_up_kmint(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16,	       \
				 ++(x))

typedef uint64_t kmint_t;
typedef uint64_t kmkey_t;
typedef uint32_t kmval_t;

#define KMVAL_LOG			5
#define KMVAL_MASK			31

// struct kmbucket_t {
// 	kmkey_t idx;
// 	kmval_t cnt;
// };

struct kmhash_t {
	kmint_t size;
	kmint_t old_size;
	kmint_t n_items;
	kmint_t n_probe;

	kmkey_t *keys;
	kmkey_t *old_keys;
	kmval_t *vals;
	kmval_t *old_vals;

	// additional storage for resize
	uint8_t *rs_flag;

	int status;
	int n_workers;
	// struct sem_wrap_t gsem;
	/* Each threads need a lock of hashtable access state (busy|idle)
	 * Using a shared semaphore (like in HeraT) causes a big data race and
	 * slow down much, much (100 times compare to HeraT)
	 */
	pthread_mutex_t *locks;
};

static inline kmkey_t __hash_int2(kmkey_t k)
{
	kmkey_t x = k;
	x = (x ^ (x >> 30)) * HM_MAGIC_1;
	x = (x ^ (x >> 27)) * HM_MAGIC_2;
	x ^= (x > 31);
	return x;
}

// struct kmhash_t *init_kmhash(kmint_t size, int n_threads);

struct kmhash_t *init_filter_kmhash(kmint_t size, int n_threads);

void kmhash_destroy(struct kmhash_t *h);

void kmhash_put(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock);

kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key);

kmint_t kmphash_get(struct kmhash_t *h, kmkey_t key);

// void kmhash_put_wrap(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock);

// void kmhash_inc_val(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock);

// kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key);


#endif
