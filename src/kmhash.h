#ifndef _KMHASH_H_
#define _KMHASH_H_

#include <pthread.h>
#include <stdint.h>

#define HM_MAGIC_1			UINT64_C(0xbf58476d1ce4e5b9)
#define HM_MAGIC_2			UINT64_C(0x94d049bb133111eb)

#define KMHASH_MAX_SIZE			UINT64_C(0x400000000)
#define KMHASH_SINGLE_RESIZE		UINT64_C(0x100000)

#define TOMB_STONE			((kmkey_t)-1)

#define __round_up_kmint(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16,	       \
				 (x) |= (x) >> 32, ++(x))

typedef uint64_t kmint_t;
typedef uint64_t kmkey_t;
typedef uint32_t kmval_t;

#define KMVAL_LOG			5
#define KMVAL_MASK			31

struct kmhash_t {
	kmint_t size;			/* current size */
	kmint_t old_size;		/* previous size */
	kmint_t n_item;			/* number of items */
	kmint_t n_probe;		/* maximum probing times */

	// informative part
	kmkey_t  *keys;			/* keys */
	uint32_t *sgts;			/* count > 1? */
	uint8_t  *adjs;			/* adjacency edges */

	// additional storage for resize
	uint8_t  *flag;

	int status;
	int n_worker;
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

kmint_t hash_get(struct kmhash_t *h, kmkey_t key);

struct kmhash_t *init_sgt_adj_kmhash(kmint_t size, int n_threads);

struct kmhash_t *init_sgt_kmhash(kmint_t size, int n_threads);

void sgt_put(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock);

void sgt_adj_put(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock);

void sgt_add_edge(struct kmhash_t *h, kmkey_t key, int c,
						pthread_mutex_t *lock);

void kmhash_destroy(struct kmhash_t *h);

void sgt_adj_hash_filter(struct kmhash_t *h);

void sgt_hash_filter(struct kmhash_t *h);

#endif
