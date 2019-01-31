#ifndef __K63_HASH_H__
#define __K63_HASH_H__

#include <stdint.h>

#include "attribute.h"

typedef struct {
	uint64_t bin[2];
} k63key_t;

#define k63_equal(x, y) ((x).bin[0] == (y).bin[0] && (x).bin[1] == (y).bin[1])

static inline uint64_t __hash_k63(k63key_t x)
{
	uint64_t k1, k2, k;
	k1 = x.bin[0]; k2 = x.bin[1];
	k1 *= HM_COFF_1; k1 = __rotl64(k1, 31); k1 *= HM_COFF_2;
	k1 = __rotl64(k1, 27); k1 = k1 * 5 + 0x52dce729;
	k2 *= HM_COFF_2; k2 = __rotl64(k2, 33); k2 *= HM_COFF_1;
	k2 = __rotl64(k2, 31); k2 = k2 * 5 + 0x38495ab5;
	k = k1 + k2;
	k = __hash_int(k);
	return k;
}

struct k63hash_t {
	/* multi-threaded singleton kmer filtered hash table */
	kmint_t size;			/* current size */
	kmint_t old_size;		/* previous size */
	kmint_t n_item;			/* number of items */
	kmint_t n_probe;		/* maximum probing times */

	/* informative part */
	k63key_t *keys;			/* keys */
	uint32_t *sgts;			/* count > 1? */
	uint8_t  *adjs;			/* adjacency edges */
	uint8_t  *flag;			/* bucket empty or not */

	int status;
	int n_worker;
	/* Each threads need a lock of hashtable access state (busy|idle)
	 * Using a shared semaphore (like in HeraT) causes a big data race and
	 * slow down much, much (100 times compare to HeraT)
	 */
	pthread_mutex_t *locks;
};

/* Turn on adj bit of a kmer */
void k63hash_add_edge(struct k63hash_t *h, k63key_t key, int c, pthread_mutex_t *lock);

/*
 * Put a kmer to hash table
 * auto-resize with doubling the adjs array
 */
void k63hash_put_adj(struct k63hash_t *h, k63key_t key, pthread_mutex_t *lock);

/*
 * Put a kmer to hash table
 * auto-resize without doubling the adjs array
 */
void k63hash_put(struct k63hash_t *h, k63key_t key, pthread_mutex_t *lock);

/*
 * Get the position of a kmer on the table
 * return KMHASH_MAX_SIZE of not found
 */
kmint_t k63hash_get(struct k63hash_t *h, k63key_t key);

/*
 * Init the hash table with at pre-allocated size items
 * Should know the number of working threads
 */
void k63hash_init(struct k63hash_t *h, kmint_t size, int n_threads, int adj_included);

/* Remove all singleton kmer and shrink the table size */
void k63hash_filter(struct k63hash_t *h, int adj_included);

/* Free all allocated memory on h but not h itself */
void k63hash_destroy(struct k63hash_t *h);

#endif
