#ifndef _KMHASH_H_
#define _KMHASH_H_

#include <pthread.h>
#include <stdint.h>

#include "attribute.h"

#define TOMB_STONE			((kmkey_t)-1)

typedef uint64_t kmkey_t;
typedef uint32_t kmval_t;

struct kmhash_t {
	/* multi-threaded singleton kmer filtered hash table */
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

/* Turn on adj bit of a kmer */
void kmhash_add_edge(struct kmhash_t *h, kmkey_t key, int c, pthread_mutex_t *lock);

/*
 * Put a kmer to hash table
 * auto-resize with doubling the adjs array
 */
void kmhash_put_adj(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock);

/*
 * Put a kmer to hash table
 * auto-resize without doubling the adjs array
 */
void kmhash_put(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock);

/*
 * Get the position of a kmer on the table
 * return KMHASH_MAX_SIZE of not found
 */
kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key);

/*
 * Init the hash table with at pre-allocated size items
 * Should know the number of working threads
 */
void kmhash_init(struct kmhash_t *h, kmint_t size, int n_threads, int adj_included);

/* Remove all singleton kmer and shrink the table size */
void kmhash_filter(struct kmhash_t *h, int adj_included);

/* Free all allocated memory on h but not h itself */
void kmhash_destroy(struct kmhash_t *h);

#endif
