#ifndef __K31_HASH_H__
#define __K31_HASH_H__

#include <stdint.h>

#include "atomic.h"
#include "asm_hash.h"

typedef uint64_t k31key_t;

#define atomic_val_CAS_k31key		atomic_val_CAS64

#define __hash_k31(x)		__hash_int((uint64_t)(x))

#define k31_exist(h, i) ((h)->keys[i] != (k31key_t)-1)

#define K31_NULL		((k31key_t)-1)

struct k31hash_t {
	/* multi-threaded singleton kmer filtered hash table */
	kmint_t size;			/* current size */
	kmint_t old_size;		/* previous size */
	kmint_t n_item;			/* number of items */
	kmint_t n_probe;		/* maximum probing times */

	/* informative part */
	k31key_t *keys;			/* keys */
	uint32_t *sgts;			/* count > 1? */
	uint8_t  *adjs;			/* adjacency edges */
	uint8_t  *flag;

	int status;
	int n_worker;
	/* Each threads need a lock of hashtable access state (busy|idle)
	 * Using a shared semaphore (like in HeraT) causes a big data race and
	 * slow down much, much (100 times compare to HeraT)
	 */
	pthread_mutex_t *locks;
};

/* Save binary hash table to file */
void save_k31hash(struct k31hash_t *h, const char *path);

/* Load saved binary hash table from file */
void load_k31hash(struct k31hash_t *h, const char *path);

/* Turn on adj bit of a kmer */
void k31hash_add_edge(struct k31hash_t *h, k31key_t key, int c, pthread_mutex_t *lock);

/*
 * Put a kmer to hash table
 * auto-resize with doubling the adjs array
 */
void k31hash_put_adj(struct k31hash_t *h, k31key_t key, pthread_mutex_t *lock);

/*
 * Put a kmer to hash table
 * auto-resize without doubling the adjs array
 */
void k31hash_put(struct k31hash_t *h, k31key_t key, pthread_mutex_t *lock);

/*
 * Get the position of a kmer on the table
 * return KMHASH_MAX_SIZE of not found
 */
kmint_t k31hash_get(struct k31hash_t *h, k31key_t key);

/*
 * Init the hash table with at pre-allocated size items
 * Should know the number of working threads
 */
void k31hash_init(struct k31hash_t *h, kmint_t size, int n_threads, int adj_included);

/* Remove all singleton kmer and shrink the table size */
void k31hash_filter(struct k31hash_t *h, int adj_included);

/* Free all allocated memory on h but not h itself */
void k31hash_destroy(struct k31hash_t *h);

/****************************** Single thread k31 hash ************************/

struct k31_idhash_t {
	kmint_t size;		/* maximum current size */
	kmint_t n_probe;	/* probe time */
	kmint_t n_item;		/* number of current items */

	k31key_t *keys;		/* keys */
	gint_t   *id;		/* unique id */
	uint8_t  *adjs;		/* adjacency character(s) */
};

void k31_convert_table(struct k31hash_t *h, struct k31_idhash_t *dict);

void k31_idhash_init(struct k31_idhash_t *dict, kmint_t size, int adj_included);

void k31_idhash_clean(struct k31_idhash_t *h);

kmint_t k31_idhash_get(struct k31_idhash_t *h, k31key_t key);

kmint_t k31_idhash_put(struct k31_idhash_t *h, k31key_t key);

#endif
