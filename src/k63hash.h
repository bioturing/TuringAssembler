#ifndef __K63_HASH_H__
#define __K63_HASH_H__

#include <stdint.h>

#include "asm_hash.h"
#include "attribute.h"
#include "utils.h"

typedef struct {
	uint64_t bin[2];
} k63key_t;

#define __k63_equal(x, y) ((x).bin[0] == (y).bin[0] && (x).bin[1] == (y).bin[1])

#define __k63_lt(x, y) ((x).bin[1] < (y).bin[1] || ((x).bin[1] == (y).bin[1] && (x).bin[0] < (y).bin[0]))

#define __k63_lshift2(k) (((k).bin[1] = ((k).bin[1] << 2) | ((k).bin[0] >> 62)), \
				((k).bin[0] <<= 2))
#define __k63_rshift2(k) (((k).bin[0] = ((k).bin[0] >> 2) | (((k).bin[1] & 0x3ull) << 62)), \
				((k).bin[1] >>= 2))

#define __k63_lshift(k, l) (((k).bin[1] = ((k).bin[1] << (l)) | ((k).bin[0] >> (64 - (l)))), \
				((k).bin[0] <<= (l)))
#define __k63_rshift(k, l) (((k).bin[0] = ((k).bin[0] >> (l)) | (((k).bin[1] & ((1ull << (l)) - 1)) << (64 - (l))), \
				((k).bin[1] >>= (l))))

#define __k63_and(k, v) ((k).bin[0] &= (v).bin[0], (k).bin[1] &= (v).bin[1])

#define __k63_revc_num(y, x, l, mask)					       \
(									       \
	(x) = (y), __k63_lshift(x, 128 - ((l) << 1)),			       \
	__reverse_bit_order64((x).bin[0]), __reverse_bit_order64((x).bin[1]),  \
	(x).bin[0] ^= (x).bin[1],					       \
	(x).bin[1] ^= (x).bin[0],					       \
	(x).bin[0] ^= (x).bin[1],					       \
	(x).bin[0] ^= 0xffffffffffffffffull, (x).bin[0] &= (mask).bin[0],      \
	(x).bin[1] ^= 0xffffffffffffffffull, (x).bin[1] &= (mask).bin[1]       \
)

#define __k63_rev_num(y, x, l)	\
(									       \
	(x) = (y), __k63_lshift(x, 128 - ((l) << 1)),			       \
	(x).bin[0] ^= (x).bin[1],					       \
	(x).bin[1] ^= (x).bin[0],					       \
	(x).bin[0] ^= (x).bin[1],					       \
	__reverse_bit_order64((x).bin[0]), __reverse_bit_order64((x).bin[1])   \
)

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

/* Save binary hash table to file */
void save_k63hash(struct k63hash_t *h, const char *path);

/* Load saved binary hash table from file */
void load_k63hash(struct k63hash_t *h, const char *path);

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

/****************************** Single thread k63 hash ************************/

struct k63_idhash_t {
	kmint_t size;		/* maximum current size */
	kmint_t n_probe;	/* probe time */
	kmint_t n_item;		/* number of current items */

	k63key_t *keys;		/* keys */
	gint_t   *id;		/* unique id */
	uint8_t  *adjs;		/* adjacency character(s) */
	uint8_t  *flag;		/* mark available key */
};

void k63_convert_table(struct k63hash_t *h, struct k63_idhash_t *dict);

void k63_idhash_init(struct k63_idhash_t *dict, kmint_t size, int adj_included);

void k63_idhash_clean(struct k63_idhash_t *h);

kmint_t k63_idhash_get(struct k63_idhash_t *h, k63key_t key);

kmint_t k63_idhash_put(struct k63_idhash_t *h, k63key_t key);

#endif
