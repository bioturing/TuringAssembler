#ifndef _KMHASH_H_
#define _KMHASH_H_

#include <pthread.h>
#include <stdint.h>

#include "semaphore_wrapper.h"

#define KMHASH_MAX_SIZE			UINT64_C(0x80000000)
#define KMHASH_SINGLE_RESIZE		UINT64_C(0x100000)

#define TOMB_STONE			((kmkey_t)-1)

typedef uint64_t kmint_t;
typedef uint64_t kmkey_t;
typedef uint32_t kmval_t;

// struct kmbucket_t {
// 	kmkey_t idx;
// 	kmval_t cnt;
// };

struct kmhash_t {
	kmint_t size;
	kmint_t old_size;
	kmint_t n_items;
	kmint_t n_probe;
	// struct kmbucket_t *bucks;
	// struct kmbucket_t *old_bucks;
	kmkey_t *keys;
	kmkey_t *old_keys;
	kmval_t *vals;
	kmval_t *old_vals;
	int status;
	int n_workers;
	struct sem_wrap_t gsem;
};

struct kmhash_t *init_kmhash(kmint_t size, int n_threads);

void kmhash_destroy(struct kmhash_t *h);

void kmhash_put_wrap(struct kmhash_t *h, kmkey_t key);

void kmhash_inc_val(struct kmhash_t *h, kmkey_t key);

kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key);


#endif
