#ifndef _KMHASH_H_
#define _KMHASH_H_

#include <pthread.h>
#include <semaphore.h>
#include <stdint.h>

#define KMHASH_IDLE			0
#define KMHASH_RESIZE			1

typedef uint32_t kmint_t;
#define KMHASH_MAX_SIZE			UINT64_C(0x80000000)
typedef uint64_t kmkey_t;
typedef uint64_t kmval_t;

struct kmbucket_t {
	kmkey_t idx;
	kmval_t cnt;
};

struct kmhash_t {
	kmint_t size;
	kmint_t old_size;
	kmint_t n_items;
	kmint_t n_probe;
	struct kmbucket_t *bucks;
	struct kmbucket_t *old_bucks;
	int status;
	int n_workers;
	pthread_mutex_t *locks;
};

struct kmresize_bundle_t {
	struct kmhash_t *h;
	int n_threads;
	int thread_no;
	pthread_barrier_t *barrier;
};

struct kmhash_t *init_kmhash(kmint_t size, int n_threads);

void kmhash_destroy(struct kmhash_t *h);

void kmhash_put_wrap(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock);

void kmhash_inc_val_wrap(struct kmhash_t *h, kmkey_t key, pthread_mutex_t *lock);

kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key);


#endif
