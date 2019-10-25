#ifndef __KMHASH_H__
#define __KMHASH_H__

#include <stdint.h>

#include "attribute.h"

#define KM_AUX_ADJ			0x1
#define KM_AUX_IDX			0x2
#define KM_AUX_POS			0x4

#define KM_IS_EMPTY(h, k) ((h)->flag[k] == KMFLAG_EMPTY)

#define KMHASH_END(h) (KMHASH_MAX_SIZE)

#define KMHASH_KEY(h, k) ((h)->keys + ((size_t)(k) * (h)->word_size))

#define KMHASH_ADJ(h, k) ((h)->adjs[k])

#define KMHASH_IDX(h, k) ((h)->idx[k])

#define KMHASH_POS(h, k) ((h)->pos[k])

struct edge_idx_t {
	gint_t idx;
	gint_t pos;
};

struct edge_data_t {
	struct edge_idx_t *e;
	gint_t n;
};

struct kmhash_t {
	kmint_t size;
	kmint_t upper_bound;
	kmint_t n_item;
	kmint_t n_probe;

	uint8_t *keys;
	uint8_t *flag;

	uint8_t *adjs;
	gint_t  *idx;
	gint_t *kmer_pos;
	struct edge_data_t *pos;

	int word_size;
	uint32_t aux_flag;

	uint8_t status;
	int n_worker;
	pthread_mutex_t *locks;
};

void kmhash_init(struct kmhash_t *h, kmint_t pre_alloc, int ksize,
					uint32_t aux_flag, int n_threads);

void kmhash_alloc_aux(struct kmhash_t *h, uint32_t aux_flag);

kmint_t kmhash_get(struct kmhash_t *h, const uint8_t *key);

kmint_t kmhash_put(struct kmhash_t *h, const uint8_t *key);

void kmhash_put_multi(struct kmhash_t *h, const uint8_t *key, pthread_mutex_t *lock);

#define kmhash_set_adj(h, k, c) ((h)->adjs[k] |= (uint8_t)1 << (c))

void kmhash_set_adj_multi(struct kmhash_t *h, const uint8_t *key, int c, pthread_mutex_t *lock);

void kmhash_set_idx_multi(struct kmhash_t *h, const uint8_t *key, gint_t id, pthread_mutex_t *lock);

void kmhash_set_pos_multi(struct kmhash_t *h, const uint8_t *key, int pos, gint_t id,
		pthread_mutex_t *lock);

void kmhash_destroy(struct kmhash_t *h);

uint64_t MurmurHash3_x64_64(const uint8_t *data, const int len);
#endif /* __KMHASH_H__ */
