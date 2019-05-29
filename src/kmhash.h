#ifndef __KMHASH_H__
#define __KMHASH_H__

#include <stdint.h>

#include "attribute.h"

#define KM_AUX_ADJ			0x1
#define KM_AUX_IDX			0x2

#define KM_IS_EMPTY(h, k) ((h)->flag[k] == KMFLAG_EMPTY)

#define KMHASH_END(h) ((h)->size)

#define KMHASH_KEY(h, k) ((h)->keys + ((size_t)(k) * (h)->word_size))

#define KMHASH_ADJ(h, k) ((h)->adjs[k])

#define KMHASH_IDX(h, k) ((h)->idx[k])

struct kmhash_t {
	kmint_t size;
	kmint_t upper_bound;
	kmint_t n_item;
	kmint_t n_probe;

	uint8_t *keys;
	uint8_t *flag;

	uint8_t *adjs;
	gint_t  *idx;

	int word_size;
	uint32_t aux_flag;
};

void kmhash_init(struct kmhash_t *h, kmint_t pre_alloc, int ksize,
							uint32_t aux_flag);

void kmhash_alloc_aux(struct kmhash_t *h, uint32_t aux_flag);

kmint_t kmhash_get(struct kmhash_t *h, const uint8_t *key);

kmint_t kmhash_put(struct kmhash_t *h, const uint8_t *key);

#define kmhash_set_adj(h, k, c) ((h)->adjs[k] |= (uint8_t)1 << (c))

void kmhash_destroy(struct kmhash_t *h);

#endif /* __KMHASH_H__ */
