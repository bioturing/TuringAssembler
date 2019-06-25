#ifndef __BARCODE_HASH_H__
#define __BARCODE_HASH_H__

#include <stdint.h>
#include "asm_hash.h"

typedef uint64_t k31key_t;
#define __hash_k31(x)		__hash_int((uint64_t)(x))

struct barcode_hash_t {
	uint32_t size;
	uint32_t n_item;
	uint32_t n_unique;

	uint64_t *keys;
	uint32_t *cnts;
};

void barcode_hash_init(struct barcode_hash_t *h, uint32_t size);

void barcode_hash_clean(struct barcode_hash_t *h);

uint32_t barcode_hash_get(struct barcode_hash_t *h, uint64_t key);

uint32_t barcode_hash_put(struct barcode_hash_t *h, uint64_t key);

uint32_t barcode_hash_add(struct barcode_hash_t *h, uint64_t key);

uint32_t barcode_hash_add_unique(struct barcode_hash_t *h, uint64_t key);

uint32_t barcode_hash_inc_count(struct barcode_hash_t *h, uint64_t key);

void barcode_hash_merge(struct barcode_hash_t *dst, struct barcode_hash_t *src);

void barcode_hash_merge_barcode(struct barcode_hash_t *dst, struct barcode_hash_t *src);

void barcode_hash_merge_readpair(struct barcode_hash_t *dst, struct barcode_hash_t *src);

void barcode_hash_clone(struct barcode_hash_t *dst, struct barcode_hash_t *src);

void barcode_hash_filter(struct barcode_hash_t *h, uint32_t thres);

void barcode_hash_destroy(struct barcode_hash_t *h);

#endif /* __BARCODE_HASH_H__ */
