#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "attribute.h"
#include "barcode_hash.h"
#include "io_utils.h"
#include "utils.h"
#include "verbose.h"

#define K31_NULL		((k31key_t)-1)

void barcode_hash_init(struct barcode_hash_t *h, uint32_t size)
{
	h->size = size - 1;
	__round_up_32(h->size);
	h->n_item = 0;
	// h->n_unique = 0;
	h->keys = malloc(h->size * sizeof(uint64_t));
	h->cnts = calloc(h->size, sizeof(uint32_t));
	uint32_t i;
	for (i = 0; i < h->size; ++i)
		h->keys[i] = (uint64_t)-1;
}

void barcode_hash_clean(struct barcode_hash_t *h)
{
	if (h != NULL) {
		if (h->keys != NULL) free(h->keys);
		if (h->cnts != NULL) free(h->cnts);
		h->keys = NULL;
		h->cnts = NULL;
		h->n_item = h->size = h->n_unique = 0;
	}
}

uint32_t barcode_hash_get(struct barcode_hash_t *h, uint64_t key)
{
	uint32_t i, last, mask, step = 0;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_int(key); i = k & mask;
	last = i;
	while (h->keys[i] != K31_NULL && h->keys[i] != key) {
		i = (i + (++step)) & mask;
		if (i == last) return BARCODE_HASH_END(h);
	}
	return h->keys[i] == key ? i : BARCODE_HASH_END(h);
}

static inline uint32_t internal_barcode_hash_put(struct barcode_hash_t *h, uint64_t key)
{
	if (h->n_item >= (uint32_t)(h->size * BARCODE_HASH_UPPER))
		return BARCODE_HASH_END(h);
	uint32_t i, last, mask, step = 0;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_int(key); i = k & mask;
	if (h->keys[i] == (uint64_t)-1) {
		h->keys[i] = key;
		++h->n_item;
		return i;
	}
	last = i;
	while (h->keys[i] != (uint64_t)-1 && h->keys[i] != key) {
		i = (i + (++step)) & mask;
		if (i == last)
			break;
	}
	if (h->keys[i] == (uint64_t)-1 || h->keys[i] == key) {
		if (h->keys[i] == (uint64_t)-1) {
			h->keys[i] = key;
			++h->n_item;
		}
		return i;
	} else {
		return BARCODE_HASH_END(h);
	}
}

static void barcode_hash_resize(struct barcode_hash_t *h)
{
	uint32_t old_size, mask, i;
	old_size = h->size;
	h->size <<= 1;
	mask = h->size - 1;
	h->keys = realloc(h->keys, h->size * sizeof(uint64_t));
	h->cnts = realloc(h->cnts, h->size * sizeof(uint32_t));
	uint8_t *flag = calloc(h->size, sizeof(uint8_t));

	for (i = old_size; i < h->size; ++i) {
		h->keys[i] = K31_NULL;
		h->cnts[i] = 0;
	}

	for (i = 0; i < old_size; ++i) {
		if (h->keys[i] != K31_NULL)
			flag[i] = KMFLAG_OLD;
	}

	for (i = 0; i < old_size; ++i) {
		if (flag[i] == KMFLAG_OLD) {
			uint64_t x = h->keys[i], xt;
			uint32_t y = h->cnts[i], yt;
			flag[i] = KMFLAG_EMPTY;
			h->keys[i] = K31_NULL;
			h->cnts[i] = 0;
			while (1) {
				uint64_t k = __hash_k31(x);
				uint32_t j = k & mask, step = 0, last;
				last = j;
				while (flag[j] != KMFLAG_EMPTY && flag[j] != KMFLAG_OLD) {
					j = (j + (++step)) & mask;
					if (j == last)
						break;
				}
				if (flag[j] == KMFLAG_EMPTY) {
					flag[j] = KMFLAG_NEW;
					h->keys[j] = x;
					h->cnts[j] = y;
					break;
				} else if (flag[j] == KMFLAG_OLD) {
					flag[j] = KMFLAG_NEW;
					xt = h->keys[j];
					yt = h->cnts[j];
					h->keys[j] = x;
					h->cnts[j] = y;
					x = xt; y = yt;
				} else {
					__ERROR("Resize barcode hash failed");
				}
			}
		}
	}
	free(flag);
}

void barcode_hash_filter(struct barcode_hash_t *h, uint32_t thres)
{
	uint32_t n_item, i, j, new_size, cur_size, mask, step;
	cur_size = h->size;
	uint8_t *flag = calloc(cur_size, sizeof(uint8_t));
	n_item = 0;
	for (i = 0; i < cur_size; ++i) {
		if (h->keys[i] != K31_NULL && h->cnts[i] > thres) {
			++n_item;
			flag[i] = KMFLAG_OLD;
		} else {
			h->keys[i] = K31_NULL;
			h->cnts[i] = 0;
			flag[i] = KMFLAG_EMPTY;
		}
	}
	new_size = n_item;
	__round_up_32(new_size);
	if (new_size > cur_size)
		new_size = cur_size;
	mask = new_size - 1;

	uint64_t x = K31_NULL, xt;
	uint32_t y = 0, yt;
	int retry = 0;
loop_refill:
	if (retry) {
		fprintf(stderr, "retry shrink #%d\n", retry);
		for (i = 0; i < cur_size; ++i) {
			if (flag[i] == KMFLAG_NEW)
				flag[i] = KMFLAG_OLD;
		}
		for (i = 0; i < cur_size; ++i) {
			if (flag[i] == KMFLAG_EMPTY) {
				h->keys[i] = x;
				h->cnts[i] = y;
				flag[i] = KMFLAG_OLD;
				break;
			}
		}
		new_size <<= 1;
		mask = new_size - 1;
	}
	for (i = 0; i < cur_size; ++i) {
		if (flag[i] == KMFLAG_OLD) {
			x = h->keys[i];
			y = h->cnts[i];
			h->keys[i] = K31_NULL;
			h->cnts[i] = 0;
			flag[i] = KMFLAG_EMPTY;
			while (1) {
				uint64_t k = __hash_k31(x);
				uint32_t j = k & mask, step = 0, last;
				last = j;
				while (flag[j] != KMFLAG_EMPTY && flag[j] != KMFLAG_OLD) {
					j = (j + (++step)) & mask;
					if (j == last)
						break;
				}
				if (flag[j] == KMFLAG_EMPTY) {
					flag[j] = KMFLAG_NEW;
					h->keys[j] = x;
					h->cnts[j] = y;
					break;
				} else if (flag[j] == KMFLAG_OLD) {
					flag[j] = KMFLAG_NEW;
					xt = h->keys[j];
					yt = h->cnts[j];
					h->keys[j] = x;
					h->cnts[j] = y;
					x = xt; y = yt;
				} else {
					if (cur_size == new_size)
						__ERROR("Resize barcode hash failed");
					++retry;
					goto loop_refill;
				}
			}
		}
	}
	h->size = new_size;
	h->n_item = n_item;
	if (new_size < cur_size) {
		h->keys = realloc(h->keys, new_size * sizeof(uint64_t));
		h->cnts = realloc(h->cnts, new_size * sizeof(uint32_t));
	}
	free(flag);
}

uint32_t barcode_hash_add(struct barcode_hash_t *h, uint64_t key)
{
	uint32_t k;
	k = internal_barcode_hash_put(h, key);
	while (k == BARCODE_HASH_END(h)) {
		barcode_hash_resize(h);
		k = internal_barcode_hash_put(h, key);
	}
	++h->cnts[k];
	return k;
}

// uint32_t barcode_hash_add_unique(struct barcode_hash_t *h, uint64_t key)
// {
// 	uint32_t k;
// 	k = internal_barcode_hash_put(h, key);
// 	while (k == BARCODE_HASH_END(h)) {
// 		barcode_hash_resize(h);
// 		k = internal_barcode_hash_put(h, key);
// 	}
// 	if (h->cnts[k] == 0) {
// 		++h->n_unique;
// 		h->cnts[k] = 1;
// 	}
// 	return k;
// }

uint32_t barcode_hash_inc_count(struct barcode_hash_t *h, uint64_t key)
{
	uint32_t k;
	// pthread_mutex_lock(h->lock);
	k = internal_barcode_hash_put(h, key);
	while (k == BARCODE_HASH_END(h)) {
		barcode_hash_resize(h);
		k = internal_barcode_hash_put(h, key);
	}
	++h->cnts[k];
	// pthread_mutex_unlock(h->lock);
	return k;
}

void barcode_hash_clone(struct barcode_hash_t *dst, struct barcode_hash_t *src)
{
	dst->size = src->size;
	dst->n_item = src->n_item;
	// dst->n_unique = src->n_unique;
	dst->keys = malloc(dst->size * sizeof(uint64_t));
	dst->cnts = NULL;
	// dst->cnts = malloc(dst->size * sizeof(uint32_t));
	memcpy(dst->keys, src->keys, dst->size * sizeof(uint64_t));
	// memcpy(dst->cnts, src->cnts, dst->size * sizeof(uint32_t));
}

// void barcode_hash_merge_barcode(struct barcode_hash_t *dst, struct barcode_hash_t *src)
// {
// 	uint32_t i, k;
// 	for (i = 0; i < src->size; ++i) {
// 		if (src->keys[i] == K31_NULL)
// 			continue;
// 		k = internal_barcode_hash_put(dst, src->keys[i]);
// 		// if (src->cnts[i] == 1)
// 		// 	dst->cnts[k] = 1;
// 	}
// }

// void barcode_hash_merge_readpair(struct barcode_hash_t *dst, struct barcode_hash_t *src)
// {
// 	uint32_t i, k;
// 	for (i = 0; i < src->size; ++i) {
// 		if (src->keys[i] == K31_NULL)
// 			continue;
// 		k = internal_barcode_hash_put(dst, src->keys[i]);
// 		// dst->cnts[k] += src->cnts[i];
// 	}
// }

void barcode_hash_merge(struct barcode_hash_t *dst, struct barcode_hash_t *src)
{
	uint32_t i, k;
	for (i = 0; i < src->size; ++i) {
		if (src->keys[i] == K31_NULL)
			continue;
		k = internal_barcode_hash_put(dst, src->keys[i]);
		// if (src->cnts[i] == 1)
		// 	dst->cnts[k] = 1;
	}
}


void barcode_hash_destroy(struct barcode_hash_t *h)
{
	free(h->keys);
	free(h->cnts);
	h->keys = NULL;
	h->cnts = NULL;
	h->size = h->n_item = h->n_unique = 0;
}
