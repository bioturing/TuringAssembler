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
	h->keys = malloc(h->size * sizeof(uint64_t));
	h->cnts = calloc(h->size, sizeof(uint32_t));
	uint32_t i;
	for (i = 0; i < h->size; ++i)
		h->keys[i] = (uint64_t)-1;
}

void barcode_hash_clean(struct barcode_hash_t *h)
{
	free(h->keys);
	free(h->cnts);
	h->keys = NULL;
	h->cnts = NULL;
	h->n_item = h->size = 0;
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

uint32_t barcode_hash_put(struct barcode_hash_t *h, uint64_t key)
{
	uint32_t k;
	// pthread_mutex_lock(h->lock);
	k = internal_barcode_hash_put(h, key);
	while (k == BARCODE_HASH_END(h)) {
		barcode_hash_resize(h);
		k = internal_barcode_hash_put(h, key);
	}
	// pthread_mutex_unlock(h->lock);
	return k;
}

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

void barcode_hash_merge(struct barcode_hash_t *dst, struct barcode_hash_t *src)
{
	uint32_t k, j;
	for (k = 0; k < src->size; ++k) {
		if (src->keys[k] == K31_NULL)
			continue;
		j = barcode_hash_put(dst, src->keys[k]);
		dst->cnts[j] += src->cnts[k];
	}
}

void barcode_hash_clone(struct barcode_hash_t *dst, struct barcode_hash_t *src)
{
	dst->size = src->size;
	dst->n_item = src->n_item;
	dst->keys = malloc(dst->size * sizeof(uint64_t));
	dst->cnts = malloc(dst->size * sizeof(uint32_t));
	memcpy(dst->keys, src->keys, dst->size * sizeof(uint64_t));
	memcpy(dst->cnts, src->cnts, dst->size * sizeof(uint32_t));
}
