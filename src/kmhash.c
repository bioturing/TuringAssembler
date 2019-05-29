#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "asm_hash.h"
#include "attribute.h"
#include "io_utils.h"
#include "kmhash.h"
#include "utils.h"
#include "verbose.h"

static inline uint64_t __kmhash_get_hash(int word_size, const uint8_t *s)
{
	uint64_t ret = 0;
	int i;
	for (i = 0; i < word_size; ++i)
		ret = (ret << 5) - ret + (uint64_t)s[i];
	return ret;
}

#define __kmhash_equal(ws, a, b) (memcmp((const char *)(a), (const char *)(b), ws) == 0)

kmint_t kmhash_get(struct kmhash_t *h, const uint8_t *key)
{
	kmint_t n_probe, mask, step, i, last;
	mask = h->size - 1;
	n_probe = h->n_probe;
	int word_size = h->word_size;
	uint64_t k = __kmhash_get_hash(word_size, key);
	i = k & mask;
	last = i;
	step = 0;
	while (h->flag[i] != KMFLAG_EMPTY &&
		!__kmhash_equal(word_size, h->keys + i * word_size, key)) {
		i = (i + (++step)) & mask;
		if (i == last) return h->size;
	}
	if (h->flag[i] != KMFLAG_EMPTY)
		assert(__kmhash_equal(word_size, h->keys + i * word_size, key));
	return (h->flag[i] == KMFLAG_EMPTY) ? h->size : i;

	// step = 0;
	// while (step <= n_probe && h->flag[i] != KMFLAG_EMPTY &&
	// 	!__kmhash_equal(word_size, h->keys + i * word_size, key)) {
	// 	++step;
	// 	i = (i + step * (step + 1) / 2) & mask;
	// }
	// return (step <= n_probe && h->flag[i] != KMFLAG_EMPTY) ? i : h->size;
}

static kmint_t internal_kmhash_put(struct kmhash_t *h, const uint8_t *key)
{
	kmint_t n_probe, mask, step, i, last;
	mask = h->size - 1;
	// n_probe = h->n_probe;
	int word_size = h->word_size;
	uint64_t k = __kmhash_get_hash(word_size, key);
	// __VERBOSE("hash_key = %lu\n", k);
	i = k & mask;
	last = i;
	step = 0;
	while (h->flag[i] != KMFLAG_EMPTY &&
		!__kmhash_equal(word_size, h->keys + i * word_size, key)) {
		i = (i + (++step)) & mask;
		if (i == last) return h->size;
	}
	if (h->flag[i] == KMFLAG_EMPTY) {
		h->flag[i] = KMFLAG_NEW;
		memcpy(h->keys + i * word_size, key, word_size);
		++h->n_item;
	}
	return i;

	// step = 0;
	// while (step <= n_probe && h->flag[i] != KMFLAG_EMPTY &&
	// 	!__kmhash_equal(word_size, h->keys + i * word_size, key)) {
	// 	++step;
	// 	i = (i + step * (step + 1) / 2) & mask;
	// }
	// if (step <= n_probe) {
	// 	if (h->flag[i] == KMFLAG_EMPTY) {
	// 		h->flag[i] = KMFLAG_NEW;
	// 		memcpy(h->keys + i * word_size, key, word_size);
	// 		++h->n_item;
	// 	}
	// 	return i;
	// }
	// return h->size;
}

static void kmhash_resize(struct kmhash_t *h)
{
	// __VERBOSE("Resize ---- size = %lu; n_item = %lu\n", h->size, h->n_item);
	kmint_t old_size, size, mask, n_probe, i;
	int word_size = h->word_size;
	old_size = h->size;
	h->size <<= 1;
	size = h->size;
	h->upper_bound = h->size * HASH_SIZE_UPPER;
	mask = size - 1;
	// h->n_probe = estimate_probe_3(h->size);
	// n_probe = h->n_probe;
	uint32_t aux_flag = h->aux_flag;

	size_t new_byte_size = h->size * h->word_size;
	h->keys = realloc(h->keys, new_byte_size);
	h->flag = realloc(h->flag, h->size * sizeof(uint8_t));
	memset(h->flag + old_size, 0, (size - old_size) * sizeof(uint8_t));
	if (aux_flag & KM_AUX_ADJ) {
		h->adjs = realloc(h->adjs, h->size * sizeof(uint8_t));
		memset(h->adjs + old_size, 0, (size - old_size) * sizeof(uint8_t));
	}
	if (aux_flag & KM_AUX_IDX)
		h->idx = realloc(h->idx, h->size * sizeof(gint_t));

	uint8_t *keys = h->keys;
	uint8_t *adjs = h->adjs;
	uint8_t *flag = h->flag;
	gint_t *idx = h->idx;

	for (i = 0; i < old_size; ++i) {
		if (flag[i] == KMFLAG_NEW)
			flag[i] = KMFLAG_OLD;
	}

	uint8_t *x, *xt;
	x = alloca(word_size);
	xt = alloca(word_size);
	uint8_t a = 0, at;
	gint_t id = 0, idt;
	// uint8_t hash_sum = 0;
	// for (i = 0; i < old_size; ++i) {
	// 	if (h->flag[i] == KMFLAG_EMPTY)
	// 		continue;
	// 	int k;
	// 	for (k = 0; k < word_size; ++k)
	// 		hash_sum ^= h->keys[i * word_size + k];
	// 	hash_sum ^= h->adjs[i];
	// }
	for (i = 0; i < old_size; ++i) {
		if (flag[i] == KMFLAG_OLD) {
			memcpy(x, keys + i * word_size, word_size);
			if (aux_flag & KM_AUX_ADJ) {
				a = adjs[i];
				adjs[i] = 0;
			}
			if (aux_flag & KM_AUX_IDX)
				id = idx[i];
			flag[i] = KMFLAG_EMPTY;
			while (1) {
				uint64_t k = __kmhash_get_hash(word_size, x);
				kmint_t j, step, last;
				j = k & mask;
				last = j;
				step = 0;
				while (flag[j] != KMFLAG_EMPTY && flag[j] != KMFLAG_OLD) {
					j = (j + (++step)) & mask;
					if (j == last)
						break;
				}
				if (flag[j] == KMFLAG_EMPTY) {
					flag[j] = KMFLAG_NEW;
					memcpy(keys + j * word_size, x, word_size);
					if (aux_flag & KM_AUX_ADJ)
						adjs[j] = a;
					if (aux_flag & KM_AUX_IDX)
						idx[j] = id;
					break;
				} else if (flag[j] == KMFLAG_OLD) {
					flag[j] = KMFLAG_NEW;
					memcpy(xt, keys + j * word_size, word_size);
					memcpy(keys + j * word_size, x, word_size);
					memcpy(x, xt, word_size);
					if (aux_flag & KM_AUX_ADJ) {
						at = adjs[j];
						adjs[j] = a;
						a = at;
					}
					if (aux_flag & KM_AUX_IDX) {
						idt = idx[j];
						idx[j] = id;
						id = idt;
					}
				} else {
					__ERROR("Resize kmhash failed");
				}
			}
		}
	}
	// uint8_t hash_sum2 = 0;
	// for (i = 0; i < size; ++i) {
	// 	if (h->flag[i] == KMFLAG_EMPTY)
	// 		continue;
	// 	int k;
	// 	for (k = 0; k < word_size; ++k)
	// 		hash_sum2 ^= h->keys[i * word_size + k];
	// 	hash_sum2 ^= h->adjs[i];
	// }
	// assert(hash_sum == hash_sum2);
	// fprintf(stderr, "hash_sum = %hhu\n", hash_sum);
	// fprintf(stderr, "hash_sum2 = %hhu\n", hash_sum2);

	// kmint_t recount = 0;
	// for (i = 0; i < h->size; ++i) {
	// 	assert(h->flag[i] == KMFLAG_EMPTY || h->flag[i] == KMFLAG_NEW);
	// 	recount += (h->flag[i] == KMFLAG_NEW);
	// }
	// assert(recount == h->n_item);
	// for (i = 0; i < h->size; ++i) {
	// 	if (h->flag[i] == KMFLAG_EMPTY)
	// 		continue;
	// 	kmint_t k = kmhash_get(h, h->keys + i * word_size);
	// 	if (k != i) {
	// 		for (int x = 0; x < word_size; ++x)
	// 			__VERBOSE(x + 1 == word_size ? "%hhu\n" : "%hhu ", h->keys[i * word_size + x]);
	// 		for (int x = 0; x < word_size; ++x)
	// 			__VERBOSE(x + 1 == word_size ? "%hhu\n" : "%hhu ", h->keys[k * word_size + x]);
	// 		__VERBOSE("k = %lu; i = %lu; h->size = %lu\n", k, i, h->size);
	// 		char *t1 = h->keys + i * word_size;
	// 		char *t2 = h->keys + k * word_size;
	// 		__VERBOSE("t1 = %lu; t2 = %lu\n", t1, t2);
	// 		__VERBOSE("strcmp = %d\n", strncmp(t1, t2, word_size));
	// 		__VERBOSE("strncmp = %d\n", __kmhash_equal(word_size, h->keys + i * word_size, h->keys + k * word_size));
	// 		__VERBOSE("flag[i] = %hhu; flag[k] = %hhu\n", h->flag[i], h->flag[k]);
	// 		assert(0);
	// 	}
	// 	assert(k == i);
	// }
}

kmint_t kmhash_put(struct kmhash_t *h, const uint8_t *key)
{
	if (h->n_item >= h->upper_bound)
		kmhash_resize(h);
	kmint_t k;
	k = internal_kmhash_put(h, key);
	int n_try = 0;
	while (k == KMHASH_END(h)) {
		++n_try;
		if (n_try > 3)
			__ERROR("Insert kmhash failed");
		kmhash_resize(h);
		k = internal_kmhash_put(h, key);
	}
	return k;
}

void kmhash_alloc_aux(struct kmhash_t *h, uint32_t aux_flag)
{
	if (aux_flag & KM_AUX_ADJ) {
		if (h->aux_flag & KM_AUX_ADJ)
			free(h->adjs);
		h->adjs = calloc(h->size, sizeof(uint8_t));
		h->aux_flag |= KM_AUX_ADJ;
	}
	if (aux_flag & KM_AUX_IDX) {
		if (h->aux_flag & KM_AUX_IDX)
			free(h->idx);
		h->idx = malloc(h->size * sizeof(gint_t));
		h->aux_flag |= KM_AUX_IDX;
	}
}

void kmhash_init(struct kmhash_t *h, kmint_t pre_alloc, int ksize, uint32_t aux_flag)
{
	h->size = pre_alloc;
	__round_up_kmint(h->size);
	h->n_probe = estimate_probe_3(h->size);
	h->n_item = 0;
	h->upper_bound = h->size * HASH_SIZE_UPPER;
	h->aux_flag = aux_flag;
	h->word_size = (ksize + 3) >> 2;
	size_t l = h->size * h->word_size;
	h->keys = malloc(l);
	h->flag = calloc(h->size, sizeof(uint8_t));
	if (aux_flag & KM_AUX_ADJ)
		h->adjs = calloc(h->size, sizeof(uint8_t));
	else
		h->adjs = NULL;
	if (aux_flag & KM_AUX_IDX)
		h->idx = malloc(h->size * sizeof(gint_t));
	else
		h->idx = NULL;
}

void kmhash_destroy(struct kmhash_t *h)
{
	free(h->keys);
	free(h->flag);
	free(h->adjs);
	free(h->idx);
	h->keys = NULL;
	h->flag = NULL;
	h->adjs = NULL;
	h->idx = NULL;
}
