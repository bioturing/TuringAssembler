#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "atomic.h"
#include "asm_hash.h"
#include "attribute.h"
#include "io_utils.h"
#include "kmhash.h"
#include "utils.h"
#include "verbose.h"

static inline uint64_t __kmhash_get_hash(const uint8_t *s, int word_size)
{
	uint64_t ret = 0;
	int i;
	for (i = 0; i < word_size; ++i)
		ret = (ret << 5) - ret + (uint64_t)s[i];
	return ret;
}

#if defined(_MSC_VER)

#define FORCE_INLINE	__forceinline

#include <stdlib.h>

#define ROTL32(x,y)	_rotl(x,y)
#define ROTL64(x,y)	_rotl64(x,y)

#define BIG_CONSTANT(x) (x)

// Other compilers

#else	// defined(_MSC_VER)

#define	FORCE_INLINE inline __attribute__((always_inline))

inline uint32_t rotl32(uint32_t x, int8_t r)
{
	return (x << r) | (x >> (32 - r));
}

inline uint64_t rotl64(uint64_t x, int8_t r)
{
	return (x << r) | (x >> (64 - r));
}

#define	ROTL32(x,y)	rotl32(x,y)
#define ROTL64(x,y)	rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)

#endif // !defined(_MSC_VER)

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

FORCE_INLINE uint32_t getblock32(const uint32_t * p, int i)
{
	return p[i];
}

FORCE_INLINE uint64_t getblock64(const uint64_t * p, int i)
{
	return p[i];
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

FORCE_INLINE uint32_t fmix32(uint32_t h)
{
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;

	return h;
}

//----------

FORCE_INLINE uint64_t fmix64(uint64_t k)
{
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xff51afd7ed558ccd);
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
	k ^= k >> 33;

	return k;
}

//-----------------------------------------------------------------------------


static inline uint64_t MurmurHash3_x64_64(const uint8_t *data, const int len)
{
	int n_blocks = len >> 4;
	uint64_t h1 = BIG_CONSTANT(0x13097);
	uint64_t h2 = BIG_CONSTANT(0x13097);

	const uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
	const uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

	const uint64_t *blocks = (const uint64_t *)(data);

	int i;
	for (i = 0; i < n_blocks; ++i) {
		uint64_t k1 = getblock64(blocks, i << 1);
		uint64_t k2 = getblock64(blocks, (i << 1) + 1);
		k1 *= c1; k1 = ROTL64(k1, 31); k1 *= c2; h1 ^= k1;
		h1 = ROTL64(h1, 27); h1 += h2; h1 = h1 * 5 + 0x52dce729;
		k2 *= c2; k2 = ROTL64(k2, 33); k2 *= c1; h2 ^= k2;
		h2 = ROTL64(h2, 31); h2 += h1; h2 = h2 * 5 + 0x38495ab5;
	}

	const uint8_t *tail = (const uint8_t *)(data + (n_blocks << 4));
	uint64_t k1 = 0;
	uint64_t k2 = 0;

	switch (len & 15) {
	case 15: k2 ^= ((uint64_t)tail[14]) << 48;
	case 14: k2 ^= ((uint64_t)tail[13]) << 40;
	case 13: k2 ^= ((uint64_t)tail[12]) << 32;
	case 12: k2 ^= ((uint64_t)tail[11]) << 24;
	case 11: k2 ^= ((uint64_t)tail[10]) << 16;
	case 10: k2 ^= ((uint64_t)tail[ 9]) << 8;
	case  9: k2 ^= ((uint64_t)tail[ 8]) << 0;
		k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

	case  8: k1 ^= ((uint64_t)tail[ 7]) << 56;
	case  7: k1 ^= ((uint64_t)tail[ 6]) << 48;
	case  6: k1 ^= ((uint64_t)tail[ 5]) << 40;
	case  5: k1 ^= ((uint64_t)tail[ 4]) << 32;
	case  4: k1 ^= ((uint64_t)tail[ 3]) << 24;
	case  3: k1 ^= ((uint64_t)tail[ 2]) << 16;
	case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;
	case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;
		k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
	}

	h1 ^= len; h2 ^= len;

	h1 += h2;
	h2 += h1;

	h1 = fmix64(h1);
	h2 = fmix64(h2);

	h1 += h2;
	h2 += h1;

	return h2;
}

#define __kmhash_equal(a, b, ws) (memcmp((const char *)(a), (const char *)(b), ws) == 0)

kmint_t kmhash_get(struct kmhash_t *h, const uint8_t *key)
{
	kmint_t n_probe, mask, step, i;
	mask = h->size - 1;
	n_probe = h->n_probe;
	int word_size = h->word_size;
	uint64_t k = MurmurHash3_x64_64(key, word_size);
	i = k & mask;
	step = 0;
	do {
		if (h->flag[i] == KMFLAG_EMPTY)
			return KMHASH_END(h);
		if (__kmhash_equal(h->keys + i * word_size, key, word_size))
			return i;
		++step;
		i = (i + step * (step + 1) / 2) & mask;
	} while (step <= n_probe);
	return KMHASH_END(h);
}

static kmint_t internal_kmhash_put_multi(struct kmhash_t *h, const uint8_t *key)
{
	kmint_t mask, step, i, last, n_probe;
	mask = h->size - 1;
	int word_size = h->word_size;
	uint64_t k = MurmurHash3_x64_64(key, word_size);
	uint8_t cur_flag;
	n_probe = h->n_probe;
	i = k & mask;
	step = 0;
	do {
		cur_flag = atomic_val_CAS8(&(h->flag[i]), KMFLAG_EMPTY,
								KMFLAG_LOADING);
		if (cur_flag == KMFLAG_EMPTY) {
			memcpy(h->keys + i * word_size, key, word_size);
			atomic_add_and_fetch_kmint(&(h->n_item), 1);
			h->flag[i] = KMFLAG_NEW;
			return i;
		} else if (cur_flag == KMFLAG_LOADING) {
			continue;
		} else {
			if (__kmhash_equal(h->keys + i * word_size, key, word_size))
				return i;
			++step;
			i = (i + step * (step + 1) / 2) & mask;
		}
	} while (step <= n_probe);
	return KMHASH_END(h);
}

static kmint_t internal_kmhash_put(struct kmhash_t *h, const uint8_t *key)
{
	kmint_t mask, step, i, last, n_probe;
	mask = h->size - 1;
	int word_size = h->word_size;
	uint64_t k = MurmurHash3_x64_64(key, word_size);
	n_probe = h->n_probe;
	i = k & mask;
	step = 0;
	do {
		if (h->flag[i] == KMFLAG_EMPTY) {
			memcpy(h->keys + i * word_size, key, word_size);
			++h->n_item;
			h->flag[i] = KMFLAG_NEW;
			return i;
		} else {
			if (__kmhash_equal(h->keys + i * word_size, key, word_size))
				return i;
			++step;
			i = (i + step * (step + 1) / 2) & mask;
		}
	} while (step <= n_probe);
	return KMHASH_END(h);
}

static void kmhash_resize(struct kmhash_t *h)
{
	__VERBOSE("Resize ---- size = %lu; n_item = %lu\n", h->size, h->n_item);
	kmint_t old_size, size, mask, n_probe, i;
	int word_size = h->word_size;
	old_size = h->size;
	h->size <<= 1;
	size = h->size;
	h->upper_bound = h->size * HASH_SIZE_UPPER;
	mask = size - 1;
	h->n_probe = estimate_probe_3(h->size);
	n_probe = h->n_probe;
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

	if (aux_flag & KM_AUX_POS) {
		h->pos = realloc(h->pos, h->size * sizeof(struct edge_data_t));
		memset(h->pos + old_size, 0, (size - old_size) * sizeof(struct edge_data_t));
	}

	uint8_t *keys = h->keys;
	uint8_t *adjs = h->adjs;
	uint8_t *flag = h->flag;
	gint_t *idx = h->idx;
	struct edge_data_t *pos = h->pos;

	for (i = 0; i < old_size; ++i) {
		if (flag[i] == KMFLAG_NEW)
			flag[i] = KMFLAG_OLD;
	}

	uint8_t *x, *xt;
	x = alloca(word_size);
	xt = alloca(word_size);
	uint8_t a = 0, at;
	gint_t id = 0, idt;
	struct edge_data_t p = (struct edge_data_t){NULL, 0}, pt;
	for (i = 0; i < old_size; ++i) {
		if (flag[i] == KMFLAG_OLD) {
			memcpy(x, keys + i * word_size, word_size);
			if (aux_flag & KM_AUX_ADJ) {
				a = adjs[i];
				adjs[i] = 0;
			}
			if (aux_flag & KM_AUX_IDX)
				id = idx[i];
			if (aux_flag & KM_AUX_POS) {
				p = pos[i];
				pos[i] = (struct edge_data_t){(void *)NULL, (gint_t)0};
			}
			flag[i] = KMFLAG_EMPTY;
			while (1) {
				uint64_t k = MurmurHash3_x64_64(x, word_size);
				kmint_t j, step, last;
				j = k & mask;
				// last = j;
				step = 0;
				while (step <= n_probe) {
					if (flag[j] == KMFLAG_EMPTY || flag[j] == KMFLAG_OLD)
						break;
					++step;
					j = (j + step * (step + 1) / 2) & mask;
				}
				// while (flag[j] != KMFLAG_EMPTY && flag[j] != KMFLAG_OLD) {
				// 	j = (j + (++step)) & mask;
				// 	if (j == last)
				// 		break;
				// }
				if (flag[j] == KMFLAG_EMPTY) {
					flag[j] = KMFLAG_NEW;
					memcpy(keys + j * word_size, x, word_size);
					if (aux_flag & KM_AUX_ADJ)
						adjs[j] = a;
					if (aux_flag & KM_AUX_IDX)
						idx[j] = id;
					if (aux_flag & KM_AUX_POS)
						pos[j] = p;
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
					if (aux_flag & KM_AUX_POS) {
						pt = pos[j];
						pos[j] = p;
						p = pt;
					}
				} else {
					__ERROR("Resize kmhash failed");
				}
			}
		}
	}
}

static void kmhash_resize_multi(struct kmhash_t *h)
{
}

void kmhash_put_multi(struct kmhash_t *h, const uint8_t *key, pthread_mutex_t *lock)
{
	kmint_t k;

	pthread_mutex_lock(lock);
	k = internal_kmhash_put_multi(h, key);
	pthread_mutex_unlock(lock);

	while (k == KMHASH_END(h)) {
		if (atomic_bool_CAS8(&(h->status), KMHASH_IDLE, KMHASH_BUSY)) {
			kmhash_resize_multi(h);
			atomic_val_CAS8(&(h->status), KMHASH_BUSY, KMHASH_IDLE);
		}

		pthread_mutex_lock(lock);
		k = internal_kmhash_put_multi(h, key);
		pthread_mutex_unlock(lock);
	}
}

kmint_t kmhash_put(struct kmhash_t *h, const uint8_t *key)
{
	// if (h->n_item >= h->upper_bound)
	// 	kmhash_resize(h);
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
	if (aux_flag & KM_AUX_POS) {
		if (h->aux_flag & KM_AUX_POS)
			free(h->pos);
		h->pos = calloc(h->size, sizeof(struct edge_data_t));
		h->aux_flag |= KM_AUX_POS;
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
	if (aux_flag & KM_AUX_POS)
		h->pos = calloc(h->size, sizeof(struct edge_data_t));
	else
		h->pos = NULL;
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
