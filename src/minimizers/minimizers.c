//
// Created by BioTuring on 2019-11-08.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <inttypes.h>

#include "assembly_graph.h"
#include "khash.h"
#include "attribute.h"
#include "utils.h"
#include "minimizers.h"

const char *bit_rep[16] = {
	[ 0] = "0000", [ 1] = "0001", [ 2] = "0010", [ 3] = "0011",
	[ 4] = "0100", [ 5] = "0101", [ 6] = "0110", [ 7] = "0111",
	[ 8] = "1000", [ 9] = "1001", [10] = "1010", [11] = "1011",
	[12] = "1100", [13] = "1101", [14] = "1110", [15] = "1111",
};

#define PRINT_BYTE(byte) printf("%s%s\n", bit_rep[(byte >> 4)& 0x0F], bit_rep[byte & 0x0F]);
#define PRINT_UINT64(byte) printf("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n", bit_rep[(byte >> 60) & 0x0F], bit_rep[(byte >> 56)& 0x0F], \
								bit_rep[(byte >> 52) & 0x0F], bit_rep[(byte >> 48) & 0x0F], \
								bit_rep[(byte >> 44) & 0x0F], bit_rep[(byte >> 40) & 0x0F], \
								bit_rep[(byte >> 36) & 0x0F], bit_rep[(byte >> 32) & 0x0F], \
								bit_rep[(byte >> 28) & 0x0F], bit_rep[(byte >> 24) & 0x0F], \
								bit_rep[(byte >> 20) & 0x0F], bit_rep[(byte >> 16) & 0x0F], \
								bit_rep[(byte >> 12) & 0x0F], bit_rep[(byte >> 8) & 0x0F],  \
								bit_rep[(byte >> 4) & 0x0F], bit_rep[(uint8_t)byte & 0x0F]);


KHASH_MAP_INIT_INT64(mm_hash, uint32_t);        /* Hash table structure of the minimizers */



#define BARCODES100M 100663320
#define BIG_CONSTANT(x) (x##LLU)

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

static inline uint64_t fmix64(uint64_t k)
{
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xff51afd7ed558ccd);
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
	k ^= k >> 33;

	return k;
}

static inline uint64_t getblock64(const uint64_t * p, int i)
{
	return p[i];
}

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

static inline char *uin32t2seq(uint32_t *s, int l)
{
	char *seq;
	int i;
	seq = calloc(l, sizeof(char));
	for (i = 0; i < l; ++i)
		seq[i] = nt4_char[__binseq_get(s, i)];
	return seq;
}

void mm_print(struct mm_db_t *db)
{
	int i, j, pad;
	uint64_t km;
	uint8_t c;
	for (i = 0 ; i < db->n; ++i) {
		km = db->mm[i];
		printf("km[%d]:", i);
		PRINT_UINT64(km);
		for (j = 0 ; j < db->k; ++j) {
			pad = (62 - 2*j); /* 64 bits/number and 2bits/char */
			c = (km & ((uint64_t)3 << pad)) >> pad;
			printf("%c", nt4_char[c]);
		}
		printf("\n");
	}
}

void mm_db_insert(struct mm_db_t *db, uint64_t km, uint32_t p)
{
	if (db->n == db->size) {
		db->mm =  realloc(db->mm, (db->size << 1) * sizeof(uint64_t));
		db->p =  realloc(db->p, (db->size << 1) * sizeof(uint64_t));
		db->size <<= 1;
	}
	db->mm[db->n] = km;
	db->p[db->n++] = p;
}
struct mm_db_t * mm_db_init()
{
	struct mm_db_t *db = calloc(1, sizeof(struct mm_db_t));
	db->n = 0;
	db->size = 8;
	db->mm = calloc(db->size, sizeof(struct mm_db_t));
	db->p = calloc(db->size, sizeof(struct mm_db_t));
	return db;
}

static inline uint64_t get_km_i(uint32_t *s, int i, int k)
{
	uint64_t km, c;
	int j;
	int pad = (32 - k - 1)*2;
	km = 0;
	for (j = 0; j < k; ++j) {
		c = (uint64_t)__binseq_get(s, i + j);
		km |= c;
		km <<= 2;
	}
	km <<= pad;
	return km;
}

#define HASH64(k) MurmurHash3_x64_64((uint8_t *)&k, sizeof(uint64_t)/sizeof(uint8_t))
#define DEBUG_PRINT printf

struct mm_db_t * mm_index_str(uint32_t *s, int k, int w, int l)
{
	struct mm_db_t *db = mm_db_init();
	db->k = k;

	int i, j, p = -1;
	uint64_t km, mm, c;
	uint64_t km_h, mm_h;
	int pad = (32 - k - 1)*2;

	mm_h = km_h = HASH64(k);
	for (i = 0; i < l - w + 1; ++i) {
		DEBUG_PRINT("[i = %d]\n", i);
		if (i + w + k - 1 >= l)
			break;
		if (p < i) {
			km = mm = get_km_i(s, i, k);
			mm_h = km_h = HASH64(mm);
			p = i;
			for (j = 0; j < w; ++j) {
				c = (uint64_t) __binseq_get(s, i + j + k - 1);
				km |= ((uint64_t) c << (pad + 2));
				km_h = HASH64(km);
				if (km_h < mm_h) {
					mm = km;
					mm_h = km_h;
					p = i + j;
					mm_db_insert(db, km, p);
				}
				km <<= 2;
			}
			DEBUG_PRINT("[1]minimizers at window %d: %d\n", i, p);

			continue;
		} else {
			c = (uint64_t) __binseq_get(s, i + w + k - 2);
			km |= ((uint64_t) c << (pad + 2));
			km_h = HASH64(km);
			if (km_h < mm_h){
				p = i + w - 1;
				mm = km;
				mm_h = km_h;
				mm_db_insert(db, km, p);
				DEBUG_PRINT("[2]minimizers at window %d: %d\n", i, p);
			}
			km <<= 2;
		}
	}
	//mm_print(db);
	return db;
}
