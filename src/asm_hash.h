#ifndef __ASM_HASH_H__
#define __ASM_HASH_H__

#include <stdint.h>

#define HM_COFF_1			UINT64_C(0x87c37b91114253d5)
#define HM_COFF_2			UINT64_C(0x4cf5ad432745937f)

#define HM_MIX_1			UINT64_C(0xbf58476d1ce4e5b9)
#define HM_MIX_2			UINT64_C(0x94d049bb133111eb)

#define KMFLAG_EMPTY			0
#define KMFLAG_OLD			1
#define KMFLAG_NEW			2
#define KMFLAG_LOADING			3

#define KMHASH_IDLE			0
#define KMHASH_BUSY			1

#define KMHASH_MAX_SIZE			UINT64_C(0x400000000)
#define KMHASH_SINGLE_RESIZE		UINT64_C(0x100000)

#define BARCODE_HASH_UPPER	0.77

typedef uint64_t kmint_t;
#define atomic_add_and_fetch_kmint	atomic_add_and_fetch64

#define __round_up_kmint(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16,	       \
				 (x) |= (x) >> 32, ++(x))

#define IDHASH_END(h) ((h)->size)

#define IDHASH_KEY(h, k) ((h)->keys[k])

#define IDHASH_ID(h, k) ((h)->id[k])

#define IDHASH_ADJ(h, k) ((h)->adjs[k])

#define KMHASH_END(h) (KMHASH_MAX_SIZE)

#define KMHASH_KEY(h, k) ((h)->keys[k])

#define BARCODE_HASH_END(h) ((h)->size)

static inline uint64_t __hash_int(uint64_t k)
{
	uint64_t x = k;
	x = (x ^ (x >> 30)) * HM_MIX_1;
	x = (x ^ (x >> 27)) * HM_MIX_2;
	x ^= (x >> 31);
	return x;
}

static inline kmint_t estimate_probe_3(kmint_t size)
{
	kmint_t s, i;
	i = s = 0;
	while (s < size) {
		++i;
		s += i * i * i * 64;
	}
	return i;
}

#endif /* __ASM_HASH_H__ */
