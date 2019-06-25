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

#define HASH_SIZE_UPPER		0.88
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
		s += i * i * i;
	}
	return i;
}
#define __reverse_bit_order64(x)					       \
(									       \
	(x) = (((x) & 0xffffffff00000000ull) >> 32) | (((x) & 0x00000000ffffffffull) << 32), \
	(x) = (((x) & 0xffff0000ffff0000ull) >> 16) | (((x) & 0x0000ffff0000ffffull) << 16), \
	(x) = (((x) & 0xff00ff00ff00ff00ull) >>  8) | (((x) & 0x00ff00ff00ff00ffull) <<  8), \
	(x) = (((x) & 0xf0f0f0f0f0f0f0f0ull) >>  4) | (((x) & 0x0f0f0f0f0f0f0f0full) <<  4), \
	(x) = (((x) & 0xccccccccccccccccull) >>  2) | (((x) & 0x3333333333333333ull) <<  2)  \
)

#endif /* __ASM_HASH_H__ */
