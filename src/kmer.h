#ifndef __KMER_H__
#define __KMER_H__

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static inline void km_shift_left(uint8_t *kmer, int len, int l, int shift)
{
	int k, i, cr_len;
	k = ((len - 1) & 3) + 1;
	cr_len = 8 - shift;
	uint8_t cr, tcr;
	cr = 0;
	for (i = 0; i < l; ++i) {
		tcr = kmer[i] >> cr_len;
		kmer[i] = (kmer[i] << shift) | cr;
		cr = tcr;
	}
	kmer[l - 1] &= ((uint8_t)1 << (k << 1)) - 1;
}

static inline void km_shift_right(uint8_t *kmer, int len, int l, int shift)
{
	int i, cr_len;
	uint8_t shift_mask = ((uint8_t)1 << (shift << 1)) - 1;
	uint8_t cr, tcr;
	cr_len = 8 - shift;
	cr = 0;
	for (i = l - 1; i >= 0; --i) {
		tcr = kmer[i] & shift_mask;
		kmer[i] = (kmer[i] >> shift) | (cr << cr_len);
		cr = tcr;
	}
}

static inline void km_shift_left2(uint8_t *kmer, int len, int l)
{
	int k, i;
	k = ((len - 1) & 3) + 1;
	uint8_t cr, tcr;
	cr = 0;
	for (i = 0; i < l; ++i) {
		tcr = kmer[i] >> 6;
		kmer[i] = (kmer[i] << 2) | cr;
		cr = tcr;
	}
	kmer[l - 1] &= ((uint8_t)1 << (k << 1)) - 1;
}

static inline void km_shift_right2(uint8_t *kmer, int len, int l)
{
	int i;
	uint8_t cr, tcr;
	cr = 0;
	for (i = l - 1; i >= 0; --i) {
		tcr = kmer[i] & 0x3;
		kmer[i] = (kmer[i] >> 2) | (cr << 6);
		cr = tcr;
	}
}

static inline void km_shift_append(uint8_t *kmer, int len, int l, uint8_t c)
{
	km_shift_left2(kmer, len, l);
	kmer[0] |= c;
}

static inline void km_shift_append_rv(uint8_t *kmer, int len, int l, uint8_t c)
{
	km_shift_right2(kmer, len, l);
	int k = (len - 1) & 3;
	kmer[l - 1] |= c << (k << 1);
}

/* get left shift right
 * get right and
 */

static inline void kedge_get_left(uint8_t *dst, uint8_t *kmer, int len, int l)
{
	int i;
	uint8_t cr;
	if (l < ((len + 4) >> 2))
		cr = kmer[l] & 0x3;
	else
		cr = 0;
	for (i = l - 1; i >= 0; --i) {
		dst[i] = (kmer[i] >> 2) | (cr << 6);
		cr = kmer[i] & 0x3;
	}
}

static inline void kedge_get_right(uint8_t *dst, uint8_t *kmer, int len, int l)
{
	int i, k;
	k = ((len - 1) & 3) + 1;
	memcpy(dst, kmer, l);
	dst[l - 1] &= ((uint8_t)1 << (k << 1)) - 1;
}

static inline int km_cmp(uint8_t *k1, uint8_t *k2, int l)
{
	int i;
	for (i = l - 1; i >= 0; --i) {
		if (k1[i] < k2[i])
			return -1;
		if (k1[i] > k2[i])
			return 1;
	}
	return 0;
}

static inline void km_get_rc(uint8_t *dst, uint8_t *kmer, int len, int l)
{
	int k, i;
	k = len & 3;
	for (i = 0; i < l; ++i) {
		dst[i] = kmer[l - i - 1];
		dst[i] = (dst[i] >> 4) | (dst[i] << 4);
		dst[i] = ((dst[i] & 0xCC) >> 2) | ((dst[i] & 0x33) << 2);
		dst[i] ^= 0xFF;
	}
	km_shift_right(dst, len, l, ((4 - k) & 3) << 1);
}

#endif /* __KMER_H__ */
