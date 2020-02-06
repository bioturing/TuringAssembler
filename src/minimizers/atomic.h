#ifndef __ATOMIC_H__
#define __ATOMIC_H__

#include <stdint.h>

#define atomic_add_and_fetch32(a, x)		__sync_add_and_fetch(a, x)
#define atomic_add_and_fetch64(a, x)		__sync_add_and_fetch(a, x)
#define atomic_sub_and_fetch64(a, x)		__sync_sub_and_fetch(a, x)

#define atomic_val_CAS8(ptr, old_val, new_val)				       \
			__sync_val_compare_and_swap(ptr, old_val, new_val)
#define atomic_val_CAS32(ptr, old_val, new_val)				       \
			__sync_val_compare_and_swap(ptr, old_val, new_val)
#define atomic_val_CAS64(ptr, old_val, new_val)				       \
			__sync_val_compare_and_swap(ptr, old_val, new_val)
#define atomic_bool_CAS8(ptr, old_val, new_val)			       \
			__sync_bool_compare_and_swap(ptr, old_val, new_val)
#define atomic_bool_CAS32(ptr, old_val, new_val)			       \
			__sync_bool_compare_and_swap(ptr, old_val, new_val)
#define atomic_bool_CAS64(ptr, old_val, new_val)			       \
			__sync_bool_compare_and_swap(ptr, old_val, new_val)

#define get_bit32(x, i) (((x) >> (i)) & 1)
#define set_bit32(x, i) ((x) | ((uint32_t)1 << (i)))
#define clear_bit32(x, i) ((x) & (~((uint32_t)1 << (i))))
/* i-th bit must be zero */
#define set_bit_val32(x, i, b) ((x) | ((uint32_t)(b) << (i)))
/* any i-th bit */
#define set_bit_var32(x, i, b) (((x) & (~((uint32_t)1 << (i)))) | ((uint32_t)(b) << (i)))

#define atomic_get_bit32(ptr, i) (((*(volatile uint32_t *)(ptr)) >> (i)) & 1)
#define atomic_get_bit8(ptr, i) (((*(volatile uint8_t *)(ptr)) >> (i)) & 1)

static inline uint32_t atomic_set_bit32(uint32_t *ptr, int i)
{
	uint32_t old_bin, new_bin, cur_bin;
	cur_bin = *(volatile uint32_t *)ptr;
	do {
		old_bin = cur_bin;
		new_bin = old_bin | ((uint32_t)1 << i);
		cur_bin = atomic_val_CAS32(ptr, old_bin, new_bin);
	} while (cur_bin != old_bin);
	return cur_bin;
}

static inline uint8_t atomic_set_bit8(uint8_t *ptr, int i)
{
	uint8_t old_bin, new_bin, cur_bin;
	cur_bin = *(volatile uint8_t *)ptr;
	do {
		old_bin = cur_bin;
		new_bin = old_bin | ((uint8_t)1 << i);
		cur_bin = atomic_val_CAS8(ptr, old_bin, new_bin);
	} while (cur_bin != old_bin);
	return cur_bin;
}

static inline uint32_t atomic_clear_bit32(uint32_t *ptr, int i)
{
	uint32_t old_bin, new_bin, cur_bin;
	cur_bin = *(volatile uint32_t *)ptr;
	do {
		old_bin = cur_bin;
		new_bin = old_bin & (~((uint32_t)1 << i));
		cur_bin = atomic_val_CAS32(ptr, old_bin, new_bin);
	} while (cur_bin != old_bin);
	return cur_bin;
}

static inline uint8_t atomic_clear_bit8(uint8_t *ptr, int i)
{
	uint8_t old_bin, new_bin, cur_bin;
	cur_bin = *(volatile uint8_t *)ptr;
	do {
		old_bin = cur_bin;
		new_bin = old_bin & (~((uint8_t)1 << i));
		cur_bin = atomic_val_CAS8(ptr, old_bin, new_bin);
	} while (cur_bin != old_bin);
	return cur_bin;
}

static inline uint32_t atomic_set_bit_val32(uint32_t *ptr, int i, uint32_t b)
{
	uint32_t old_bin, new_bin, cur_bin;
	cur_bin = *(volatile uint32_t *)ptr;
	do {
		old_bin = cur_bin;
		new_bin = old_bin | (b << i);
		cur_bin = atomic_val_CAS32(ptr, old_bin, new_bin);
	} while (cur_bin != old_bin);
	return cur_bin;
}

static inline uint8_t atomic_set_bit_val8(uint8_t *ptr, int i, uint8_t b)
{
	uint8_t old_bin, new_bin, cur_bin;
	cur_bin = *(volatile uint8_t *)ptr;
	do {
		old_bin = cur_bin;
		new_bin = old_bin | (b << i);
		cur_bin = atomic_val_CAS8(ptr, old_bin, new_bin);
	} while (cur_bin != old_bin);
	return cur_bin;
}

static inline uint32_t atomic_set_bit_var32(uint32_t *ptr, int i, uint32_t b)
{
	uint32_t old_bin, new_bin, cur_bin;
	cur_bin = *(volatile uint32_t *)ptr;
	do {
		old_bin = cur_bin;
		new_bin = (old_bin & (~((uint32_t)1 << i))) | (b << i);
		cur_bin = atomic_val_CAS32(ptr, old_bin, new_bin);
	} while (cur_bin != old_bin);
	return cur_bin;
}

static inline uint8_t atomic_set_bit_var8(uint8_t *ptr, int i, uint8_t b)
{
	uint8_t old_bin, new_bin, cur_bin;
	cur_bin = *(volatile uint8_t *)ptr;
	do {
		old_bin = cur_bin;
		new_bin = (old_bin & (~((uint8_t)1 << i))) | (b << i);
		cur_bin = atomic_val_CAS8(ptr, old_bin, new_bin);
	} while (cur_bin != old_bin);
	return cur_bin;
}

#endif /* __ATOMIC_H__ */
