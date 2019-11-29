#ifndef __KHASH_OPERATIONS__
#define __KHASH_OPERATIONS__
#include <assert.h>
#include "khash.h"

#define KHASH_MAP_OPERATIONS(name, key_type, val_type)\
\
static int kh_##name##_key_exist(khash_t(name) *h, key_type k)\
{\
	khiter_t it = kh_get(name, h, k);\
	return it != kh_end(h);\
}\
\
static val_type kh_##name##_get_val(khash_t(name) *h, key_type k)\
{\
	khiter_t it = kh_get(name, h, k);\
	assert(it != kh_end(h));\
	return kh_val(h, it);\
}\
static void kh_##name##_set_val(khash_t(name) *h, key_type k, val_type v)\
{\
	khiter_t it = kh_get(name, h, k);\
	if (it == kh_end(h)){\
		int ret;\
		it = kh_put(name, h, k, &ret);\
	}\
	kh_val(h, it) = v;\
}

#define KHASH_SET_OPERATIONS(name, key_type)\
\
static int kh_##name##_key_exist(khash_t(name) *h, key_type k)\
{\
	khiter_t it = kh_get(name, h, k);\
	return it != kh_end(h);\
}\
\
static void kh_##name##_add_key(khash_t(name) *h, key_type k)\
{\
	khiter_t it = kh_get(name, h, k);\
	if (it == kh_end(h)){\
		int ret;\
		it = kh_put(name, h, k, &ret);\
	}\
}

KHASH_SET_INIT_INT(set_int);
KHASH_MAP_INIT_INT(int_int, int);

void put_in_set(khash_t(set_int) *h, int k);
void erase_from_set(khash_t(set_int) *h, int k);
int check_in_set(khash_t(set_int) *h, int k);
void put_in_map(khash_t(int_int) *h, int k, int v);
void increase_in_map(khash_t(int_int) *h, int k, int v);
int check_in_map(khash_t(int_int) *h, int k);
int get_in_map(khash_t(int_int) *h, int k);
#endif
