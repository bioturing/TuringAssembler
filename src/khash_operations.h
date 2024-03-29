#ifndef __KHASH_OPERATIONS__
#define __KHASH_OPERATIONS__
#include <assert.h>
#include "khash.h"

#define KHASH_MAP_OPERATIONS(name, key_type, val_type)\
\
static int kh_##name##_exist(khash_t(name) *h, key_type k)\
{\
	khiter_t it = kh_get(name, h, k);\
	return it != kh_end(h);\
}\
\
static val_type kh_##name##_get(khash_t(name) *h, key_type k)\
{\
	khiter_t it = kh_get(name, h, k);\
	assert(it != kh_end(h));\
	return kh_val(h, it);\
}\
static val_type kh_##name##_try_get(khash_t(name) *h, key_type k, val_type dump)\
{\
	khiter_t it = kh_get(name, h, k);\
	if (it == kh_end(h))\
		return dump;\
	return kh_val(h, it);\
}\
static void kh_##name##_set(khash_t(name) *h, key_type k, val_type v)\
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
static int kh_##name##_exist(khash_t(name) *h, key_type k)\
{\
	khiter_t it = kh_get(name, h, k);\
	return it != kh_end(h);\
}\
static void kh_##name##_insert(khash_t(name) *h, key_type k)\
{\
	khiter_t it = kh_get(name, h, k);\
	if (it == kh_end(h)){\
		int ret;\
		it = kh_put(name, h, k, &ret);\
	}\
}\
static void kh_##name##_erase(khash_t(name) *h, key_type k)\
{\
	khiter_t it = kh_get(name, h, k);\
	assert(it != kh_end(h));\
	kh_del(name, h, it);\
}


KHASH_SET_INIT_INT(set_int);
KHASH_SET_OPERATIONS(set_int, int);

KHASH_MAP_INIT_INT(int_int, int);
KHASH_MAP_OPERATIONS(int_int, int, int);

void put_in_set(khash_t(set_int) *h, int k);
void erase_from_set(khash_t(set_int) *h, int k);
int check_in_set(khash_t(set_int) *h, int k);
void put_in_map(khash_t(int_int) *h, int k, int v);
void increase_in_map(khash_t(int_int) *h, int k, int v);
int check_in_map(khash_t(int_int) *h, int k);
int get_in_map(khash_t(int_int) *h, int k);
#endif
