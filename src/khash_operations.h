#ifndef __KHASH_OPERATIONS__
#define __KHASH_OPERATIONS__
#include "khash.h"

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
