#include "khash_operations.h"
#include "log.h"

void put_in_set(khash_t(set_int) *h, int k)
{
	int ret;
	kh_put(set_int, h, k, &ret);
}

void erase_from_set(khash_t(set_int) *h, int k)
{
	khiter_t it = kh_get(set_int, h, k);
	if (it == kh_end(h))
		log_error("Key is not in set, something went wrong!");
	kh_del(set_int, h, it);
}

int check_in_set(khash_t(set_int) *h, int k)
{
	return kh_get(set_int, h, k) != kh_end(h);
}

void put_in_map(khash_t(int_int) *h, int k, int v)
{
	int ret;
	khiter_t it = kh_put(int_int, h, k, &ret);
	kh_val(h, it) = v;
}

int check_in_map(khash_t(int_int) *h, int k)
{
	return kh_get(int_int, h, k) != kh_end(h);
}

void increase_in_map(khash_t(int_int) *h, int k, int v)
{
	khiter_t it = kh_get(int_int, h, k);
	kh_val(h, it) += v;
}

int get_in_map(khash_t(int_int) *h, int k)
{
	khiter_t it = kh_get(int_int, h, k);
	return kh_val(h, it);
}


