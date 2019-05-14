#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(_MSC_VER)
#include <time.h>
#include <windows.h>
#include <getopt.h>
#include <BaseTsd.h>
#else
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#endif /* _MSC_VER */

#include "utils.h"

/* no use of verbse/log/error here, using stderr if necessary */

// __KHASH_IMPL(mer26, kh_inline, khint64_t, void *, 1, kh_int64_hash_func, kh_int64_hash_equal)
// __KHASH_IMPL(int64, kh_inline, khint64_t, khint64_t, 1, kh_int64_hash_func, kh_int64_hash_equal)

int8_t nt4_table[256] = {
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4
};

char *nt4_char = "ACGTN";
char *rev_nt4_char = "TGCAN";

char *str_concate(const char *str1, const char *str2)
{
	/* TODO: put this to depricated */
	size_t len1 = strlen(str1);
	size_t len2 = strlen(str2);
	char *str3 = malloc(len1 + len2 + 1);
	strcpy(str3, str1);
	strcpy(str3 + len1, str2);
	str3[len1 + len2] = '\0';
	return str3;
}

char *get_rev(const char *seq, int len)
{
	if (seq == NULL)
		return NULL;

	int i, k;
	char *ret = malloc(len + 1);
	for (i = 0, k = len - 1; i < len; ++i, --k)
		ret[i] = seq[k];
	ret[len] = '\0';
	return ret;
}

char *get_rev_complement(const char *seq, int len)
{
	if (seq == NULL)
		return NULL;

	char *ret = malloc(len + 1);
	int i, k;
	for (i = 0, k = len - 1; i < len; ++i, --k)
		ret[i] = rev_nt4_char[nt4_table[(int)seq[k]]];
	ret[len] = '\0';
	return ret;
}

double realtime()
{
	struct timeval tp;
#if defined(_MSC_VER)
	_gettimeofday(&tp,  NULL);
#else
	gettimeofday(&tp,  NULL);
#endif
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

int64_t seq2num(const char *seq, int len)
{
	int64_t ret = 0;
	int i;
	for (i = 0; i < len; ++i)
		ret = ret * 5 + nt4_table[(int)seq[i]];
	return ret;
}

char *num2seq(int64_t num, int len)
{
	char *ret = malloc(len + 1);
	int i;
	for (i = 0; i < len; ++i) {
		ret[len - i - 1] = nt4_char[num % 5];
		num /= 5;
	}
	ret[len] = '\0';
	return ret;
}

int *binary_search(int *list, int n, int value)
{
	int l = 0, r = n - 1;
	while (l != r) {
		int m = (l + r) / 2;
		if (value > list[m]) 
			l = m + 1;
		else
			r = m;
	}
	assert(list[l] == value);
	if (list[l] != value)
		return list-1;
	return &list[l];
}

void unique(int *listV, int *n_v)
{
	int new_n_v = 0;
	for (int i = 0; i < *n_v; ) {
		int j = i;
		while (j < *n_v && listV[i] == listV[j]) j++;
		listV[new_n_v++] = listV[i];
		i = j;
	}
	*n_v = new_n_v;
}

