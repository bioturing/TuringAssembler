#ifndef __KMER_HASH__
#define __KMER_HASH__
#include "read_list.h"
#include "helper.h"
#include "khash.h"
KHASH_MAP_INIT_INT64(kmer_int, int);

khash_t(kmer_int) *get_kmer_hash(char *r1_path, char *r2_path, int ksize);
void count_hash_from_read(char *seq, int ksize, khash_t(kmer_int) *h);
int count_kmer_on_seq(khash_t(kmer_int) *h, char *seq, int ksize);
void print_kmer_count_on_seq(khash_t(kmer_int) *h, char *seq, int ksize);
int kmer_check(char *first, char *second, int overlap_ksize, int check_ksize,
		khash_t(kmer_int) *kmer_count);
int count_zero_kmer_map(char *first, char *second, int overlap_ksize,
		int check_ksize, khash_t(kmer_int) *kmer_count);
int count_max_consecutive_zero_kmer(char *first, char *second, int overlap_ksize,
		int check_ksize, khash_t(kmer_int) *kmer_count);
#endif
