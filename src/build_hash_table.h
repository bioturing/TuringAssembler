//
// Created by che on 11/11/2019.
//

#ifndef SKIPPING_BUILD_HASH_TABLE_H
#define SKIPPING_BUILD_HASH_TABLE_H

#include "khash.h"
KHASH_MAP_INIT_INT(pair_kmer_count, int);

struct kmer_pair_iterator_bundle_t {
	struct dqueue_t *q;
	char prefix[MAX_PATH];
	int64_t sm;
};
void build_pair_kmer_table(struct opt_proc_t *opt);
void get_seq(char *seq, int start, int len, uint8_t *res);
void copy_seq32_seq8(uint32_t *seq, int start, uint8_t *res, int start_res ,int len);
void print_u8_seq(uint8_t *a, int len);
uint8_t *compress_seq(char *a);

#endif //SKIPPING_BUILD_HASH_TABLE_H
