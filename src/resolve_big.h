//
// Created by che on 31/10/2019.
//

#ifndef SKIPPING_RESOLVE_BIG_H
#define SKIPPING_RESOLVE_BIG_H

#include "build_hash_table.h"
#include <stdio.h>


KHASH_MAP_INIT_INT64(big_kmer_count, int);
KHASH_SET_INIT_INT64(union_barcode);


struct km_count_bundle_t {
	khash_t(big_kmer_count) *kmer_count_table;
	int ksize;
	pthread_mutex_t lock;

};

struct partition_information {
	int total_len;
	int n_barcodes;
	int *map_node;
	khash_t(union_barcode) *union_barcodes;
	struct asm_graph_t *sub_graph;
	FILE *output_gfa_file;
};
void km_count_bundle_destroy(struct km_count_bundle_t *b);
void resolve_1_2(struct asm_graph_t *g, struct opt_proc_t *opt);
#endif //SKIPPING_RESOLVE_BIG_H
