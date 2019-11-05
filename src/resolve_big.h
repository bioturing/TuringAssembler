//
// Created by che on 31/10/2019.
//

#ifndef SKIPPING_RESOLVE_BIG_H
#define SKIPPING_RESOLVE_BIG_H


KHASH_MAP_INIT_INT64(big_kmer_count, int);
KHASH_SET_INIT_INT64(union_barcode);

int resolve_using_big_kmer(struct asm_graph_t *g, int i_e, const int *partition,
						   khash_t(big_kmer_count) **partititon_kmer_count);

void partition_graph(struct read_path_t *ori_read, struct asm_graph_t *g, int *partition, int n_threads,
					 int mmem, khash_t(big_kmer_count) ***res_partition_kmer_count, int *n_partition);

struct km_count_bundle_t {
	khash_t(big_kmer_count) **kmer_count_table;
	int ksize;
};
void km_count_bundle_destroy(struct km_count_bundle_t *b);
#endif //SKIPPING_RESOLVE_BIG_H
