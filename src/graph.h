#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "attribute.h"

struct graph_t {
	// Node is 1-based but store as 0-based
	int n_v, n_k, n_e;

	// kmer list info
	int *kmer_count;
	int *chain_head;
	uint64_t *chain_kmer;

	int *kmer_chain_id;

	// adjacency list
	int *fhead, *rhead;
	int *fadj, *radj;
};


int16_t *get_edges(struct opt_count_t *opt, khash_t(kvert) *h);

void reduce_graph(struct opt_count_t *opt, khash_t(kvert) *h, int16_t *e);

#endif
