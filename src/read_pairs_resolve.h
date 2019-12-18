#ifndef __READ_PAIRS_RESOLVE__
#define __READ_PAIRS_RESOVLE__
#include <stdio.h>
#include <stdlib.h>
#include "assembly_graph.h"
#include "khash.h"
#include "khash_operations.h"
#include "log.h"
#define MIN_READ_PAIR_MAPPED 25
KHASH_MAP_INIT_INT64(long_int, int);
KHASH_MAP_OPERATIONS(long_int, uint64_t, int);

struct read_pair_cand_t{
	int n;
	int *cand;
};

void get_read_pairs_count(struct asm_graph_t *g, char *path,
		struct read_pair_cand_t *rp_cand);

void merge_sort_by_edge_length(int *edges, int l, int r);

void get_long_contigs(struct opt_proc_t *opt);
#endif
