#ifndef __READ_PAIRS_RESOLVE__
#define __READ_PAIRS_RESOVLE__
#include <stdio.h>
#include <stdlib.h>
#include "assembly_graph.h"
#include "khash.h"
#include "khash_operations.h"
#include "log.h"
#include "cluster_molecules.h"
#include "verbose.h"
#include "helper.h"
#define MIN_READ_PAIR_MAPPED_HARD 25
#define MIN_READ_PAIR_MAPPED_SOFT 2

struct read_pair_cand_t{
	int n;
	int *cand;
	int *score;
};

void get_read_pairs_count(struct asm_graph_t *g, char *path,
		struct read_pair_cand_t *rp_cand);

void merge_sort_by_edge_length(int *edges, int l, int r);

void extend_by_read_pairs(struct asm_graph_t *g, int s, float unit_cov,
		struct read_pair_cand_t *rp_cand, int **path, int *n_path);

int get_next_cand(struct asm_graph_t *g, float unit_cov, struct read_pair_cand_t *rp_cand,
		int *path, int n_path);

void get_long_contigs(struct opt_proc_t *opt);
#endif
