#ifndef __SPLIT_MOLECULES__
#define __SPLIT_MOLECULES__
#include "assembly_graph.h"
#include "complex_resolve.h"

struct edge_ordering_t{
	int *edges;
	int n;
};

void get_edges_order(struct asm_graph_t *g, int *edges, int n_e,
		struct edge_ordering_t *order);
int get_edges_order_dfs(struct asm_graph_t *g, int e, int p, int total, int *path,
		int depth, khash_t(int_int) *multi, struct edge_ordering_t *order);
#endif
