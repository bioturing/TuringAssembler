#ifndef __SPLIT_MOLECULES__
#define __SPLIT_MOLECULES__
#include "assembly_graph.h"
#include "complex_resolve.h"
#include "khash.h"

struct line_vertex_t{
	int deg_in;
	int *parents;
	int deg_out;
	int *children;
};
KHASH_MAP_INIT_INT(edge_line, struct line_vertex_t *);

struct line_graph_t{
	int n_v;
	int *line_v;
	khash_t(edge_line) *vertices;
};

void init_line_graph(struct line_graph_t *lig, int n_e, int *edges);
void construct_line_graph(struct asm_graph_t *g, struct line_graph_t *lig);
void get_edges_in_radius(struct asm_graph_t *g, int e, khash_t(set_int) *nearby);
void get_edges_in_radius_dfs(struct asm_graph_t *g, int e, int len,
		khash_t(set_int) *visited);
void add_line_edge(struct line_graph_t *lig, int v, int u);
//struct edge_ordering_t{
//	int *edges;
//	int n;
//};
//
//void get_edges_order(struct asm_graph_t *g, int *edges, int n_e,
//		struct edge_ordering_t *order);
//int get_edges_order_dfs(struct asm_graph_t *g, int e, int p, int total, int *path,
//		int depth, khash_t(int_int) *multi, struct edge_ordering_t *order);
#endif
