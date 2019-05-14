#ifndef SCAFFOLD_EDGE_H
#define SCAFFOLD_EDGE_H
#include "assembly_graph.h"

struct contig_edge {
	int src, des;
	int rv_src, rv_des;
	float score0;
};

struct matrix_score *get_score_edges_matrix(struct asm_graph_t *g, int i0, int i1, int n_bucks, float avg_bin_hash);
int get_score_big_small(int i0, int i1, struct asm_graph_t *g, float avg_bin_hash);
int check_replicate_contig_edge(struct asm_graph_t *g, int i0, int i1, const int n_bucks, float threshold, float avg_bin_hash);
struct bucks_score get_score_edges_res(int i0, int i1, struct asm_graph_t *g, const int n_bucks, float avg_bin_hash) ;
void add_contig_edge(struct asm_graph_t *g,struct contig_edge *listE, int pos, int src, int des, float score0);
void unique_edge(struct contig_edge *listE, int *n_e);
void build_V_from_E(struct contig_edge *listE, int n_e, int **listV, int *n_v);
uint32_t equal_contig_edge(struct contig_edge *e0,struct contig_edge *e1);
uint32_t better_contig_edge(struct contig_edge *e0, struct contig_edge *e1);
int less_contig_edge(const void *e0,const void *e1);
#endif

