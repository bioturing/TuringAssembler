#ifndef SCAFFOLD_EDGE_H
#define SCAFFOLD_EDGE_H
#include "assembly_graph.h"

struct scaffold_edge{
	int src, des;
	int rv_src, rv_des;
	float score0;
};

struct candidate_edge {
	int src, des;
	float score;
};

struct scaffold_edge *new_scaffold_edge(int src, int des, float score);
struct matrix_buck_score *get_score_edges_matrix(struct asm_graph_t *g, int i0, int i1, int n_bucks,
        float avg_bin_hash, struct opt_proc_t *opt);
int get_score_big_small(int i0, int i1, struct asm_graph_t *g, float avg_bin_hash);
int check_replicate_scaffold_edge(struct asm_graph_t *g, int i0, int i1, 
		const int n_bucks, float threshold, float avg_bin_hash, struct opt_proc_t *opt);
struct bucks_score get_score_edges_res(int i0, int i1, struct asm_graph_t *g, const int n_bucks, 
		float avg_bin_hash, struct opt_proc_t *opt);
void add_scaffold_edge(struct asm_graph_t *g,struct scaffold_edge *listE, int pos, int src, int des, 
		float score0, struct opt_proc_t *opt);
void unique_edge(struct scaffold_edge *listE, int *n_e);
void build_V_from_E(struct scaffold_edge *listE, int n_e, int **listV, int *n_v);
uint32_t equal_scaffold_edge(struct scaffold_edge *e0,struct scaffold_edge *e1);
uint32_t better_scaffold_edge(struct scaffold_edge *e0, struct scaffold_edge *e1);
int less_scaffold_edge(const void *e0,const void *e1);
#endif

