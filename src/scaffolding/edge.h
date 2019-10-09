#ifndef SCAFFOLD_EDGE_H
#define SCAFFOLD_EDGE_H
#include "assembly_graph.h"
#include "scaffolding/score.h"

struct edges_score_type{
	int n_edge;
	struct scaffold_edge *list_edge;
};
void find_edge_from(struct edges_score_type *edges_score, int i_contig, int *n_edge_adj, 
		struct scaffold_edge **list_edge);

void sort_edges_score(struct edges_score_type *edges_score);

void append_edge_score(struct edges_score_type *edges_score, struct scaffold_edge *edge);
void print_edge_score(struct edges_score_type *edges_score);

struct scaffold_edge{
	int src, des;
	struct pair_contigs_score score;
};

struct scaffold_edge *new_scaffold_edge(int src, int des, struct pair_contigs_score *score);
struct matrix_buck_score *get_score_edges_matrix(struct asm_graph_t *g, int i0, int i1, int n_bucks,
        float avg_bin_hash, struct opt_proc_t *opt);
int get_score_big_small(int i0, int i1, struct asm_graph_t *g, float avg_bin_hash);
int check_replicate_scaffold_edge(struct asm_graph_t *g, int i0, int i1, 
		const int n_bucks, float threshold, float avg_bin_hash, struct opt_proc_t *opt);
struct pair_contigs_score *get_score_edges_res(int i0, int i1, struct asm_graph_t *g, 
		float avg_bin_hash, struct opt_proc_t *opt);
void add_scaffold_edge(struct asm_graph_t *g,struct scaffold_edge *listE, int pos, int src, int des, 
		float score0, struct opt_proc_t *opt);
void unique_edge(struct scaffold_edge *listE, int *n_e);
void build_V_from_E(struct scaffold_edge *listE, int n_e, int **listV, int *n_v);
uint32_t equal_scaffold_edge(struct scaffold_edge *e0,struct scaffold_edge *e1);
uint32_t better_scaffold_edge(struct scaffold_edge *e0, struct scaffold_edge *e1);
int less_scaffold_edge(const void *e0,const void *e1);
void destroy_edges_score_type(struct edges_score_type *edges_score);
#endif

