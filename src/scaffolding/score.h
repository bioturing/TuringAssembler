#ifndef __SCAFFOLDING_SCORE_H__
#define __SCAFFOLDING_SCORE_H__
#include "scaffolding/edge.h"

struct bucks_score {
	float score;
};

struct matrix_buck_score{
	int n_bucks;
	float *A;
};

struct edges_score_type{
	int n_edge;
	struct scaffold_edge *list_edge;
};

void find_edge_from(struct edges_score_type *edges_score, int i_contig, int *n_edge_adj, 
		struct scaffold_edge **list_edge);

void sort_edges_score(struct edges_score_type *edges_score);

void destroy_matrix_score(struct matrix_buck_score *x);
int detect_anomal_diagonal(struct matrix_buck_score *score, float threshold);
void append_edge_score(struct edges_score_type *edges_score, struct scaffold_edge *edge);
void print_edge_score(struct edges_score_type *edges_score);
#endif /* __SCAFFOLDING_SCORE_H__ */

