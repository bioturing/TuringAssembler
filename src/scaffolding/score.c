#include "scaffolding/score.h"
#include "scaffolding/edge.h"
#include "scaffolding/compare.h"
#include <stdlib.h>
#include <assert.h>
#include "verbose.h"
#include "scaffolding/global_params.h"
void destroy_matrix_score(struct matrix_buck_score *x)
{
	free(x->A);
	free(x);
}

int detect_anomal_diagonal(struct matrix_buck_score *score, float threshold)
{
	int n_bucks = score->n_bucks;
	float sum_dia, sum_all, avg_dia, avg_all;
	sum_dia = 0;
	sum_all = 0;
	avg_dia = 0;
	avg_all = 0;
	int count_dia = 0, count_all = 0;
	for (int i = 0; i < n_bucks; i++) {
		for (int j = 0; j < n_bucks; j++) {
			if (i == j) {
				float sc = score->A[i*n_bucks+j];
				if (sc > -0.000001){
					sum_dia += sc;
					count_dia++;
				}
			} else {
				float sc = score->A[i*n_bucks+j];
				if (sc > -0.000001){
					sum_all += sc;
					count_all++;
				}
			}
		}
	}
	avg_all = sum_all / count_all;
	avg_dia = sum_dia / count_dia;
	return (sum_all + sum_dia > threshold && avg_dia > avg_all * 2);
}

void append_edge_score(struct edges_score_type *edges_score, struct scaffold_edge *edge)
{
	edges_score->n_edge++;
	edges_score->list_edge = realloc(edges_score->list_edge, edges_score->n_edge * 
			sizeof(struct scaffold_edge));
	edges_score->list_edge[edges_score->n_edge-1] = *edge;
}

struct scaffold_edge *find_lower_bound_from(struct edges_score_type *edges_score,
		int i_contig)
{
	int l = 0, r = edges_score->n_edge - 1;
	while (l != r) {
		int mid = (l + r+1)/2;
		if (edges_score->list_edge[mid].src > i_contig) 
			r = mid-1;
		else
			l = mid;
	}
	return &edges_score->list_edge[l];
}

struct scaffold_edge *find_upper_bound_from(struct edges_score_type *edges_score,
		int i_contig)
{
	int l = 0, r = edges_score->n_edge - 1;
	while (l != r) {
		int mid = (l + r)/2;
		if (edges_score->list_edge[mid].src > i_contig) 
			r = mid;
		else
			l = mid + 1;
	}
	return &edges_score->list_edge[l];
}

void find_edge_from(struct edges_score_type *edges_score, int i_contig, int *n_edge_adj, 
		struct scaffold_edge **list_edge)
{
	struct scaffold_edge *start_pos = find_lower_bound_from(edges_score, i_contig);
	struct scaffold_edge *end_pos = find_upper_bound_from(edges_score, i_contig);
	*n_edge_adj = end_pos - start_pos;
	*list_edge = start_pos;
}

void sort_edges_score(struct edges_score_type *edges_score)
{
	qsort(edges_score->list_edge, edges_score->n_edge, sizeof(struct scaffold_edge), ascending_edge); 
}

void print_edge_score(struct edges_score_type *edges_score) 
{
	for (int i = 0; i < edges_score->n_edge; i++) {
		struct scaffold_edge *edge = &edges_score->list_edge[i];
		VERBOSE_FLAG(0, "edges score %d %d %f\n", edge->src, edge->des, edge->score0);
	}
}

