#include "assembly_graph.h"
#include <stdlib.h>
#include "verbose.h"
#include <math.h>
#include <assert.h>
#include "scaffolding/global_params.h"
#include "scaffolding/buck.h"
#include "scaffolding/edge.h"
#include "scaffolding/contig.h"
#include "scaffolding/compare.h"
#include "scaffolding/score.h"
#include "utils.h"

struct pair_contigs_score *get_score_l_l_mat(struct asm_graph_t *g, int i0, int i1,
 			float avg_bin_hash, struct opt_proc_t *opt)
{
	int rev_i0 = g->edges[i0].rc_id;
	assert(rev_i0 < g->n_e);
	struct asm_edge_t *rev_e0 = &g->edges[rev_i0], *e1 = &g->edges[i1];
	int e1_len = get_edge_len(e1), e0_len = get_edge_len(rev_e0);
	VERBOSE_FLAG(3, "len e0, e1: %d %d \n" , e0_len, e1_len);

	float cov_rev_e0 = __get_edge_cov(rev_e0, g->ksize);
	float cov_e1 = __get_edge_cov(e1, g->ksize);

	struct pair_contigs_score *score = calloc(1, sizeof(struct pair_contigs_score));
	score->bc_score = get_share_barcode(&rev_e0->barcodes, &e1->barcodes, cov_rev_e0, cov_e1, avg_bin_hash);
	//todo @huu not hardcode
//	if (score->bc_score < 0.3){
//		return score;
//	}
	score->m_score = get_share_mate(g, i0, i1);
	score->m2_score = get_share_mate_2(g, i0, i1);
	return score;
}

int get_score_big_small(int i0, int i1, struct asm_graph_t *g, float avg_bin_hash) 
{
	//todo @huu
	assert(0);
}

int check_replicate_scaffold_edge(struct asm_graph_t *g, int i0, int i1, 
		const int n_bucks, float threshold, float avg_bin_hash, struct opt_proc_t *opt)
{
	return 0;
//	int rev_i0 = g->edges[i0].rc_id, rev_i1 = g->edges[i1].rc_id;
//	float s0 = get_score_l_l_mat(g, i0, i1, n_bucks, avg_bin_hash, opt);
//	float s1 = get_score_l_l_mat(g, i0, rev_i1, n_bucks, avg_bin_hash, opt);
//	float s2 = get_score_l_l_mat(g, rev_i0, i1, n_bucks, avg_bin_hash, opt);
//	float s3 = get_score_l_l_mat(g, rev_i0, rev_i1, n_bucks, avg_bin_hash, opt);
//	return res;
}

struct pair_contigs_score *get_score_edges_res(int i0, int i1, struct asm_graph_t *g, 
		float avg_bin_hash, struct opt_proc_t *opt) 
{
	struct pair_contigs_score *pair_score = calloc(1, sizeof(struct pair_contigs_score));
	pair_score->bc_score = -1;

	struct asm_edge_t *e0 = &g->edges[i0], *e1 = &g->edges[i1];
	if (is_very_short_contig(e0) || is_very_short_contig(e1))
		return pair_score;
	//todo @huu get_score_l_s
	//todo @huu get_score_s_l

	return get_score_l_l_mat(g, i0, i1, avg_bin_hash, opt);
}

void unique_edge(struct scaffold_edge *listE, int *n_e)
{
	int new_n_e = 0;
	for (int i = 0; i < *n_e; ) {
		int j = i;
		struct scaffold_edge *best = calloc(1, sizeof(struct scaffold_edge));
		best->score.bc_score = -1;
		while (j < *n_e && equal_scaffold_edge(&listE[j], &listE[i])){
			if (better_scaffold_edge(&listE[j], best)){
				best = &listE[j];
			}
			++j;
		}
		listE[new_n_e++] = *best;
		i = j;
	}
	*n_e = new_n_e;
}

void build_V_from_E(struct scaffold_edge *listE, int n_e, int **listV, int *n_v)
{
	*listV = NULL; 
	*n_v = 0;
	for (int i = 0; i < n_e; i++) {
		VERBOSE_FLAG(3, "listE edges: %d %d BBBB\n", listE[i].src, listE[i].des);

		++(*n_v);
		int t = (*n_v)*sizeof(int);
		*listV = realloc(*listV, t);
		(*listV)[(*n_v)-1] = listE[i].src;

		++(*n_v);
		*listV = realloc(*listV, (*n_v)*(sizeof(int)));
		(*listV)[(*n_v)-1] = listE[i].des;
	}
	qsort(*listV, *n_v, sizeof(int), ascending_unint32);
	for (int i = 0; i < *n_v; i++) {
		VERBOSE_FLAG(3, "list V: %d ", (*listV)[i]);
	}
	unique(*listV, n_v);
}

uint32_t equal_scaffold_edge(struct scaffold_edge *e0,struct scaffold_edge *e1)
{
	return (e0->src == e1->src && e0->des == e1->des);
}

uint32_t better_scaffold_edge(struct scaffold_edge *e0, struct scaffold_edge *e1)
{
	return (e0->score.bc_score > e1->score.bc_score);
}

int less_scaffold_edge(const void *e0,const void *e1)
{
	struct scaffold_edge * v0 = (struct scaffold_edge *) e0;
	struct scaffold_edge * v1 = (struct scaffold_edge *) e1;
	return (v0->src > v1->src || (v0->src == v1->src && v0->des > v1->des));
}

struct scaffold_edge *new_scaffold_edge(int src, int des, struct pair_contigs_score *score)
{
	struct scaffold_edge* new_edge = calloc(1, sizeof(struct scaffold_edge));
	new_edge->src = src;
	new_edge->des = des;
	new_edge->score = *score;
	return new_edge;
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
		int mid = (l + r)/2;
		if (edges_score->list_edge[mid].src < i_contig) 
			l = mid+1;
		else
			r = mid;
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
		VERBOSE_FLAG(0, "edges score %d %d %f %d m2score %d\n", edge->src, edge->des, edge->score.bc_score, edge->score.m_score, edge->score.m2_score);
	}
}

