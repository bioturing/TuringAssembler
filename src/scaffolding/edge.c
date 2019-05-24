#include "assembly_graph.h"
#include <stdlib.h>
#include "verbose.h"
#include "math.h"
#include "scaffolding/global_params.h"
#include "scaffolding/buck.h"
#include "scaffolding/edge.h"
#include "scaffolding/contig.h"
#include "scaffolding/compare.h"
#include "scaffolding/score.h"

struct matrix_score *get_score_edges_matrix(struct asm_graph_t *g, int i0, int i1, int n_bucks,
 			float avg_bin_hash, struct opt_proc_t *opt)
{
	struct matrix_score *score = NULL;
	score = realloc(score, sizeof(struct matrix_score));
	int rev_i0 = g->edges[i0].rc_id;
	assert(rev_i0 < g->n_e);
	struct asm_edge_t *rev_e0 = &g->edges[rev_i0], *e1 = &g->edges[i1];
	int n0_bucks = (get_edge_len(rev_e0) + g->bin_size-1) / g->bin_size;
	int n1_bucks = (get_edge_len(e1) + g->bin_size-1) / g->bin_size;
	int e1_len = get_edge_len(e1), e0_len = get_edge_len(rev_e0);
	VERBOSE_FLAG(3, "len e0, e1: %d %d \n" , e0_len, e1_len);
	assert(n0_bucks > n_bucks);
	assert(n1_bucks > n_bucks);

	score->n_bucks = n_bucks;
	score->A = NULL;
	score->A = realloc(score->A, n_bucks * n_bucks * sizeof(float));
	// check_bucks_A[i] && check_bucks_B[i] de improve performance
	VERBOSE_FLAG(2, "cov %f binsize %d ksize %d", global_genome_coverage, g->bin_size, g->ksize);
	int *check_bucks_A = NULL, *check_bucks_B = NULL;
	check_bucks_A = realloc(check_bucks_A, n_bucks * sizeof(int));
	check_bucks_B = realloc(check_bucks_B, n_bucks * sizeof(int));
	for (int i = 0; i < score->n_bucks; ++i) {
		check_bucks_A[i] = check_qualify_buck(g, rev_e0, i, avg_bin_hash, opt);
		check_bucks_B[i] = check_qualify_buck(g, e1, i, avg_bin_hash, opt);
	}
	float cov_rev_e0 = __get_edge_cov(rev_e0, g->ksize);
	float cov_e1 = __get_edge_cov(e1, g->ksize);
	for (int i = 0; i < n_bucks; ++i) {
		for (int j = 0; j < n_bucks; ++j) {
			if (check_bucks_A[i] && check_bucks_B[j]) {
				score->A[i*n_bucks+j] = get_score_bucks(&rev_e0->bucks[i],
						&e1->bucks[j], cov_rev_e0, cov_e1);
			}
			else 
				score->A[i*n_bucks+j] = -1;
		}
	}

	free(check_bucks_A);
	free(check_bucks_B);
	return score;
}

int get_score_big_small(int i0, int i1, struct asm_graph_t *g, float avg_bin_hash) 
{
	struct asm_edge_t *e0 = &g->edges[i0];
	struct asm_edge_t *e1 = &g->edges[i1];
	int n_bucks = global_thres_n_buck_big_small;
	int n0_bucks = (get_edge_len(e0) + g->bin_size-1) / g->bin_size;
	int n1_bucks = (get_edge_len(e1) + g->bin_size-1) / g->bin_size;
	assert(n0_bucks >= n1_bucks);
	int score;
	score = 0;

	assert(n0_bucks >= n_bucks && n1_bucks >= n_bucks);
	float res = 0;
	float maxtmp = 0;
	float cov_e0 = __get_edge_cov(e0, g->ksize);
	float cov_e1 = __get_edge_cov(e1, g->ksize);
	for (int i = MAX(0, n0_bucks-20); i < n0_bucks-1; ++i) {
		float left_value = 0, right_value = 0, count_left = 0, count_right = 0;
		for (int j = 0; j < 3; ++j) {
			float tmp = 0;
			tmp = get_score_bucks(&e0->bucks[i],& e1->bucks[j], cov_e0, cov_e1);
			if (tmp > global_thres_bucks_score) {
				left_value += tmp;
				count_left++;
			}
		}
		for (int j = n1_bucks - 4; j < n1_bucks-1; ++j) {
			float tmp = 0;
			tmp = get_score_bucks(&e0->bucks[i], &e1->bucks[j], cov_e0, cov_e1);
			if (tmp > global_thres_bucks_score) {
				right_value += tmp;
				count_right++;
			}
			VERBOSE_FLAG(3, "%f ", tmp);
		}
		left_value /= count_left;
		right_value /= count_right;
		VERBOSE_FLAG(3, "%f %f", left_value , right_value);
		if (fabsf(left_value - right_value) > (0.05 * MAX(left_value, right_value))) {
			VERBOSE_FLAG(2, "%f", fabsf(left_value - right_value));
			if (left_value > right_value)
				score++;
			else 
				score--;
		}
	}
	return score;
}

int check_replicate_contig_edge(struct asm_graph_t *g, int i0, int i1, 
		const int n_bucks, float threshold, float avg_bin_hash, struct opt_proc_t *opt)
{
	int rev_i0 = g->edges[i0].rc_id, rev_i1 = g->edges[i1].rc_id;
	struct matrix_score *s0 = get_score_edges_matrix(g, i0, i1, n_bucks, avg_bin_hash, opt);
	struct matrix_score *s1 = get_score_edges_matrix(g, i0, rev_i1, n_bucks, avg_bin_hash, opt);
	struct matrix_score *s2 = get_score_edges_matrix(g, rev_i0, i1, n_bucks, avg_bin_hash, opt);
	struct matrix_score *s3 = get_score_edges_matrix(g, rev_i0, rev_i1, n_bucks, avg_bin_hash, opt);
	int res = (detect_anomal_diagonal(s0, threshold) || detect_anomal_diagonal(s1, threshold) || detect_anomal_diagonal(s2, threshold) || detect_anomal_diagonal(s3, threshold));
	destroy_matrix_score(s0);
	destroy_matrix_score(s1);
	destroy_matrix_score(s2);
	destroy_matrix_score(s3);
	return res;
}

struct bucks_score get_score_edges_res(int i0, int i1, struct asm_graph_t *g, const int n_bucks, 
		float avg_bin_hash, struct opt_proc_t *opt) 
{
	struct matrix_score *mat_score = get_score_edges_matrix(g, i0, i1, n_bucks, avg_bin_hash, opt);
	mat_score->A[0] = -1;
	struct asm_edge_t *e0 = &g->edges[i0], *e1 = &g->edges[i1];
	float res = 0;
	int count = 0;
	for (int i = 0; i < mat_score->n_bucks; ++i) {
		for (int j = 0; j < mat_score->n_bucks; ++j) {
			float tmp = mat_score->A[i * n_bucks + j];
			VERBOSE_FLAG(3, "%f ", tmp);
			if (tmp >= -0.000001) {
				count++;
				res += tmp;
			}
		}
		VERBOSE_FLAG(3, "#\n");
	}
	struct bucks_score res_score;
	res_score.score = res/count;
	destroy_matrix_score(mat_score);
	return res_score;
}

void add_contig_edge(struct asm_graph_t *g,struct contig_edge *listE, int pos, int src, int des, 
		float score0, struct opt_proc_t *opt)
{
	assert(src >= 0 && des >= 0 && src < g->n_e && des < g->n_e);
	struct contig_edge *e = calloc(1, sizeof(struct contig_edge));
	e->src = src;
	e->des = des;
	e->score0 = score0;
	e->rv_src = 0;
	e->rv_des = 0;
	normalize_min_index(g, e);
	listE[pos] = *e;
}

void unique_edge(struct contig_edge *listE, int *n_e)
{
	int new_n_e = 0;
	for (int i = 0; i < *n_e; ) {
		int j = i;
		struct contig_edge *best = calloc(1, sizeof(struct contig_edge));
		best->score0 = -1;
		while (j < *n_e && equal_contig_edge(&listE[j], &listE[i])){
			if (better_contig_edge(&listE[j], best)){
				best = &listE[j];
			}
			++j;
		}
		listE[new_n_e++] = *best;
		i = j;
	}
	*n_e = new_n_e;
}

void build_V_from_E(struct contig_edge *listE, int n_e, int **listV, int *n_v)
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

uint32_t equal_contig_edge(struct contig_edge *e0,struct contig_edge *e1)
{
	return (e0->src == e1->src && e0->des == e1->des);
}

uint32_t better_contig_edge(struct contig_edge *e0, struct contig_edge *e1)
{
	return (e0->score0 > e1->score0);
}

int less_contig_edge(const void *e0,const void *e1)
{
	struct contig_edge * v0 = (struct contig_edge *) e0;
	struct contig_edge * v1 = (struct contig_edge *) e1;
	return (v0->src > v1->src || (v0->src == v1->src && v0->des > v1->des));
}

