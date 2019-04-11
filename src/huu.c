#include <stdlib.h>
#include "assembly_graph.h"
#include "k31hash.h"
#include "pthread.h"
#include "verbose.h"
#include <string.h>
#include "algorithm.h"

const int global_thres_length = 10000;
const float global_thres_score = 0.015;
const int global_n_buck = 6;

struct bucks_score {
	float score;
};

uint32_t min(uint32_t a, uint32_t b)
{
	if (a < b) return a;
	else return b;
}

uint32_t roundint(float x)
{ 
	return (x)>=0?(int)((x)+0.5):(int)((x)-0.5);
}

float get_score_bucks(struct barcode_hash_t *buck0,struct barcode_hash_t *buck1) 
{
	const uint32_t thres_cnt = 30;
	uint32_t cnt0 = 0, cnt1 = 0, res2 = 0;
	__VERBOSE("bucksize %d ", buck1->size);
	uint32_t O = 0;
	

	for (uint32_t i = 0; i < buck1->size; ++i) {
		if (buck1->cnts[i] != (uint32_t)(-1) && buck1->cnts[i] >= thres_cnt) {
			cnt1++;
		}
	}

	for (uint32_t i = 0; i < buck0->size; ++i) {
		O = 0;
		if ((buck0->keys[i]) != (uint64_t)(-1) && buck0->cnts[i] >= thres_cnt) {
			cnt0++;
			uint32_t tmp = barcode_hash_get(buck1, buck0->keys[i]);
			if (tmp != BARCODE_HASH_END(buck1) && buck1->cnts[tmp] >= thres_cnt) {
				res2++;
			}
		}
	}
//	__VERBOSE("%d ", res2);
	return 1.0 * res2 / (cnt0 + cnt1);
}

uint32_t abssub(uint32_t a, uint32_t b) {
	if (a>b) 
		return a-b;
	else 
		return b-a;
}

struct matrix_score{
	int n_bucks;
	float *A;
};

struct matrix_score *get_score_edges_matrix(struct asm_graph_t *g, uint32_t i0, uint32_t i1, int n_bucks)
{
	struct matrix_score *score = NULL;
	score = realloc(score, sizeof(struct matrix_score));
	uint32_t rev_i0 = g->edges[i0].rc_id;
	struct asm_edge_t *rev_e0 = &g->edges[rev_i0], *e1 = &g->edges[i1];

	score->n_bucks = n_bucks;
	score->A = NULL;
	score->A = realloc(score->A, n_bucks * n_bucks * sizeof(float));
	for (int i = 0; i < score->n_bucks; ++i) {
		for (int j = 0; j < score->n_bucks; ++j) {
			score->A[i*score->n_bucks+j] = get_score_bucks(&rev_e0->bucks[i], &e1->bucks[j]);
		}
	}

	return score;
}

int detect_anomal_diagonal(struct matrix_score *score, float threshold)
{
	int n = score->n_bucks;
	float sum_dia, sum_all, avg_dia, avg_all;
	sum_dia = 0;
	sum_all = 0;
	avg_dia = 0;
	avg_all = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				sum_dia += score->A[i*n+j];
			} else {
				sum_all += score->A[i*n+j];
			}
		}
	}
	avg_all = sum_all / (n*n-n);
	avg_dia = sum_dia / n;
	return (sum_all + sum_dia > threshold && avg_dia > avg_all * 2);
}

int check_replicate_contig_edge(struct asm_graph_t *g, uint32_t i0, uint32_t i1, const int n_bucks, float threshold)
{
	uint32_t rev_i0 = g->edges[i0].rc_id, rev_i1 = g->edges[i1].rc_id;
	struct matrix_score *s0 = get_score_edges_matrix(g, i0, i1, n_bucks);
	struct matrix_score *s1 = get_score_edges_matrix(g, i0, rev_i1, n_bucks);
	struct matrix_score *s2 = get_score_edges_matrix(g, rev_i0, i1, n_bucks);
	struct matrix_score *s3 = get_score_edges_matrix(g, rev_i0, rev_i1, n_bucks);
	return (detect_anomal_diagonal(s0, threshold) || detect_anomal_diagonal(s1, threshold) || detect_anomal_diagonal(s2, threshold) || detect_anomal_diagonal(s3, threshold));
}

struct bucks_score get_score_edges_res(uint32_t i0, uint32_t i1, struct asm_graph_t *g, const int n_bucks) {
	struct matrix_score *mat_score = get_score_edges_matrix(g, i0, i1, n_bucks);
	struct asm_edge_t *e0 = &g->edges[i0], *e1 = &g->edges[i1];
	uint32_t n0_bucks = (get_edge_len(e0) + g->bin_size-1) / g->bin_size;
	uint32_t n1_bucks = (get_edge_len(e1) + g->bin_size-1) / g->bin_size;
	float res = 0;
	__VERBOSE("edges: %d %d\n", i0, i1);
	for (int i = 0; i < mat_score->n_bucks; ++i) {
		for (int j = 0; j < mat_score->n_bucks; ++j) {
			float tmp = mat_score->A[i * n_bucks + j];
			if (i!=0 || j!=0) {
			} else {
				tmp = 0;
			}
//			__VERBOSE("%f ", tmp);
			tmp /= (mat_score->n_bucks * mat_score->n_bucks-1);
			res += tmp;
		}
//		__VERBOSE("\n");
	}
//	__VERBOSE("\n");
	struct bucks_score res_score;
	res_score.score = res;
//	__VERBOSE("%f\n", res_score.score);
	return res_score;
}


void listContig(struct asm_graph_t *g, FILE *out_file) {
	const uint32_t thres_len_e = global_thres_length; 
	const int n_bucks = global_n_buck;
	float thres_score = global_thres_score; 
	uint32_t *listE = NULL;
	uint32_t n_e=0;
	for (uint32_t e = 0; e < g->n_e; ++e) {
		uint32_t len = get_edge_len(&g->edges[e]);
		if (len > thres_len_e) {
			++n_e;
			listE = realloc(listE, n_e*sizeof(uint32_t));
			listE[n_e-1] = e; 
		}
	}
	__VERBOSE("n_e: %d\n", n_e);
	for (uint32_t i = 0; i < n_e; i++) {
		__VERBOSE("%d\n",i);
		uint32_t e0 = listE[i];
		for (uint32_t i1 = 0; i1 < n_e; i1++) {
			uint32_t e1 = listE[i1];
			struct bucks_score score = get_score_edges_res(e0, e1, g, n_bucks);
			uint32_t check = 0;
			if (e0 == e1){
				float cvr = get_genome_coverage(g);
				if (__get_edge_cov(&g->edges[e0], g->ksize)/cvr > 1.8) {
					check = 1;
				}
			} else if (e0 == g->edges[e1].rc_id) {
			} else {
				check = 1;
			}
			if (check) {
				if (score.score > thres_score) {
					if (!check_replicate_contig_edge(g, e0, e1, n_bucks, thres_score*(n_bucks * n_bucks - 1))) {
						fprintf(out_file, "score: %f center:%f edge: %d %d\n", score.score, 0.0000001, e0, e1); 
					}
				}
			}
		}
	}
	fclose(out_file);
}

struct contig_edge {
	uint32_t src, des;
	uint32_t rv_src, rv_des;
	float score0, score1;
};

void swap(uint32_t *a, uint32_t *b)
{ 
	uint32_t tmp;
	tmp = *a; 
	*a = *b; 
	*b = tmp;
}

void normalize_min_index(struct asm_graph_t *g, struct contig_edge *e)
{
	uint32_t rc_id_src = g->edges[e->src].rc_id;
	if (rc_id_src < e->src) {
		e->src = rc_id_src;
		e->rv_src ^= 1;
	}
	uint32_t rc_id_des = g->edges[e->des].rc_id;
	if (rc_id_des < e->des) {
		e->des = rc_id_des ;
		e->rv_des ^= 1;
	}
	if (e->src > e->des) {
		swap(&e->src, &e->des);
		swap(&e->rv_src, &e->rv_des);
		e->rv_src ^= 1 ;
		e->rv_des ^= 1;
	}
}

void normalize_one_dir(struct asm_graph_t *g, struct contig_edge *e)
{
	if (e->rv_src == 1) {
		e->src = g->edges[e->src].rc_id;
		e->rv_src = 0;
	}
	if (e->rv_des == 1) {
		e->des = g->edges[e->des].rc_id;
		e->rv_des = 0;
	}
}

void add_contig_edge(struct asm_graph_t *g,struct contig_edge *listE, uint32_t pos, uint32_t src, uint32_t des, float score0, float score1)
{
	struct contig_edge e;
	e.src = src;
	e.des = des;
	e.score0 = score0;
	e.score1 = score1;
	e.rv_src = 0;
	e.rv_des = 0;
	normalize_min_index(g, &e);
	listE[pos] = e;
}

int less_contig_edge(const void *e0,const void *e1)
{
	struct contig_edge * v0 = (struct contig_edge *) e0;
	struct contig_edge * v1 = (struct contig_edge *) e1;
	return (v0->src > v1->src || (v0->src == v1->src && v0->des > v1->des));
}

int less_uint32(const void *e0,const void *e1)
{
	return *(uint32_t*)(e0) > *(uint32_t*)(e1);
}

uint32_t equal_contig_edge(struct contig_edge *e0,struct contig_edge *e1)
{
	return (e0->src == e1->src && e0->des == e1->des);
}

uint32_t better_contig_edge(struct contig_edge *e0, struct contig_edge *e1)
{
	return (e0->score0 > e1->score0);
}

void unique_edge(struct contig_edge *listE, uint32_t *n_e)
{
	uint32_t new_n_e = 0;
	for (uint32_t i = 0; i < *n_e; ) {
		uint32_t j = i;
		struct contig_edge best;
		best.score0 = 0;
		while (equal_contig_edge(&listE[j], &listE[i])){
			if (better_contig_edge(&listE[j], &best)){
				best = listE[j];
			}
			++j;
		}
		listE[new_n_e++] = best;
		i = j;
	}
	*n_e = new_n_e;
}

uint32_t *binary_search(uint32_t *list, uint32_t n, uint32_t value)
{
	uint32_t l = 0, r = n - 1;
	while (l != r) {
		uint32_t m = (l + r) / 2;
		if (value > list[m]) 
			l = m + 1;
		else
			r = m;
	}
	if (list[l] != value)
		return list-1;
	return &list[l];
}

void unique_vertex(uint32_t *listV, uint32_t *n_v)
{
	uint32_t new_n_v = 0;
	for (uint32_t i = 0; i < *n_v; ) {
		uint32_t j = i;
		while (listV[i] == listV[j]) j++;
		listV[new_n_v++] = listV[i];
		i = j;
	}
	*n_v = new_n_v;
}

// **listV
void build_V_from_E(struct contig_edge *listE, uint32_t n_e, uint32_t **listV, uint32_t *n_v)
{
	*listV = NULL; 
	*n_v = 0;
	for (uint32_t i = 0; i < n_e; i++) {
		__VERBOSE("%d %d BBBB\n", listE[i].src, listE[i].des);

		++(*n_v);
		uint32_t t = (*n_v)*sizeof(uint32_t);
		*listV = realloc(*listV, t);
		(*listV)[(*n_v)-1] = listE[i].src;

		++(*n_v);
		*listV = realloc(*listV, (*n_v)*(sizeof(uint32_t)));
		(*listV)[(*n_v)-1] = listE[i].des;
	}
	for (uint32_t i = 0; i < *n_v; i++) {
		__VERBOSE("%d ", (*listV)[i]);
	}
	__VERBOSE("\n");
	qsort(*listV, *n_v, sizeof(uint32_t), less_uint32);
	for (uint32_t i = 0; i < *n_v; i++) {
		__VERBOSE("%d ", (*listV)[i]);
	}
	__VERBOSE("\n");
	qsort(*listV, *n_v, sizeof(uint32_t), less_uint32);
	unique_vertex(*listV, n_v);
//	__VERBOSE("listV %d\n", listV);
}

void print_seq(FILE *fp, uint32_t index, char *seq, uint32_t len, uint32_t cov)
{
	fprintf(fp, ">SEQ_%lld_length_%lld_count_%llu\n", (long long)index,
		(long long)len, (long long unsigned)cov);
	gint_t k = 0;
	char *buf = alloca(81);
//	printf("%s\n", seq);
	while (k < len) {
		gint_t l = __min(80, len - k);
		memcpy(buf, seq + k, l);
		buf[l] = '\0';
		fprintf(fp, "%s\n", buf);
		k += l;
	}
	while (k < len) {
		gint_t l = __min(80, len - k);
		memcpy(buf, seq + k, l);
		buf[l] = '\0';
		fprintf(fp, "%s\n", buf);
		k += l;
	}
}

void find_hamiltonian_contig_edge(FILE *out_file, struct asm_graph_t *g, struct contig_edge *listE_ori, uint32_t n_e){
	struct contig_edge *list_one_dir_E = calloc(2*n_e, sizeof(struct contig_edge));
	for (uint32_t i = 0; i < n_e; i++) {
		list_one_dir_E[i] = listE_ori[i];
		normalize_one_dir(g, list_one_dir_E+i);
	}
	for (uint32_t i = 0; i < n_e; i++) {
		list_one_dir_E[i + n_e] = listE_ori[i];
		list_one_dir_E[i + n_e].src = g->edges[list_one_dir_E[i].des].rc_id;
		list_one_dir_E[i + n_e].des = g->edges[list_one_dir_E[i].src].rc_id;
	}
	n_e *=2;
	uint32_t *listV = NULL , n_v; 
	build_V_from_E(list_one_dir_E, n_e, &listV, &n_v);
	__VERBOSE("n_v: %d \n", n_v);
	for (uint32_t i = 0; i < n_v; i++) {
		__VERBOSE("%d ", listV[i]);
	}

	uint32_t m = 0;
	uint32_t *E = calloc(n_e * 2, sizeof(uint32_t)), *next = calloc(n_e * 2, sizeof(uint32_t));
	int *head = calloc(n_v, sizeof(uint32_t));
	for (uint32_t i = 0; i < n_v; i++) {
		head[i] = -1;
	}
	for (uint32_t i = 0; i < n_e; i++) {
		int u = list_one_dir_E[i].src;
		u = binary_search(listV, n_v, u) - listV;
		int v = list_one_dir_E[i].des;
		v = binary_search(listV, n_v, v) - listV;
		__VERBOSE("edge: %d %d \n", u, v);
		if (u == -1 || v == -1){
			__VERBOSE("%d %d ERRRRRR\n", list_one_dir_E[i].src, list_one_dir_E[i].des);
		}
		E[m] = v;
		next[m] = head[u];
		head[u] = m;
		m++;
	}
	__VERBOSE("\nE:\n"); 
	__VERBOSE("\n"); 
	uint32_t *res = calloc(n_v, sizeof(uint32_t)), n_res=0;
	algo_find_hamiltonian(out_file, g,E, m, head, next, n_v, res, &n_res, listV);
	// todo free all pointer
}

void print_gfa_from_E(struct asm_graph_t *g, struct contig_edge *listE, uint32_t n_e)
{
	struct contig_edge *list_one_dir_E = calloc(2*n_e, sizeof(struct contig_edge));
	for (uint32_t i = 0; i < n_e; i++) {
		list_one_dir_E[i] = listE[i];
		normalize_min_index(g, list_one_dir_E+i);
	}
	uint32_t *listV = NULL , n_v; 
	build_V_from_E(list_one_dir_E, n_e, &listV, &n_v);
	FILE *fp = fopen("gr" , "w");
	for (uint32_t i = 0; i < n_v; i++) {
		fprintf(fp,"S\t%d\tAAA\tKC:i:1\n", listV[i]);
//		__VERBOSE("S\t%d\tAAA\tKC:i:1\n", g->edges[listV[i]].rc_id);
	}

	for (uint32_t i = 0; i < n_e; i++) {
		fprintf(fp, "L\t%d\t%c\t%d\t%c\t45M\n", 
			list_one_dir_E[i].src, list_one_dir_E[i].rv_src == 0?'+':'-', list_one_dir_E[i].des, list_one_dir_E[i].rv_des == 0?'+':'-');
		fprintf(fp, "L\t%d\t%c\t%d\t%c\t45M\n", 
			list_one_dir_E[i].des, list_one_dir_E[i].rv_des == 0?'-':'+', list_one_dir_E[i].src, list_one_dir_E[i].rv_src == 0?'-':'+');
	}
}

void build_graph_2(FILE *fp, FILE *out_file, struct asm_graph_t *g)
{
	for (uint32_t i = 0 ; i < 5; ++i) {
		uint32_t t = fscanf(fp, "%*[^\n]\n");
	}
	float score, center;
	uint32_t src, des;
	struct contig_edge *listE = NULL;
	uint32_t n_e=0;
	while (fscanf(fp,"score: %f center:%f edge: %d %d\n", &score, &center, &src, &des) !=EOF){
//		__VERBOSE("%f", score);
		n_e += 1;
		listE = realloc(listE, n_e*sizeof(struct contig_edge));
		add_contig_edge(g, listE, n_e-1, src, des, score, center);
	}
//	__VERBOSE("n_e: %d\n", n_e);
	qsort(listE, n_e, sizeof(struct contig_edge), less_contig_edge);
	unique_edge(listE, &n_e);

	// finding hamiltonian path
	find_hamiltonian_contig_edge(out_file, g, listE, n_e);

//	print_gfa_from_E(g, listE, n_e);
	// print as gfa graph
//	__VERBOSE("n_v: %d \n",n_v);
//	__VERBOSE("new_n_v: %d \n",n_v);
}

