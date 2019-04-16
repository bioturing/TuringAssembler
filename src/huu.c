#include <stdlib.h>
#include "assembly_graph.h"
#include "k31hash.h"
#include "pthread.h"
#include "verbose.h"
#include <string.h>
#include "algorithm.h"
#include <assert.h>

const float global_thres_score = 0.015;
const int global_thres_length = 10000;
const int global_thres_length_min = 5000;
const int global_thres_n_buck_big_small = 5;
int const global_n_buck = 6;

struct bucks_score {
	float score;
};

int min(int a, int b)
{
	if (a < b) return a;
	else return b;
}

int max(int a, int b)
{
	if (a<b)
		return b;
	return a;
}


uint32_t roundint(float x)
{ 
	return (x)>=0?(int)((x)+0.5):(int)((x)-0.5);
}

int get_amount_hole(struct asm_graph_t *g, struct asm_edge_t *e, uint32_t b)
{
	int res = 0, l = b * g->bin_size, r = (b+1) * g->bin_size, sum_holes = 0;
	for (uint32_t i = 0; i < e->n_holes; ++i){
		int pos = e->p_holes[i] + sum_holes;
		if (pos >= r){
			break;
		}
		if (pos >= l) {
			res += min(r-pos, e->l_holes[i]);
		}
		sum_holes += e->p_holes[i];
	}
	return res;
}

float get_score_bucks(struct asm_graph_t *g, struct asm_edge_t * e0, struct asm_edge_t *e1, uint32_t b0, uint32_t b1, float avg_bin_hash)
{
	const uint32_t thres_cnt = 30;
	struct barcode_hash_t *buck0 = &e0->bucks[b0], *buck1 = &e1->bucks[b1];
	uint32_t cnt0 = 0, cnt1 = 0, res2 = 0, cntss = 0;
	if (get_amount_hole(g, e0, b0)  > 0.7*g->bin_size) {
		__VERBOSE("hole \n");
		return -1;
	}
	if (get_amount_hole(g, e1, b1)  > 0.7*g->bin_size){
		__VERBOSE("hole \n");
		return -1;
	}
	for (uint32_t i = 0; i < buck1->size; ++i) {
		cntss++;
		if (buck1->cnts[i] != (uint32_t)(-1)) {
			if (buck1->cnts[i] >= thres_cnt) 
				cnt1++;
		}
	}

	for (uint32_t i = 0; i < buck0->size; ++i) {
		cntss++;
		if ((buck0->keys[i]) != (uint64_t)(-1)){
			cntss++;
			if (buck0->cnts[i] >= thres_cnt) {
				cnt0++;
				uint32_t tmp = barcode_hash_get(buck1, buck0->keys[i]);
				if (tmp != BARCODE_HASH_END(buck1) && buck1->cnts[tmp] >= thres_cnt) {
					res2++;
				}
			}
		}
	}
	if (avg_bin_hash/2 > cnt0 || cnt0 > avg_bin_hash * 2
		||avg_bin_hash/2 > cnt1 || cnt1 > avg_bin_hash * 2) {
//		__VERBOSE("extreme case\n");
		return -1;
	}
//	if (res2 == 0) __VERBOSE("res2==0 %d %d\n", cnt0 , cnt1);
	return 1.0 * res2 / (cnt0 + cnt1);
}

struct matrix_score{
	int n_bucks;
	float *A;
};

struct matrix_score *get_score_edges_matrix(struct asm_graph_t *g, uint32_t i0, uint32_t i1, int n_bucks, float avg_bin_hash)
{
	struct matrix_score *score = NULL;
	score = realloc(score, sizeof(struct matrix_score));
	uint32_t rev_i0 = g->edges[i0].rc_id;
	struct asm_edge_t *rev_e0 = &g->edges[rev_i0], *e1 = &g->edges[i1];
	int n0_bucks = (get_edge_len(rev_e0) + g->bin_size-1) / g->bin_size;
	int n1_bucks = (get_edge_len(e1) + g->bin_size-1) / g->bin_size;
	assert(n0_bucks > n_bucks);
	assert(n1_bucks > n_bucks);

	score->n_bucks = n_bucks;
	score->A = NULL;
	score->A = realloc(score->A, n_bucks * n_bucks * sizeof(float));
	for (int i = 0; i < score->n_bucks; ++i) {
		for (int j = 0; j < score->n_bucks; ++j) {
			score->A[i*score->n_bucks+j] = get_score_bucks(g, rev_e0, e1, i, j, avg_bin_hash);
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
	int count_dia = 0, count_all = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				float sc = score->A[i*n+j];
				if (sc > -0.000001){
					sum_dia += sc;
					count_dia++;
				}
			} else {
				float sc = score->A[i*n+j];
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

int check_replicate_contig_edge(struct asm_graph_t *g, uint32_t i0, uint32_t i1, const int n_bucks, float threshold, float avg_bin_hash)
{
	uint32_t rev_i0 = g->edges[i0].rc_id, rev_i1 = g->edges[i1].rc_id;
	struct matrix_score *s0 = get_score_edges_matrix(g, i0, i1, n_bucks, avg_bin_hash);
	struct matrix_score *s1 = get_score_edges_matrix(g, i0, rev_i1, n_bucks, avg_bin_hash);
	struct matrix_score *s2 = get_score_edges_matrix(g, rev_i0, i1, n_bucks, avg_bin_hash);
	struct matrix_score *s3 = get_score_edges_matrix(g, rev_i0, rev_i1, n_bucks, avg_bin_hash);
	return (detect_anomal_diagonal(s0, threshold) || detect_anomal_diagonal(s1, threshold) || detect_anomal_diagonal(s2, threshold) || detect_anomal_diagonal(s3, threshold));
}

struct bucks_score get_score_edges_res(uint32_t i0, uint32_t i1, struct asm_graph_t *g, const int n_bucks, float avg_bin_hash) {
	struct matrix_score *mat_score = get_score_edges_matrix(g, i0, i1, n_bucks, avg_bin_hash);
	mat_score->A[0] = -1;
	struct asm_edge_t *e0 = &g->edges[i0], *e1 = &g->edges[i1];
	float res = 0;
	int count = 0;
	for (int i = 0; i < mat_score->n_bucks; ++i) {
		for (int j = 0; j < mat_score->n_bucks; ++j) {
			float tmp = mat_score->A[i * n_bucks + j];
//			__VERBOSE("%f ", tmp);
			if (tmp >= -0.000001) {
				count++;
				res += tmp;
			}
		}
//		__VERBOSE("\n");
	}
	struct bucks_score res_score;
	res_score.score = res/count;
//	__VERBOSE("score %d %d\n", i0, i1); 
	return res_score;
}

float abssss(float x)
{
	if (x < 0)
		x = -x;
	return x;
}

int get_score_big_small(int i0, int i1, struct asm_graph_t *g, float avg_bin_hash) {
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
	for (int i = max(0, n0_bucks-20); i < n0_bucks-1; ++i) {
		float left_value = 0, right_value = 0, count_left = 0, count_right = 0;
		for (int j = 0; j < 3; ++j) {
			float tmp = 0;
			tmp = get_score_bucks(g, e0, e1, i, j, avg_bin_hash);
			if (tmp > -0.000001) {
				left_value += tmp;
				count_left++;
			}
		}
		for (int j = n1_bucks - 4; j < n1_bucks-1; ++j) {
			float tmp = 0;
			tmp = get_score_bucks(g, e0, e1, i, j, avg_bin_hash);
			if (tmp > -0.000001) {
				right_value += tmp;
				count_right++;
			}
//			__DEBUG_VERBOSE("%f ", tmp);
		}
		if ( i0 == 3005 && i1== 5321){
			__VERBOSE("caseee %f %f\n",  left_value, right_value);
		}
		left_value /= count_left;
		right_value /= count_right;
//		__VERBOSE("%f %f", left_value , right_value);
		if (abssss(left_value - right_value) > 0.05 * max(left_value, right_value)) {
//			__VERBOSE("%f", abssss(left_value - right_value));
			if (left_value > right_value)
				score++;
			else 
				score--;
		}
	}
	return score;
}

void check_contig(struct asm_graph_t *g, float avg_bin_hash) {
	int cmp(const void *i, const void *j)
	{
		uint32_t x = *(uint32_t *)i;
		uint32_t y = *(uint32_t *)j;
		return get_edge_len(&g->edges[x]) > get_edge_len(&g->edges[y]);
	}
	uint32_t *listE = NULL;
	uint32_t n_e= g->n_e;
	listE = realloc(listE, n_e * sizeof(uint32_t));
	for (uint32_t e = 0; e < n_e; ++e) {
		listE[e] = e;
	}
	qsort(listE, n_e, sizeof(4), cmp);

	FILE * f = fopen("list_pair.txt","r");
	uint32_t a, b, asdf;
	char *s = NULL;
	s = realloc(s, 100);
	while (fscanf(f, "%d %d %d\n", &a, &b, &asdf) != EOF) {
		float score = get_score_big_small(a, b, g, avg_bin_hash);
		__DEBUG_VERBOSE("%d %d edge length: %d %d score:%f %f\n", a, b, get_edge_len(&g->edges[a]), get_edge_len(&g->edges[b]), score, 0.000000);
	}
	fclose(f);
}

float count_bin_hash(struct asm_graph_t *g, struct barcode_hash_t *buck)
{
	const uint32_t thres_cnt = 30;
	uint32_t cnt = 0;
	
	for (uint32_t i = 0; i < buck->size; ++i) {
		if (buck->cnts[i] != (uint32_t)(-1)) {
			if (buck->cnts[i] >= thres_cnt) 
				cnt++;
		}
	}
	return cnt;
}

float get_avg_bin_hash(struct asm_graph_t *g)
{
	int count = 0;
	long long sum = 0;
	for (uint32_t i = 0; i < g->n_e; ++i) {
		int n_bucks = (get_edge_len(&g->edges[i]) + g->bin_size-1) / g->bin_size;
		for (int j = 0; j < n_bucks - 1; j++){
			int tmp = count_bin_hash(g, &g->edges[i].bucks[j]);
			if (tmp > 0) {
				count++;
				sum += tmp;
			}
		}
	}
	return 1.0*sum/count;
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
	fprintf(out_file, "n_v: %d\n", n_e);
	for (uint32_t i = 0; i < n_e ; i++) {
		uint32_t e = listE[i];
		fprintf(out_file, "vertex:%d\n", e);
	}

	float avg_bin_hash = get_avg_bin_hash(g);
	__VERBOSE("avg_bin_hash %f", avg_bin_hash);
	__VERBOSE("n_e: %d\n", n_e);
	for (uint32_t i = 0; i < n_e; i++) {
		__VERBOSE("%d\n",i);
		uint32_t e0 = listE[i];
		for (uint32_t i1 = 0; i1 < n_e; i1++) {
			uint32_t e1 = listE[i1];
			struct bucks_score score = get_score_edges_res(e0, e1, g, n_bucks, avg_bin_hash);
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
					if (!check_replicate_contig_edge(g, e0, e1, n_bucks, thres_score*(n_bucks * n_bucks - 1), avg_bin_hash)) {
						fprintf(out_file, "score: %f edge: %d %d\n", score.score, e0, e1); 
					}
				}
			}
			__VERBOSE("score %f edge: %d %d\n", score.score, e0, e1); 
		}
	}
	fclose(out_file);
}

struct contig_edge {
	uint32_t src, des;
	uint32_t rv_src, rv_des;
	float score0;
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

void add_contig_edge(struct asm_graph_t *g,struct contig_edge *listE, uint32_t pos, uint32_t src, uint32_t des, float score0)
{
	struct contig_edge e;
	e.src = src;
	e.des = des;
	e.score0 = score0;
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
	qsort(*listV, *n_v, sizeof(uint32_t), less_uint32);
	for (uint32_t i = 0; i < *n_v; i++) {
		__VERBOSE("%d ", (*listV)[i]);
	}
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

void algo_find_hamiltonian(FILE *out_file, struct asm_graph_t *g, uint32_t *E, uint32_t n_e, int *head, uint32_t *next, uint32_t n_v, uint32_t *res, uint32_t *n_res, uint32_t *listV, float avg_bin_hash)
{
	void print_contig(int index, uint32_t n_contig, uint32_t *list_contig)
	{
		char *seq = NULL, *total_seq = NULL, *NNN = NULL;
		NNN = calloc(1000, sizeof(char));
		for (uint32_t i = 0; i < 1000; i++) 
			NNN[i] = 'N';
		uint32_t seq_len = 0, total_len = 0;
		for(uint32_t i = 0; i < n_contig; i++) {
			uint32_t e = list_contig[i];
			uint32_t len = dump_edge_seq(&seq, &seq_len, &g->edges[e]);
			total_seq = realloc(total_seq, (total_len + len) * sizeof(char));
			memcpy(total_seq+total_len, seq, len);
			total_len += len;
		}
		print_seq(out_file, index, total_seq, total_len, 1);
		for(uint32_t i = 0; i < n_contig; i++) {
			uint32_t e = list_contig[i];
			__VERBOSE("%d ", e);
		}
	}

	void insert_short_contigs(uint32_t n_big_contigs, uint32_t *big_contigs, uint32_t n_insert, int *arr_insert, 
			int n_short, int *arr_i_short, int *mark_short)
	{
		//todo insert to the beginning
		__VERBOSE("big contig length %d:", n_big_contigs);
		for (int i = 0; i < n_big_contigs; i++){
			__VERBOSE("%d ", big_contigs[i]);
		}
		for (uint32_t i = 1; i < n_insert - 1; i++) {
			__VERBOSE("i insert: %d\n", i);
			int max_score = -1, pos = -1;
			for (int j = 0; j < n_short; j++) if (mark_short[j] == 0) {
				int score0 = get_score_big_small(big_contigs[i-1], arr_i_short[j], g, avg_bin_hash);
				int score1 = get_score_big_small(g->edges[big_contigs[i]].rc_id, g->edges[arr_i_short[j]].rc_id, g, avg_bin_hash);
				int score = score0 + score1;
				if (score > max_score){
					max_score = score;
					pos = j;
				}
			}
			__VERBOSE("max score: %d \n",max_score);
			//todo check if normal score large enough
			if (max_score > 15) {
				arr_insert[i] = arr_i_short[pos];
				mark_short[pos] = 1;
			}
		}
		//todo insert to the end
	}

	void add_insert_except_m1(int *n_arr, uint32_t **new_arr, uint32_t ele)
	{
		if (ele == -1) 
			return;
		*new_arr = realloc(*new_arr, (*n_arr+1) * sizeof(uint32_t));
		(*new_arr)[*n_arr] = ele;
		++(*n_arr);
	}

	void merge_big_and_small(uint32_t *best_n_res, uint32_t **best_res, uint32_t n_insert, int *arr_insert)
	{
		uint32_t *new_arr = NULL;
		int n_arr = 0;
		add_insert_except_m1(&n_arr, &new_arr, arr_insert[0]);
		for (uint32_t i = 0; i < *best_n_res; i++) {
			add_insert_except_m1(&n_arr, &new_arr, (*best_res)[i]); 
			add_insert_except_m1(&n_arr, &new_arr, arr_insert[i]); 
		}
		*best_n_res = n_arr;
		*best_res = new_arr;
		__VERBOSE("after merge contig length %d:", *best_n_res);
		for (int i = 0; i < *best_n_res; i++){
			__VERBOSE("%d ", (*best_res)[i]);
		}
	}

	int thres_len = global_thres_length, thres_len_min = global_thres_length_min;

	int *arr_i_short = NULL, n_arr_short = 0;
	int *mark_short = NULL;
	for (uint32_t e = 0; e < g->n_e; ++e) {
		int len = get_edge_len(&g->edges[e]);
		if (len < thres_len && len > thres_len_min) {
			++n_arr_short;
			arr_i_short = realloc(arr_i_short, n_arr_short * sizeof(uint32_t));
			mark_short = realloc(mark_short, n_arr_short * sizeof(uint32_t));
			arr_i_short[n_arr_short-1] = e;
			mark_short[n_arr_short-1] = 0;
		}
	}

	int *mark = calloc(n_v, sizeof(uint32_t));
	for (uint32_t i = 0; i < n_v; i++) {
		float cvr = get_genome_coverage(g);
		float cov_times = (__get_edge_cov(&g->edges[listV[i]], g->ksize)/cvr) ;
		mark[i] = roundint(cov_times);
		__VERBOSE("%d " , mark[i]);
	}
	uint32_t count = 0;
	for(uint32_t ii = 0; ii < n_v; ii++){
		uint32_t *best_res = calloc(n_v, sizeof(uint32_t)), best_n_res = 0 ;
		for (uint32_t i = 0; i < n_v; i++) if (mark[i] > 0) {
//			__VERBOSE("dfs from %d %d\n", i, mark[i]);
			dfs_hamiltonian(i,  1, E, head, next, mark, n_v, res, best_res, &best_n_res);
		}
		if (best_n_res == 0) 
			break;

		for (uint32_t i = 0; i < best_n_res; i++) {
			mark[best_res[i]]--;
			mark[best_res[i]^1]--;
			//__VERBOSE("%d %d \n", listV[best_res[i]], mark[listV[best_res[i]]] );
		}
		int *arr_insert = NULL, n_insert = best_n_res + 1;
		arr_insert = realloc(arr_insert, n_insert * sizeof(uint32_t));
		for (int i = 0; i < n_insert; i++) arr_insert[i] = -1;
		for (uint32_t i = 0; i < best_n_res; i++) 
			best_res[i] = listV[best_res[i]];
		insert_short_contigs(best_n_res, best_res, n_insert, arr_insert, n_arr_short, arr_i_short, mark_short);
		merge_big_and_small(&best_n_res, &best_res, n_insert, arr_insert);

		print_contig(ii, best_n_res, best_res);
		__VERBOSE("best n %d\n", best_n_res);
//		__VERBOSE("\nmark\n");
		count++;
	}

	for (int i = 0; i < n_arr_short; i++) if (mark[i] == 0) {
		uint32_t e = arr_i_short[i]; 
		uint32_t seq_len = 0;
		char *seq = NULL;
		uint32_t len = dump_edge_seq(&seq, &seq_len, &g->edges[e]);
		print_seq(out_file, count, seq, seq_len, 1); 
		count++;
	}

	for (int e = 0; e < g->n_e; e++){
		int len  = get_edge_len(&g->edges[e]);
		if (len < thres_len_min) {
			uint32_t seq_len = 0;
			char *seq = NULL;
			uint32_t len = dump_edge_seq(&seq, &seq_len, &g->edges[e]);
			print_seq(out_file, count, seq, seq_len, 1); 
			count++;
		}
	}

	fclose(out_file);
}

void find_hamiltonian_contig_edge(FILE *out_file, struct asm_graph_t *g, struct contig_edge *listE_ori, uint32_t n_e, uint32_t n_v, uint32_t *listV, float avg_bin_hash){
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
	algo_find_hamiltonian(out_file, g,E, m, head, next, n_v, res, &n_res, listV, avg_bin_hash);
	// todo free all pointer
}

void print_gfa_from_E(struct asm_graph_t *g, struct contig_edge *listE, uint32_t n_e, uint32_t *listV, uint32_t n_v)
{
	struct contig_edge *list_one_dir_E = calloc(2*n_e, sizeof(struct contig_edge));
	for (uint32_t i = 0; i < n_e; i++) {
		list_one_dir_E[i] = listE[i];
		normalize_min_index(g, list_one_dir_E+i);
	}
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
	uint32_t n_v, *listV = NULL;
	fscanf(fp, "n_v: %d\n", &n_v);
	listV = realloc(listV , n_v * sizeof(uint32_t));
	for (uint32_t i = 0; i < n_v; i++) {
		fscanf(fp, "vertex:%d\n", &listV[i]);
	}
	float score;
	uint32_t src, des;
	struct contig_edge *listE = NULL;
	uint32_t n_e=0;
	float avg_bin_hash = get_avg_bin_hash(g);
	while (fscanf(fp,"score: %f edge: %d %d\n", &score, &src, &des) !=EOF){
		n_e += 1;
		listE = realloc(listE, n_e*sizeof(struct contig_edge));
		add_contig_edge(g, listE, n_e-1, src, des, score);
	}
	qsort(listE, n_e, sizeof(struct contig_edge), less_contig_edge);
	unique_edge(listE, &n_e);

	find_hamiltonian_contig_edge(out_file, g, listE, n_e, n_v, listV, avg_bin_hash);
}

