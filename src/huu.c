#include <stdlib.h>
#include "assembly_graph.h"
#include "k31hash.h"
#include "pthread.h"
#include "verbose.h"
#include <string.h>
#include "algorithm.h"
#include <assert.h>
#include "compare.h"
#include "contig_graph.h"
#include "math.h"

#define LIST_GLOBAL_PARAMS \
	X(float, global_thres_bucks_score, -1)\
	X(int, global_thres_count_kmer , -1)\
	X(int, global_thres_length , -1)\
	X(int, global_thres_length_min , -1)\
	X(int, global_thres_n_buck_big_small , -1)\
	X(int, global_n_buck , -1)\
	X(float, global_genome_coverage, -1)\
	X(int, global_molecule_length, -1)\
	X(float, global_thres_coefficent, -1); 

// constant for logging
int log_debug = 1;
int log_bin_score = 0;
int log_share_barcode = 0;
int log_global_var = 0;
int log_hole = 0;
int log_outliner = 0;
int log_beautify = 0;
int log_edge_score = 0;
int log_build_contig = 0;
int log_build_V = 0;
int log_count_kmer = 0;
int log_check_contig = 1;
int log_insert_small_contig = 0;
int log_hamiltonian = 0;
int log_graph = 1;
int log_assert = 1;
int log_distance_score = 1;

#define X(type, name, default_value) type name=default_value;
LIST_GLOBAL_PARAMS
#undef X

void check_global_params()
{
#define X(type, name, default_value) assert((name) != default_value);
LIST_GLOBAL_PARAMS
#undef X
}

#define __VERBOSE_FLAG(flag, ...) if (flag) { __VERBOSE("["#flag"] "); __VERBOSE(__VA_ARGS__);}

#define __COPY_ARR(src, des, n) for(int i = 0; i < n; i++) des[i] = src[i]

struct bucks_score {
	float score;
};

struct matrix_score{
	int n_bucks;
	float *A;
};

void destroy_matrix_score(struct matrix_score *x)
{
	free(x->A);
	free(x);
}

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

int roundint(float x)
{ 
	return (x)>=0?(int)((x)+0.5):(int)((x)-0.5);
}

float get_global_thres_score(struct asm_graph_t *g)
{
	float cvr = get_genome_coverage(g);
	float res = 1.0 * g->bin_size / global_molecule_length * global_thres_coefficent;
	__VERBOSE_FLAG(log_global_var, "global thres score: %f\n", res);
	return res;
}

int get_global_count_kmer(struct asm_graph_t *g)
{
	int *arr_count = NULL, n_arr = 0, *count_count = NULL; 
	 int res = -1;
	for (int i = 0; i < g->n_e; i++) {
		int n_bucks = (get_edge_len(&g->edges[i]) + g->bin_size-1) / g->bin_size;
		for (int j = 0; j < n_bucks-1; j++){
			struct barcode_hash_t buck = g->edges[i].bucks[j];
			for (uint32_t l = 0; l < buck.n_item; l++){
				arr_count = realloc(arr_count, (n_arr + 1) * sizeof(int));
				arr_count[n_arr] = buck.cnts[l];
				n_arr++;
			}
		}
	}
	int size_count = 1000;
	count_count = realloc(count_count, size_count*sizeof(int)); 
	for(int i = 0; i < n_arr; i++) {
		__VERBOSE_FLAG(log_assert, "%d %d", arr_count[i], size_count);
		assert(arr_count[i] >  size_count);
		count_count[arr_count[i]]++;
	}
	int max_count = 0;
	for (int i = 1; i < 10 ; i++){
		max_count = max(max_count, count_count[i]);
	}
	__VERBOSE_FLAG(log_count_kmer, "max_count %d \n", max_count);
	for(int i = 10; i < size_count; i++) {
		__VERBOSE_FLAG(log_count_kmer, "%count hash ",  count_count[i]);
		if (count_count[i] > 3 * max_count){
			res = i;
			break;
		}
	}
	__VERBOSE_FLAG(log_count_kmer, "global thres count kmer: %d\n", res);
	__VERBOSE_FLAG(log_assert, "%d\n", res);
	assert(res == -1);
	return res;
}

void init_global_params(struct asm_graph_t *g, float huu_1_score)
{
	global_thres_length = 10000;
	global_thres_length_min = 5000;
	global_thres_n_buck_big_small = 5;
	global_n_buck = 6;
	global_molecule_length = 20000;
	global_thres_count_kmer =  1;//get_global_count_kmer(g);
	if (huu_1_score != -1) {
		global_thres_coefficent = huu_1_score;
	} else {
		global_thres_coefficent = 0.2;
	}
	global_genome_coverage = get_genome_coverage(g);
	global_thres_bucks_score = get_global_thres_score(g);
}

int get_amount_hole(struct asm_graph_t *g, struct asm_edge_t *e, int b)
{
	int res = 0, l = b * g->bin_size, r = (b+1) * g->bin_size, sum_holes = 0;
//	for (int i = 0; i < e->n_holes; ++i){
//		__VERBOSE_FLAG(log_hole, "holeee %d %d %d\n" , e->seq_len, e->l_holes[i], e->p_holes[i]);
//	}
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

int check_qualify_buck(struct asm_graph_t *g, struct asm_edge_t *e, int b, float avg_bin_hash)
{
	if (get_amount_hole(g, e, b)  > 0.7*g->bin_size) {
		__VERBOSE_FLAG(log_hole, "NNNNN size is to big ");
		return 0;
	}
	int cnt = 0, cov = global_genome_coverage, normal_count = (g->bin_size - g->ksize +1) * cov;

	struct barcode_hash_t *buck = &e->bucks[b];
	for (uint32_t i = 0; i < buck->size; ++i) {
		if (buck->cnts[i] != (uint32_t)(-1)) {
			cnt += buck->cnts[i];
		}
	}
	assert(normal_count != 0);
	if  (cnt > 2 * normal_count || cnt < normal_count * 0.5) 
	{
		__VERBOSE_FLAG(log_outliner, "count hash is abnormal: %d %d\n", cnt, normal_count);
		return 0;
	}
	return 1;
}

int check_count_hash_buck(struct asm_graph_t *g, struct asm_edge_t *e, struct barcode_hash_t *buck,  float avg_bin_hash)
{
	int cnt = 0, cov = global_genome_coverage, normal_count = (g->bin_size - g->ksize +1) * cov;
	for (uint32_t i = 0; i < buck->size; ++i) {
		if (buck->cnts[i] != (uint32_t)(-1)) {
			cnt += buck->cnts[i];
		}
	}
	assert(normal_count != 0);
	if  (cnt > 2 * normal_count || cnt < normal_count * 0.5) 
	{
		__VERBOSE_FLAG(log_outliner, "count hash is abnormal: %d %d\n", cnt, normal_count);
		return 0;
	}
	return 1;
}

float get_score_bucks(struct barcode_hash_t *buck0, struct barcode_hash_t *buck1)
{
	const int thres_cnt = global_thres_count_kmer;
	int cnt0 = 0, cnt1 = 0, res2 = 0, cntss = 0;

	for (uint32_t i = 0; i < buck1->size; ++i) {
		cntss++;
		if (buck1->cnts[i] != (uint32_t)(-1)) {
			if (buck1->cnts[i] >= (uint32_t)thres_cnt) 
				cnt1+= buck1->cnts[i];
		}
	}

	for (uint32_t i = 0; i < buck0->size; ++i) {
		cntss++;
		if ((buck0->keys[i]) != (uint64_t)(-1)){
			cntss++;
			if (buck0->cnts[i] >= (uint32_t)thres_cnt) {
				cnt0 += buck0->cnts[i];
				uint32_t tmp = barcode_hash_get(buck1, buck0->keys[i]);
				if (tmp != BARCODE_HASH_END(buck1) && buck1->cnts[tmp] >= (uint32_t)thres_cnt) {
					res2+= min(buck0->cnts[i] , buck1->cnts[tmp]);
				}
			}
		}
	}
	__VERBOSE_FLAG(log_share_barcode, "res %d cnt0 %d cnt1 %d \n", res2, cnt0, cnt1 );
	if (min(cnt0, cnt1) == 0) 
		return 0;
	return 1.0 * res2 / min(cnt0 , cnt1);
}

struct matrix_score *get_score_edges_matrix(struct asm_graph_t *g, int i0, int i1, int n_bucks, float avg_bin_hash)
{
	struct matrix_score *score = NULL;
	score = realloc(score, sizeof(struct matrix_score));
	int rev_i0 = g->edges[i0].rc_id;
	assert(rev_i0 < g->n_e);
	struct asm_edge_t *rev_e0 = &g->edges[rev_i0], *e1 = &g->edges[i1];
	int n0_bucks = (get_edge_len(rev_e0) + g->bin_size-1) / g->bin_size;
	int n1_bucks = (get_edge_len(e1) + g->bin_size-1) / g->bin_size;
	assert(n0_bucks > n_bucks);
	assert(n1_bucks > n_bucks);

	score->n_bucks = n_bucks;
	score->A = NULL;
	score->A = realloc(score->A, n_bucks * n_bucks * sizeof(float));
	// check_bucks_A[i] && check_bucks_B[i] de improve performance
	__VERBOSE_FLAG(log_global_var, "cov %f binsize %d ksize %d", global_genome_coverage, g->bin_size, g->ksize);
	int *check_bucks_A = NULL, *check_bucks_B = NULL;
	check_bucks_A = realloc(check_bucks_A, n_bucks * sizeof(int));
	check_bucks_B = realloc(check_bucks_B, n_bucks * sizeof(int));
	for (int i = 0; i < score->n_bucks; ++i) {
		check_bucks_A[i] = check_qualify_buck(g, rev_e0, i, avg_bin_hash);
		check_bucks_B[i] = check_qualify_buck(g, e1, i, avg_bin_hash);
	}
	for (int i = 0; i < n_bucks; ++i) {
		for (int j = 0; j < n_bucks; ++j) {
			if (check_bucks_A[i] && check_bucks_B[j]) 
				score->A[i*n_bucks+j] = get_score_bucks(&rev_e0->bucks[i], &e1->bucks[j]);
			else 
				score->A[i*n_bucks+j] = -1;
		}
	}
	__VERBOSE_FLAG(log_beautify, "`````````````\n");

	free(check_bucks_A);
	free(check_bucks_B);
	return score;
}

int detect_anomal_diagonal(struct matrix_score *score, float threshold)
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

int check_replicate_contig_edge(struct asm_graph_t *g, int i0, int i1, const int n_bucks, float threshold, float avg_bin_hash)
{
	int rev_i0 = g->edges[i0].rc_id, rev_i1 = g->edges[i1].rc_id;
	struct matrix_score *s0 = get_score_edges_matrix(g, i0, i1, n_bucks, avg_bin_hash);
	struct matrix_score *s1 = get_score_edges_matrix(g, i0, rev_i1, n_bucks, avg_bin_hash);
	struct matrix_score *s2 = get_score_edges_matrix(g, rev_i0, i1, n_bucks, avg_bin_hash);
	struct matrix_score *s3 = get_score_edges_matrix(g, rev_i0, rev_i1, n_bucks, avg_bin_hash);
	int res = (detect_anomal_diagonal(s0, threshold) || detect_anomal_diagonal(s1, threshold) || detect_anomal_diagonal(s2, threshold) || detect_anomal_diagonal(s3, threshold));
	destroy_matrix_score(s0);
	destroy_matrix_score(s1);
	destroy_matrix_score(s2);
	destroy_matrix_score(s3);
	return res;
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
	for (int i = max(0, n0_bucks-20); i < n0_bucks-1; ++i) {
		float left_value = 0, right_value = 0, count_left = 0, count_right = 0;
		for (int j = 0; j < 3; ++j) {
			float tmp = 0;
			tmp = get_score_bucks(&e0->bucks[i],& e1->bucks[j]);
			if (tmp > global_thres_bucks_score) {
				left_value += tmp;
				count_left++;
			}
		}
		for (int j = n1_bucks - 4; j < n1_bucks-1; ++j) {
			float tmp = 0;
			tmp = get_score_bucks(&e0->bucks[i], &e1->bucks[j]);
			if (tmp > global_thres_bucks_score) {
				right_value += tmp;
				count_right++;
			}
			__VERBOSE_FLAG(log_edge_score, "%f ", tmp);
		}
		left_value /= count_left;
		right_value /= count_right;
		__VERBOSE_FLAG(log_edge_score, "%f %f", left_value , right_value);
		if (fabsf(left_value - right_value) > 0.05 * max(left_value, right_value)) {
			__VERBOSE_FLAG(log_edge_score, "%f", fabsf(left_value - right_value));
			if (left_value > right_value)
				score++;
			else 
				score--;
		}
	}
	return score;
}

float count_bin_hash(struct asm_graph_t *g, struct barcode_hash_t *buck)
{
	const int thres_cnt = global_thres_count_kmer;
	int cnt = 0;
	
	for (uint32_t i = 0; i < buck->size; ++i) {
		if (buck->cnts[i] != (uint32_t)(-1)) {
			if (buck->cnts[i] >= (uint32_t)thres_cnt)
				cnt++;
		}
	}
	return cnt;
}

float get_avg_bin_hash(struct asm_graph_t *g)
{
	int count = 0;
	long long sum = 0;
	for (int i = 0; i < g->n_e; ++i) {
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

float get_score_multiple_buck(struct asm_graph_t *g, struct asm_edge_t *e, struct barcode_hash_t *b_left, struct barcode_hash_t *b_right)
{
	float avg_bin_hash = get_avg_bin_hash(g);
	float res = 0;
	int count = 0;
	for (int i = 0; i < 3; i++) if (check_count_hash_buck(g, e, &b_left[i], avg_bin_hash)) {
		for(int j = 0; j < 3; j++) if (check_count_hash_buck(g, e, &b_right[j], avg_bin_hash)) {
			count ++;
			float t = get_score_bucks(&b_left[i], &b_right[j]);
//			__VERBOSE_FLAG(log_check_contig, "score %f\n", t);
			res += t;
		}
	}
	if (count == 0)
		return 0;
	return res/count;
}

static pthread_mutex_t lock_merge = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t lock_id = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t lock_write_file = PTHREAD_MUTEX_INITIALIZER;

struct params{
	struct asm_graph_t *g;
	int i;
};

void *process(void *data)
{
	struct params *pa = (struct params *) data;
	do{
		int iiii;
		pthread_mutex_lock(&lock_id);
		if (pa->i == (pa->g)->n_e){
			pthread_mutex_unlock(&lock_id);
			break;
		}
		else {
			iiii = pa->i++;
			__VERBOSE("start %d\n", iiii);
		//	__VERBOSE("huu %d\n", i);
		}
		pthread_mutex_unlock(&lock_id);
		struct asm_graph_t *g = pa->g;
		struct asm_edge_t *e = &g->edges[iiii];
		int n_bucks = (get_edge_len(&g->edges[iiii]) + g->bin_size-1) / g->bin_size;
		if (n_bucks < 40 )
			continue;
		char *file_name = calloc(20, 1);
		sprintf(file_name, "logfile_%d", iiii);
		FILE *f = fopen(file_name, "w");
		for (int j = 0; j < min(20, n_bucks - 4); j+=5){
			for (int j1 = j + 4; j1 < n_bucks -4; j1++) {
	//			__VERBOSE("%d %d\n", j, j1);
				float tmp  = get_score_multiple_buck(pa->g, e, &e->bucks[j], &e->bucks[j1]);
				fprintf(f, "distance %d %f\n", j1 - j, tmp);
	//			__VERBOSE_FLAG(log_distance_score, "distance %d %f \n", j1 - j, tmp);
			}
		}
		fclose(f);
		pthread_mutex_lock(&lock_id);
		__VERBOSE("done %d\n", iiii);
		pthread_mutex_unlock(&lock_id);
	} while (1);
	pthread_exit(NULL);
	return NULL;
}

void check_contig(struct asm_graph_t *g, float avg_bin_hash, float huu_1_score) 
{
	int cmp(const void *i, const void *j)
	{
		int x = *(int *)i;
		int y = *(int *)j;
		return get_edge_len(&g->edges[x]) > get_edge_len(&g->edges[y]);
	}
	init_global_params(g, huu_1_score);
	check_global_params(g);
	__VERBOSE_FLAG(log_check_contig, "check contig\n");
	for (int i = 0; i < g->n_e; i++) {
		__VERBOSE("edge %d leng %d ", i, get_edge_len(&g->edges[i]));
		for (uint32_t j = 0; j < g->edges[i].n_holes; j++){
			__VERBOSE("hole %d ", g->edges[i].p_holes[j]);
		}
		__VERBOSE("\n");

	}

//	int n_threads = 5;
//	pthread_t *thr = (pthread_t *)calloc(n_threads, sizeof(pthread_t));
//	pthread_attr_t attr;
//	pthread_attr_init(&attr);
//	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
//	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
//
//	int count_thread = 0;
//	struct params *para =  calloc(1, sizeof(struct params));
//	para->g = g;
//	para->i = 0;
//	for (int i = 0; i < n_threads; ++i)
//		pthread_create(&thr[i], &attr, process, para);
//	for (int i = 0; i < n_threads; ++i)
//		pthread_join(thr[i], NULL);
//	free(thr);
//	pthread_attr_destroy(&attr);
}

struct bucks_score get_score_edges_res(int i0, int i1, struct asm_graph_t *g, const int n_bucks, float avg_bin_hash) 
{
	struct matrix_score *mat_score = get_score_edges_matrix(g, i0, i1, n_bucks, avg_bin_hash);
	mat_score->A[0] = -1;
	struct asm_edge_t *e0 = &g->edges[i0], *e1 = &g->edges[i1];
	float res = 0;
	int count = 0;
	for (int i = 0; i < mat_score->n_bucks; ++i) {
		for (int j = 0; j < mat_score->n_bucks; ++j) {
			float tmp = mat_score->A[i * n_bucks + j];
			__VERBOSE_FLAG(log_bin_score, "%f ", tmp);
			if (tmp >= -0.000001) {
				count++;
				res += tmp;
			}
		}
		__VERBOSE_FLAG(log_bin_score, "#\n");
	}
	struct bucks_score res_score;
	res_score.score = res/count;
	destroy_matrix_score(mat_score);
	return res_score;
}

struct params_check_edge {
	struct asm_graph_t *g;
	int i, j;
	int *listE, n_bucks, n_listE;
	float avg_bin_hash, thres_score;
	FILE *out_file;
};

void *process_check_edge(void *data)
{
	struct params_check_edge *pa_check_edge = (struct params_check_edge *) data;
	int *listE = pa_check_edge->listE;
	struct asm_graph_t *g = pa_check_edge->g;
	int n_bucks = pa_check_edge->n_bucks;
	float avg_bin_hash = pa_check_edge->avg_bin_hash;
	float thres_score = pa_check_edge->thres_score;
	FILE* out_file = pa_check_edge->out_file;
	int n_listE = pa_check_edge->n_listE;
	do {
		int i, j;
		pthread_mutex_lock(&lock_id);
		if (pa_check_edge->i == n_listE) {
			pthread_mutex_unlock(&lock_id);
			break;
		} else {
			i = pa_check_edge->i;
			j = pa_check_edge->j++;
			if (pa_check_edge->j == n_listE){
				pa_check_edge->i++;
				pa_check_edge->j = 0;
			}
		}
		assert(i < n_listE && j < n_listE);
		pthread_mutex_unlock(&lock_id);
		int e0 = listE[i];
		int e1 = listE[j];
		__VERBOSE_FLAG(log_build_contig, "edge: %d %d\n", e0, e1); 
		__VERBOSE_FLAG(log_assert, "edges: %d %d\n", e0, e1); 
		assert(e1 < g->n_e && e0 < g->n_e);
		struct bucks_score score = get_score_edges_res(e0, e1, g, n_bucks, avg_bin_hash);
		int check = 0;
		if (e0 == e1){
			float cvr = global_genome_coverage;
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
					pthread_mutex_lock(&lock_write_file);
					fprintf(out_file, "score: %f edge: %d %d\n", score.score, e0, e1); 
					pthread_mutex_unlock(&lock_write_file);
				}
			}
		}
		__VERBOSE_FLAG(log_build_contig, "score %f\n", score.score); 
	} while (1); 
	pthread_exit(NULL);
	return NULL;
}

void build_list_contig(struct asm_graph_t *g, FILE *out_file, float huu_1_score) 
{
	init_global_params(g, huu_1_score);
	check_global_params(g);
	const int thres_len_e = global_thres_length; 
	const int n_bucks = global_n_buck;
	float thres_score = global_thres_bucks_score; 
	int *listE = NULL;
	int n_e=0;

	for (int e = 0; e < g->n_e; ++e) {
		int len = get_edge_len(&g->edges[e]);
		if (len > thres_len_e) {
			++n_e;
			listE = realloc(listE, n_e*sizeof(int));
			listE[n_e-1] = e; 
		}
	}
	fprintf(out_file, "n_v: %d\n", n_e);
	for (int i = 0; i < n_e ; i++) {
		int e = listE[i];
		fprintf(out_file, "vertex:%d\n", e);
	}

	float avg_bin_hash = get_avg_bin_hash(g);
	__VERBOSE_FLAG(log_global_var, "avg_bin_hash %f", avg_bin_hash);
	__VERBOSE_FLAG(log_global_var, "n_e: %d\n", n_e);

	int n_threads = 56;
	pthread_t *thr = (pthread_t *)calloc(n_threads, sizeof(pthread_t));
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	struct params_check_edge *para =  calloc(1, sizeof(struct params_check_edge));
	para->g = g;
	para->i = 0;
	para->j = 0;
	para->listE = listE;
	para->n_bucks = n_bucks;
	para->avg_bin_hash = avg_bin_hash;
	para->thres_score = thres_score;
	para->out_file = out_file;
	para->n_listE = n_e;
	for (int i = 0; i < n_threads; ++i)
		pthread_create(&thr[i], &attr, process_check_edge, para);
	for (int i = 0; i < n_threads; ++i)
		pthread_join(thr[i], NULL);
	free(thr);
	pthread_attr_destroy(&attr);
	free(listE);
}

void swap(int *a, int *b)
{ 
	int tmp;
	tmp = *a; 
	*a = *b; 
	*b = tmp;
}

void normalize_min_index(struct asm_graph_t *g, struct contig_edge *e)
{
	assert(e->src < g->n_e && e->des < g->n_e);
	int rc_id_src = g->edges[e->src].rc_id;
	if (rc_id_src < e->src) {
		e->src = rc_id_src;
		e->rv_src ^= 1;
	}
	int rc_id_des = g->edges[e->des].rc_id;
	if (rc_id_des < e->des) {
		e->des = rc_id_des ;
		e->rv_des ^= 1;
	}
//	if (e->src > e->des) {
//		swap(&e->src, &e->des);
//		swap(&e->rv_src, &e->rv_des);
//		e->rv_src ^= 1 ;
//		e->rv_des ^= 1;
//	}
}

void normalize_one_dir(struct asm_graph_t *g, struct contig_edge *e)
{
	assert(e->src < g->n_e && e->des < g->n_e);
	if (e->rv_src == 1) {
		e->src = g->edges[e->src].rc_id;
		e->rv_src = 0;
	}
	if (e->rv_des == 1) {
		e->des = g->edges[e->des].rc_id;
		e->rv_des = 0;
	}
}

void add_contig_edge(struct asm_graph_t *g,struct contig_edge *listE, int pos, int src, int des, float score0)
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

void unique_edge(struct contig_edge *listE, int *n_e)
{
	int new_n_e = 0;
	for (int i = 0; i < *n_e; ) {
		int j = i;
		struct contig_edge best;
		best.score0 = 0;
		while (j < *n_e && equal_contig_edge(&listE[j], &listE[i])){
			if (better_contig_edge(&listE[j], &best)){
				best = listE[j];
			}
			++j;
		}
//		for (int jj = i; jj < j; jj++) if (listE[jj].score0 > best.score0*0.8) {
//			listE[new_n_e++] = listE[jj];
//		}
		listE[new_n_e++] = best;
		i = j;
	}
	*n_e = new_n_e;
}

int *binary_search(int *list, int n, int value)
{
	int l = 0, r = n - 1;
	while (l != r) {
		int m = (l + r) / 2;
		if (value > list[m]) 
			l = m + 1;
		else
			r = m;
	}
	assert(list[l] == value);
	if (list[l] != value)
		return list-1;
	return &list[l];
}

void unique_vertex(int *listV, int *n_v)
{
	int new_n_v = 0;
	for (int i = 0; i < *n_v; ) {
		int j = i;
		while (j < *n_v && listV[i] == listV[j]) j++;
		listV[new_n_v++] = listV[i];
		i = j;
	}
	*n_v = new_n_v;
}

void build_V_from_E(struct contig_edge *listE, int n_e, int **listV, int *n_v)
{
	*listV = NULL; 
	*n_v = 0;
	for (int i = 0; i < n_e; i++) {
		__VERBOSE_FLAG(log_build_V, "%d %d BBBB\n", listE[i].src, listE[i].des);

		++(*n_v);
		int t = (*n_v)*sizeof(int);
		*listV = realloc(*listV, t);
		(*listV)[(*n_v)-1] = listE[i].src;

		++(*n_v);
		*listV = realloc(*listV, (*n_v)*(sizeof(int)));
		(*listV)[(*n_v)-1] = listE[i].des;
	}
	qsort(*listV, *n_v, sizeof(int), less_uint32);
	for (int i = 0; i < *n_v; i++) {
		__VERBOSE_FLAG(log_build_V, "%d ", (*listV)[i]);
	}
	unique_vertex(*listV, n_v);
//	__VERBOSE_FLAG("listV %d\n", listV);
}

void print_seq(FILE *fp, int index, char *seq, int len, int cov)
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

gint_t dump_edge_seq_reduce_N(char **seq, uint32_t *m_seq, struct asm_edge_t *e)
{
	for (uint32_t i = 0; i < e->n_holes; i++) {
		if (e->l_holes[i] > 1000)
			e->l_holes[i] = 1000;
	}
	return dump_edge_seq(seq, m_seq, e);
}

void algo_find_hamiltonian(FILE *out_file, struct asm_graph_t *g, float *E, int n_v, int *res, int *n_res, int *listV, float avg_bin_hash)
{
	void print_contig(int index, int n_contig, int *list_contig)
	{
		char *seq = NULL, *total_seq = NULL, *NNN = NULL;
		const int len_NNN = 300;
		NNN = calloc(len_NNN, sizeof(char));
		for (int i = 0; i < len_NNN; i++) 
			NNN[i] = 'N';
		uint32_t seq_len = 0;
 		int total_len = 0;
		for(int i = 0; i < n_contig; i++) {
			int e = list_contig[i];
			int len_of_contig = dump_edge_seq_reduce_N(&seq, &seq_len, &g->edges[e]);
			total_seq = realloc(total_seq, (total_len + len_of_contig + len_NNN) * sizeof(char));
			memcpy(total_seq + total_len, seq, len_of_contig);
			memcpy(total_seq + total_len + len_of_contig, NNN, len_NNN);
			total_len += len_of_contig + len_NNN;
		}
		print_seq(out_file, index, total_seq, total_len, 1);
		for(int i = 0; i < n_contig; i++) {
			int e = list_contig[i];
		}
		free(seq);
		free(total_seq);
		free(NNN);
	}

	// deprecated due to low performance
	void insert_short_contigs(int n_big_contigs, int *big_contigs, int n_insert, int *arr_insert, 
			int n_short, int *arr_i_short, int *mark_short)
	{
		//todo @huu insert to the beginning
		__VERBOSE_FLAG(log_insert_small_contig, "big contig length %d:", n_big_contigs);
		for (int i = 0; i < n_big_contigs; i++){
			__VERBOSE_FLAG(log_insert_small_contig, "%d ", big_contigs[i]);
		}
		for (int i = 1; i < n_insert - 1; i++) {
			__VERBOSE_FLAG(log_insert_small_contig, "i insert: %d\n", i);
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
			__VERBOSE_FLAG(log_insert_small_contig, "max score: %d \n",max_score);
			//todo @huu check if normal score large enough
			if (max_score > 15) {
				arr_insert[i] = arr_i_short[pos];
				mark_short[pos] = 1;
			}
		}
		//todo @huu insert to the end
	}

	void add_insert_except_m1(int *n_arr, int **new_arr, int ele)
	{
		if (ele == -1)
			return;
		*new_arr = realloc(*new_arr, (*n_arr+1) * sizeof(int));
		(*new_arr)[*n_arr] = ele;
		++(*n_arr);
	}

	void merge_big_and_small(int *best_n_res, int **best_res, int n_insert, int *arr_insert)
	{
		int *new_arr = NULL;
		int n_arr = 0;
		add_insert_except_m1(&n_arr, &new_arr, arr_insert[0]);
		for (int i = 0; i < *best_n_res; i++) {
			add_insert_except_m1(&n_arr, &new_arr, (*best_res)[i]); 
			add_insert_except_m1(&n_arr, &new_arr, arr_insert[i]); 
		}
		free(*best_res);
		*best_n_res = n_arr;
		*best_res = new_arr;
		__VERBOSE_FLAG(log_insert_small_contig, "after merge contig length %d:", *best_n_res);
		for (int i = 0; i < *best_n_res; i++){
			__VERBOSE_FLAG(log_insert_small_contig, "%d ", (*best_res)[i]);
		}
	}

	void find_longest_path_from_node(int x, int *best_n_hamiltonian_path, int *best_hamiltonian_path, int *remain_unvisited)
	{
		__VERBOSE_FLAG(log_hamiltonian, "find longest path from node");
		remain_unvisited[x]--;
		remain_unvisited[x^1]--;
		assert(*best_n_hamiltonian_path == 0);
		*best_n_hamiltonian_path = 1;
		best_hamiltonian_path[0] = x;
		while (1) {
			int *list_adj = NULL, count_adj = 0;
			int last_pos = best_hamiltonian_path[*best_n_hamiltonian_path-1];
			for (int i = 0; i < n_v; i++) if (E[last_pos * n_v + i] != 0 && remain_unvisited[i]){
				assert(remain_unvisited[i] >0);
				list_adj = realloc(list_adj, (count_adj+1) * sizeof(int));
				list_adj[count_adj] = i;
				count_adj++;
			}
			if (count_adj == 0)
				break;
			int best_n_local_path = 0 , *best_local_path = calloc(n_v, sizeof(int)), best_add_len = 0;
			for (int i_adj = 0; i_adj < count_adj; i_adj++) {
				int adj = list_adj[i_adj];
				__VERBOSE_FLAG(log_hamiltonian, "dfs from %d %d %d\n", i_adj, adj, remain_unvisited[adj]);
				int *cur_path = calloc(n_v, sizeof(int));
				dfs_hamiltonian(
					adj, 1, count_adj, list_adj, n_v, E, remain_unvisited,
					cur_path, &best_add_len, best_hamiltonian_path + *best_n_hamiltonian_path, listV
				);
				free(cur_path);
			}
			assert(best_add_len != 0);
			__VERBOSE_FLAG(log_assert, "best n hamiltonian path n_v best add len %d %d %d\n", *best_n_hamiltonian_path, n_v, best_add_len);
			for (int i = 0; i < *best_n_hamiltonian_path + best_add_len ; i++) {
				__VERBOSE("%d ", best_hamiltonian_path[i]);
			}
			__VERBOSE("#\n\n");
			for (int i = 0 ; i < n_v; i++) __VERBOSE("%d ", remain_unvisited[i]);
			for (int i_path = *best_n_hamiltonian_path; i_path < *best_n_hamiltonian_path + best_add_len; i_path++){
				assert(remain_unvisited[best_hamiltonian_path[i_path]] >0);
				assert(remain_unvisited[best_hamiltonian_path[i_path]^1] >0);
				remain_unvisited[best_hamiltonian_path[i_path]]--;
				remain_unvisited[best_hamiltonian_path[i_path]^1]--;
			}
			*best_n_hamiltonian_path += best_add_len;
			assert(*best_n_hamiltonian_path <= n_v);
			free(best_local_path);
			for (int i = 0 ; i < n_v; i++) __VERBOSE("%d ", remain_unvisited[i]);
			free(list_adj);
		}
		for (int i = 0; i < *best_n_hamiltonian_path; i++){
			remain_unvisited[best_hamiltonian_path[i]]++;
			remain_unvisited[best_hamiltonian_path[i]^1]++;
		}
	}

	void dfs_find_connected_component(int root, int *remain_unvisited, int
			*connected_component)
	{
		connected_component[root] = 1;
		for (int i = 0; i < n_v; i++) if (remain_unvisited[i] && connected_component[i] == 0 
				&& ((E[root * n_v + i]) || (E[i * n_v + root]))) {
			dfs_find_connected_component(i, remain_unvisited, connected_component);
		}
	}

	void iter_find_longest_path(int *remain_unvisited, int *best_n_hamiltonian_path, int *best_hamiltonian_path)
	{
		int *save_remain_unvisited = calloc(n_v, sizeof(int));
		__VERBOSE("begin of iter find longest");
		for(int i = 0; i< n_v; i++) {
			__VERBOSE("%d ", remain_unvisited[i]);
		}
		__COPY_ARR(remain_unvisited, save_remain_unvisited, n_v);
		__VERBOSE_FLAG(log_hamiltonian, "iter find longest path");
		int *connected_component = calloc(n_v, sizeof(int));
		for (int i = 0; i < n_v; i++) if (remain_unvisited[i]) {
			dfs_find_connected_component(i, remain_unvisited, connected_component);
			break;
		}
		for (int i = 0; i < n_v; i++) if (remain_unvisited[i] ){ //&& connected_component[i]) {
			assert(remain_unvisited[i] > 0);
			int *hamiltonian_path = calloc(n_v, sizeof(int)), n_hamiltonian_path = 0 ;
			
			find_longest_path_from_node(i, &n_hamiltonian_path, hamiltonian_path, remain_unvisited);
			if (n_hamiltonian_path > *best_n_hamiltonian_path) {
				*best_n_hamiltonian_path = n_hamiltonian_path;
				__VERBOSE("new best hamilton iter %d %d\n", *best_n_hamiltonian_path, n_hamiltonian_path);
				for(int i_path = 0; i_path < n_hamiltonian_path; i_path++) {
					__VERBOSE("x");
					best_hamiltonian_path[i_path] = hamiltonian_path[i_path];
					__VERBOSE("%d ", best_hamiltonian_path[i_path]);
				}
			}
			free(hamiltonian_path);
		}
		__VERBOSE("best hamilton iter\n");
		for(int i = 0 ; i < *best_n_hamiltonian_path; i++) 
			__VERBOSE("%d ", best_hamiltonian_path[i]);
		for (int i = 0 ; i < n_v; i++){
			assert(save_remain_unvisited[i] ==  remain_unvisited[i]);
		}
		free(save_remain_unvisited);
	}
	
//#############################  BEGIN OF FUNCTION  ###########################
	int thres_len = global_thres_length, thres_len_min = global_thres_length_min;
	__VERBOSE_FLAG(log_hamiltonian, "xxxxxxx %d ", thres_len_min);

	int *arr_i_short = NULL, n_arr_short = 0;
	int *mark_short = NULL;
	__VERBOSE_FLAG(log_hamiltonian, "g->n_e %ld", g->n_e);
	for (int e = 0; e < g->n_e; ++e) {
		int len = get_edge_len(&g->edges[e]);
		__VERBOSE_FLAG(log_hamiltonian, "len %d\n", len);
		if (len < thres_len && len > thres_len_min) {
			++n_arr_short;
			arr_i_short = realloc(arr_i_short, n_arr_short * sizeof(int));
			mark_short = realloc(mark_short, n_arr_short * sizeof(int));
			arr_i_short[n_arr_short-1] = e;
			mark_short[n_arr_short-1] = 0;
		}
	}

	int *remain_unvisited = calloc(n_v, sizeof(int));
	for (int i = 0; i < n_v; i++) {
		float cvr = global_genome_coverage;
		__VERBOSE_FLAG(log_debug, "list v %d", listV[i]);
		float cov_times = (__get_edge_cov(&g->edges[listV[i]], g->ksize)/cvr) ;
		remain_unvisited[i] = roundint(cov_times);
		__VERBOSE_FLAG(log_hamiltonian, "%d " , remain_unvisited[i]);
	}
	int count = 0;
	while (1){
		int *best_hamiltonian_path = calloc(n_v, sizeof(int)), best_n_hamiltonian_path = 0 ;
		iter_find_longest_path(remain_unvisited, &best_n_hamiltonian_path, best_hamiltonian_path);
		__VERBOSE_FLAG(log_hamiltonian, "best_n_hamiltonian_path %d\n", best_n_hamiltonian_path);
		for(int i = 0; i < best_n_hamiltonian_path; i++) {
			__VERBOSE_FLAG(log_hamiltonian, "remain %d %d\n", best_hamiltonian_path[i], remain_unvisited[best_hamiltonian_path[i]]);
			remain_unvisited[best_hamiltonian_path[i]]--;
			remain_unvisited[best_hamiltonian_path[i]^1]--;
			assert(remain_unvisited[best_hamiltonian_path[i]]>=0);
			assert(remain_unvisited[best_hamiltonian_path[i]^1]>=0);
		}
		if (best_n_hamiltonian_path == 0){
			free(best_hamiltonian_path);
			break;
		}

		int *arr_insert = NULL, n_insert = best_n_hamiltonian_path + 1;
		arr_insert = realloc(arr_insert, n_insert * sizeof(int));
		for (int i = 0; i < n_insert; i++) arr_insert[i] = -1;
		for (int i = 0; i < best_n_hamiltonian_path; i++) 
			best_hamiltonian_path[i] = listV[best_hamiltonian_path[i]];
//		insert_short_contigs(best_n_hamiltonian_path, best_hamiltonian_path, n_insert, arr_insert, n_arr_short, arr_i_short, mark_short);
		merge_big_and_small(&best_n_hamiltonian_path, &best_hamiltonian_path, n_insert, arr_insert);

		print_contig(count, best_n_hamiltonian_path, best_hamiltonian_path);
		__VERBOSE_FLAG(log_hamiltonian, "contig path ");
		for (int i = 0; i < best_n_hamiltonian_path; i++) 
			__VERBOSE_FLAG(log_hamiltonian, "%d " , best_hamiltonian_path[i]);
		__VERBOSE_FLAG(log_hamiltonian, "best n %d\n", best_n_hamiltonian_path);
		count++;
		free(best_hamiltonian_path);
		free(arr_insert);
	}

	__VERBOSE_FLAG(log_hamiltonian, "n arr short %d\n", n_arr_short);
	for (int i = 0; i < n_arr_short; i++) if (remain_unvisited[i] == 0) {
		__VERBOSE_FLAG(log_hamiltonian, "printf short edge \n");
		int e = arr_i_short[i]; 
		uint32_t seq_len = 0;
		char *seq = NULL;
		int len = dump_edge_seq_reduce_N(&seq, &seq_len, &g->edges[e]);
		print_seq(out_file, count, seq, seq_len, 1); 
		count++;
		free(seq);
	}

	for (int e = 0; e < g->n_e; e++){
		int len  = get_edge_len(&g->edges[e]);
		__VERBOSE_FLAG(log_hamiltonian, "len very short %d %d\n", len, thres_len_min); 
		if (len < thres_len_min && len > 1000) {
			__VERBOSE_FLAG(log_hamiltonian, "printf very short edge \n");
			uint32_t seq_len = 0;
			char *seq = NULL;
			int len = dump_edge_seq_reduce_N(&seq, &seq_len, &g->edges[e]);
			print_seq(out_file, count, seq, seq_len, 1); 
			count++;
			free(seq);
		}
	}

	free(remain_unvisited);
	free(arr_i_short);
	free(mark_short);
}

void find_hamiltonian_contig_edge(FILE *out_file, struct asm_graph_t *g, struct contig_edge *listE_ori, int n_e, int n_v, int *listV, float avg_bin_hash)
{
	struct contig_edge *list_one_dir_E = calloc(2*n_e, sizeof(struct contig_edge));
	for (int i = 0; i < n_e; i++) {
		list_one_dir_E[i] = listE_ori[i];
		normalize_one_dir(g, list_one_dir_E+i);
	}
	__VERBOSE_FLAG(log_hamiltonian, "in graph when find hamiltonian n_e n_v: %d %d\n", n_e, n_v);
	__VERBOSE_FLAG(log_hamiltonian, "listV:\n");
	for (int i = 0; i < n_v; i++) {
		__VERBOSE_FLAG(log_hamiltonian, "%d ", listV[i]);
	}

	int m = 0;
	float *E = calloc(n_v * n_v, sizeof(float));
	for (int i = 0; i < n_e; i++) {
		int u = list_one_dir_E[i].src;
		u = binary_search(listV, n_v, u) - listV;
		int v = list_one_dir_E[i].des;
		v = binary_search(listV, n_v, v) - listV;
		__VERBOSE_FLAG(log_hamiltonian, "edge: %d %d \n", listV[u], listV[v]);
		if (u == -1 || v == -1){
			__VERBOSE_FLAG(log_hamiltonian, "%d %d ERRRRRR\n", list_one_dir_E[i].src, list_one_dir_E[i].des);
		}
		E[u * n_v + v] = list_one_dir_E[i].score0;
	}
	int *res = calloc(n_v, sizeof(int)), n_res=0;
	__VERBOSE_FLAG(log_hamiltonian, "n_v %d \n", n_v);
	for(int i = 0; i < n_v ; i++) {
		__VERBOSE_FLAG(log_hamiltonian, "list v i %d ", listV[i]);
	}
	algo_find_hamiltonian(out_file, g,E, n_v, res, &n_res, listV, avg_bin_hash);
	free(E);
	free(list_one_dir_E);
	free(res);
}

void print_gfa_from_E(struct asm_graph_t *g, int n_e, struct contig_edge *listE, int n_v, int *listV, FILE *out_graph)
{
	struct contig_edge *list_one_dir_E = calloc(n_e, sizeof(struct contig_edge));
	for (int i = 0; i < n_e; i++) {
		list_one_dir_E[i] = listE[i];
		normalize_min_index(g, list_one_dir_E+i);
	}
	for (int i = 0; i < n_v; i++) {
		struct asm_edge_t *e = &g->edges[listV[i]];
		char *seq = NULL;
		uint32_t seq_len = 0;
		dump_edge_seq_reduce_N(&seq, &seq_len, e);
		fprintf(out_graph,"S\t%d\t%s\tKC:i:%lu\n", listV[i], seq, e->count);
	}

	for (int i = 0; i < n_e; i++) {
		fprintf(out_graph, "L\t%d\t%c\t%d\t%c\t45M\n", 
			list_one_dir_E[i].src, list_one_dir_E[i].rv_src == 0?'+':'-', list_one_dir_E[i].des, list_one_dir_E[i].rv_des == 0?'+':'-');
		fprintf(out_graph, "L\t%d\t%c\t%d\t%c\t45M\n", 
			list_one_dir_E[i].des, list_one_dir_E[i].rv_des == 0?'-':'+', list_one_dir_E[i].src, list_one_dir_E[i].rv_src == 0?'-':'+');
	}
}

void connect_contig(FILE *fp, FILE *out_file, FILE *out_graph, struct asm_graph_t *g, float huu_1_score)
{
	init_global_params(g, huu_1_score);
	check_global_params(g);
	int n_v, *listV = NULL;
	fscanf(fp, "n_v: %d\n", &n_v);
	listV = realloc(listV , n_v * sizeof(int));
	for (int i = 0; i < n_v; i++) {
		fscanf(fp, "vertex:%d\n", &listV[i]);
	}
	float score;
	int src, des;
	struct contig_edge *listE = NULL;
	int n_e=0;
	float avg_bin_hash = get_avg_bin_hash(g);
	while (fscanf(fp,"score: %f edge: %d %d\n", &score, &src, &des) !=EOF) {
		n_e += 1;
		listE = realloc(listE, n_e*sizeof(struct contig_edge));
		add_contig_edge(g, listE, n_e-1, src, des, score);
	}
	qsort(listE, n_e, sizeof(struct contig_edge), less_contig_edge);
	__VERBOSE_FLAG(log_graph, "after sort");
	for (int i = 0; i < n_e; i++){
		__VERBOSE_FLAG(log_graph, "%d %d %d %d\n",listE[i].src, listE[i].des, listE[i].rv_src, listE[i].rv_des);
	}
	unique_edge(listE, &n_e);
	print_gfa_from_E(g, n_e, listE, n_v, listV, out_graph);
	__VERBOSE_FLAG(log_graph, "after unique");
	for (int i = 0; i < n_e; i++){
		__VERBOSE_FLAG(log_graph, "%d %d %d %d\n",listE[i].src, listE[i].des, listE[i].rv_src, listE[i].rv_des);
	}

	find_hamiltonian_contig_edge(out_file, g, listE, n_e, n_v, listV, avg_bin_hash);
	free(listE);
	free(listV);
}

