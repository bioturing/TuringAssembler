#include <stdlib.h>
#include "assembly_graph.h"
#include "k31hash.h"
#include "pthread.h"
#include "verbose.h"
#include <string.h>
#include <assert.h>
#include "math.h"
#include "attribute.h"
#include "utils.h"
#include "scaffolding/compare.h"
#include "scaffolding/score.h"
#include "scaffolding/global_params.h"
#include "scaffolding/bin_hash.h"
#include "scaffolding/edge.h"
#include "scaffolding/buck.h"
#include "scaffolding/contig.h"
#include "scaffolding/algorithm.h"
#include "scaffolding/output.h"

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
		for (int j = 0; j < MIN(20, n_bucks - 4); j+=5){
			for (int j1 = j + 4; j1 < n_bucks -4; j1++) {
				float tmp  = get_score_multiple_buck(pa->g, e, &e->bucks[j], &e->bucks[j1]);
				fprintf(f, "distance %d %f\n", j1 - j, tmp);
			}
		}
		fclose(f);
		pthread_mutex_lock(&lock_id);
		pthread_mutex_unlock(&lock_id);
	} while (1);
	pthread_exit(NULL);
	return NULL;
}

void check_contig(struct asm_graph_t *g, float avg_bin_hash) 
{
	int cmp(const void *i, const void *j)
	{
		int x = *(int *)i;
		int y = *(int *)j;
		return get_edge_len(&g->edges[x]) > get_edge_len(&g->edges[y]);
	}
	init_global_params(g);
	check_global_params(g);
	VERBOSE_FLAG(1, "check contig\n");
	FILE *f = fopen("list_edge.txt", "r");
	int n_bucks = global_n_buck;
	for (int i =0 ; i < 100; i++) {
		int x, y;
		fscanf(f, "%d %d\n", &x,&y);
		VERBOSE_FLAG(0, "x,y : %d %d\n", x, y);
		get_score_edges_matrix(g, x, y, n_bucks, avg_bin_hash);
	}
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
		pthread_mutex_unlock(&lock_id);
		int e0 = listE[i];
		int e1 = listE[j];
		VERBOSE_FLAG(3, "edge: %d %d\n", e0, e1); 
		VERBOSE_FLAG(3, "edges: %d %d\n", e0, e1); 
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
		VERBOSE_FLAG(3, "score %f\n", score.score); 
	} while (1); 
	pthread_exit(NULL);
	return NULL;
}

void build_list_contig(struct asm_graph_t *g, FILE *out_file, struct opt_build_t *opt) 
{
	init_global_params(g);
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

	float avg_bin_hash = get_avg_unique_bin_hash(g);
	VERBOSE_FLAG(1, "avg_bin_hash %f", avg_bin_hash);
	VERBOSE_FLAG(1, "n_e: %d\n", n_e);

	int n_threads = opt->n_threads;
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


void find_hamiltonian_contig_edge(FILE *out_file, struct asm_graph_t *g, struct contig_edge *listE_ori, int n_e, int n_v, int *listV, float avg_bin_hash)
{
	struct contig_edge *list_one_dir_E = calloc(2*n_e, sizeof(struct contig_edge));
	for (int i = 0; i < n_e; i++) {
		list_one_dir_E[i] = listE_ori[i];
		normalize_one_dir(g, list_one_dir_E+i);
	}
	VERBOSE_FLAG(1, "in graph when find hamiltonian n_e n_v: %d %d\n", n_e, n_v);
	VERBOSE_FLAG(3, "listV:\n");
	for (int i = 0; i < n_v; i++) {
		VERBOSE_FLAG(3, "%d ", listV[i]);
	}

	int m = 0;
	float *E = calloc(n_v * n_v, sizeof(float));
	for (int i = 0; i < n_e; i++) {
		int u = list_one_dir_E[i].src;
		u = binary_search(listV, n_v, u) - listV;
		int v = list_one_dir_E[i].des;
		v = binary_search(listV, n_v, v) - listV;
		VERBOSE_FLAG(3, "edge: %d %d \n", listV[u], listV[v]);
		if (u == -1 || v == -1){
			VERBOSE_FLAG(3, "%d %d ERRRRRR\n", list_one_dir_E[i].src, list_one_dir_E[i].des);
		}
		E[u * n_v + v] = list_one_dir_E[i].score0;
	}
	int *res = calloc(n_v, sizeof(int)), n_res=0;
	abc();
	algo_find_hamiltonian(out_file, g, E, n_v, res, &n_res, listV, avg_bin_hash);
//	algo_find_hamiltonian(FILE *out_file, struct asm_graph_t *g, float *E, int n_v, int *res, int *n_res, int *listV, float avg_bin_hash)
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

void connect_contig(FILE *fp, FILE *out_file, FILE *out_graph, struct asm_graph_t *g)
{
	init_global_params(g);
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
	float avg_bin_hash = get_avg_unique_bin_hash(g);
	while (fscanf(fp,"score: %f edge: %d %d\n", &score, &src, &des) !=EOF) {
		n_e += 1;
		listE = realloc(listE, n_e*sizeof(struct contig_edge));
		add_contig_edge(g, listE, n_e-1, src, des, score);
	}
	qsort(listE, n_e, sizeof(struct contig_edge), less_contig_edge);
	unique_edge(listE, &n_e);
	print_gfa_from_E(g, n_e, listE, n_v, listV, out_graph);

	find_hamiltonian_contig_edge(out_file, g, listE, n_e, n_v, listV, avg_bin_hash);
	free(listE);
	free(listV);
}

