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
#include "scaffolding/khash.h"

static pthread_mutex_t lock_merge = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t lock_id = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t lock_table = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t lock_append_edges = PTHREAD_MUTEX_INITIALIZER;
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

struct params_check_edge {
	struct asm_graph_t *g;
	int i, j;
	int *list_contig, n_bucks, n_list_contig;
	int n_candidate_edges;
	struct candidate_edge *list_candidate_edges;
	float avg_bin_hash, thres_score;
	FILE *out_file;
};

void *process_check_edge(void *data)
{
	struct params_check_edge *pa_check_edge = (struct params_check_edge *) data;
	struct candidate_edge *list_candidate_edges = pa_check_edge->list_candidate_edges;
	struct asm_graph_t *g = pa_check_edge->g;
	int n_bucks = pa_check_edge->n_bucks;
	float avg_bin_hash = pa_check_edge->avg_bin_hash;
	float thres_score = pa_check_edge->thres_score;
	FILE* out_file = pa_check_edge->out_file;
	int n_candidate_edges = pa_check_edge->n_candidate_edges;
	do {
		int i;
		pthread_mutex_lock(&lock_id);
		if (pa_check_edge->i == n_candidate_edges) {
			pthread_mutex_unlock(&lock_id);
			break;
		} else {
			i = pa_check_edge->i++;
		}
		pthread_mutex_unlock(&lock_id);
		int e0 = list_candidate_edges[i].src;
		int e1 = list_candidate_edges[i].des;
		VERBOSE_FLAG(0, "edge: %d %d\n", e0, e1); 
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

void get_long_contig(struct asm_graph_t *g, int min_length, int *n_contig, int **list_contig)
{
	for (int i = 0; i < g->n_e; ++i) {
		int len = get_edge_len(&g->edges[i]);
		if (len > min_length) {
			++*n_contig;
			assert(*n_contig > 0);
			*list_contig = realloc(*list_contig, (*n_contig)*sizeof(int));
			(*list_contig)[*n_contig-1] = i; 
		}
	}
}

void print_list_long_contig(FILE *f, int n_long_contig, int *list_long_contig)
{
	fprintf(f, "n_long_contig: %d\n", n_long_contig);
	for (int i = 0; i < n_long_contig ; i++) {
		int e = list_long_contig[i];
		fprintf(f, "vertex:%d\n", e);
	}

}

pthread_attr_t init_thread_attr()
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	return attr;
}

struct list_position{
	int n_pos;
	int *i_contig, *i_bin, *count;
} ;

KHASH_MAP_INIT_INT64(big_table, struct list_position*)

struct params_build_candidate_edges{
	struct asm_graph_t *g;
	int i;
	int n_candidate_edges;
	struct candidate_edge *list_candidate_edges;
	int n_long_contig, *list_long_contig;
	khash_t(big_table) *big_table;
};

void find_local_nearby_contig(int e, struct params_build_candidate_edges *params, int *n_local_edges, 
		struct candidate_edge **list_local_edges)
{
	for (int i = 0; i < params->n_long_contig; i++){
		VERBOSE_FLAG(0, "%d %ld\n ", i, (*n_local_edges+1) * sizeof(struct candidate_edge));

		VERBOSE_FLAG(0, "a\n");
		*list_local_edges = realloc(*list_local_edges, (*n_local_edges+1) *
				sizeof(struct candidate_edge));
		VERBOSE_FLAG(0, "b\n");
		struct candidate_edge *new_candidate_edge = calloc(1, sizeof(struct candidate_edge));
		new_candidate_edge->src = e;
		new_candidate_edge->des = params->list_long_contig[i];

		(*list_local_edges)[*n_local_edges] = *new_candidate_edge;
		++*n_local_edges;
	}
}

struct params_build_big_table {
	int i,j;
	struct asm_graph_t *g;
	khash_t(big_table) *big_table;
};

void *process_build_big_table(void *data)
{
	void next_index(struct params_build_big_table *params, struct asm_graph_t *g)
	{
		int len = 0;
		struct asm_edge_t *e = &g->edges[params->i];
		int n_bucks = (get_edge_len(e) + 1000-1) / 1000;
		(params->j)++;
		if (params->j == n_bucks) {
			params->i++;
			params->j = 0;
		}
		if (params->j == 3)
			params->j = MAX(n_bucks - 3, 3);
	}
	
	// ________________________________________BEGIN___________________________________
	struct params_build_big_table *params = (struct params_build_big_table *) data;
	struct asm_graph_t *g = params->g;
	khash_t(big_table) *big_table = params->big_table;

	do {
		pthread_mutex_lock(&lock_id);
		next_index(params, g);
		if (params->i == g->n_e){
			pthread_mutex_unlock(&lock_id);
			break;
		}
		int i = params->i;
		int j = params->j;
		int new_i_contig = i;
		int new_i_bin = j;
		int new_count;
		pthread_mutex_unlock(&lock_id);
		
		struct asm_edge_t *e = &g->edges[i];
		struct barcode_hash_t *buck = &e->bucks[j];
		for (int l = 0; l < buck->size; l++){
			if (buck->cnts[l] > 0){
				new_count = buck->cnts[l];
				uint64_t barcode  = buck->keys[l];
				pthread_mutex_lock(&lock_table);
				khint_t k = kh_get(big_table, big_table, barcode);
				if (k == kh_end(big_table)) {
					int tmp;
					k = kh_put(big_table, big_table, barcode, &tmp);
					VERBOSE_FLAG(3, "tmp is %d\n", tmp);
					assert(tmp == 1);
				}
				assert(k != kh_end(big_table));
				struct list_position *pos = kh_value(big_table, k);
				if (pos == NULL){
					pos = calloc(1, sizeof(struct list_position));
					kh_value(big_table,k) = pos;
				} 
				int v = pos->n_pos;
				pos->i_contig = realloc(pos->i_contig, (v+1) *
						sizeof(int));
				pos->i_contig[v] = new_i_contig;
				pos->i_bin = realloc(pos->i_bin, (v+1) *
						sizeof(int));
				pos->i_bin[v] = new_i_bin;
				pos->count = realloc(pos->count, (v+1) * sizeof(int));
				pos->count[v] = new_count;
				pos->n_pos++;
				pthread_mutex_unlock(&lock_table);
			}
   		}
		// todo @huu lock at only this entry pthread_mutex_lock(
	} while (1);
}

khash_t(big_table) *build_big_table(struct asm_graph_t *g, struct opt_proc_t *opt)
{
	pthread_t *thr = (pthread_t *)calloc(opt->n_threads, sizeof(pthread_t));
	pthread_attr_t attr = init_thread_attr();

	struct params_build_big_table *params_build_table = calloc(1, sizeof(struct
				params_build_big_table));
	params_build_table->i = 0;
	params_build_table->j = -1;
	params_build_table->g = g;
	params_build_table->big_table = kh_init(big_table);
	// todo @huu auto resize
	kh_resize(big_table, params_build_table->big_table, 100000000);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(&thr[i], &attr, process_build_big_table, params_build_table);
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(thr[i], NULL);

	VERBOSE_FLAG(1, "build done");
	khash_t(big_table) *big_table = params_build_table->big_table ;

	//test
	struct barcode_hash_t buck = g->edges[152].bucks[0];
	int *count = calloc(2000, 4);
	for (int jj = 0; jj < buck.size; jj++) {
		uint64_t barcode;
			if (buck.cnts[jj] > 0) {
				barcode = buck.keys[jj];
			}else {
				continue;
			}
		khint_t k = kh_get(big_table, big_table, barcode);
		assert(k != kh_end(big_table));
		struct list_position *pos = kh_value(big_table, k);
		assert(pos != NULL);
		for (int i = 0; i < pos->n_pos; i++) {
			count[pos->i_contig[i]]+= pos->count[i];
			VERBOSE_FLAG(3, "%d %d", pos->i_contig[i], pos->i_bin[i]);
			VERBOSE_FLAG(3, "%d ", pos->count[i]);
		}
		VERBOSE_FLAG(3, "\n");
	}
	int sum = 0;
	for (int i = 0; i < 1000; i++) {
		VERBOSE_FLAG(3, "%d %d \n", i, count[i]);
		sum += count[i];
	};
	VERBOSE_FLAG(1, "sum %d\n", sum);

	//end test
	return big_table;
}

void *process_build_candidate_edges(void *data)
{
	struct params_build_candidate_edges *params_candidate = 
		(struct params_build_candidate_edges*) data;
	int *list_long_contig = params_candidate->list_long_contig;
	do {
		int i;
		pthread_mutex_lock(&lock_id);
		if (params_candidate->i == params_candidate->n_long_contig) {
			pthread_mutex_unlock(&lock_id);
			break;
		} else {
			i = params_candidate->i++;
		}
		pthread_mutex_unlock(&lock_id);
		int e = list_long_contig[i];
		int n_local_edges = 0; 
		struct candidate_edge *list_local_edges = NULL; 
		find_local_nearby_contig(e, params_candidate,
				&n_local_edges,
				&list_local_edges);
		pthread_mutex_lock(&lock_append_edges);
		int n = params_candidate->n_candidate_edges;
		params_candidate->list_candidate_edges =
			realloc(params_candidate->list_candidate_edges, (n + n_local_edges) *
					sizeof(struct candidate_edge));
		COPY_ARR(list_local_edges, params_candidate->list_candidate_edges+n,
				n_local_edges);
		params_candidate->n_candidate_edges += n_local_edges; 
		VERBOSE_FLAG(1, "n local edge %d", n_local_edges);
		pthread_mutex_unlock(&lock_append_edges);
	} while(1);
	pthread_exit(NULL);
}

void build_list_edges(struct asm_graph_t *g, FILE *out_file, struct opt_proc_t *opt) 
{
	init_global_params(g);
	check_global_params(g);
	int n_long_contig = 0, *list_long_contig = NULL;
	get_long_contig(g, global_thres_length, &n_long_contig, &list_long_contig);
	print_list_long_contig(out_file, n_long_contig, list_long_contig);
	float avg_bin_hash = get_avg_unique_bin_hash(g);

	VERBOSE_FLAG(1, "avg_bin_hash %f", avg_bin_hash);
	VERBOSE_FLAG(1, "n_long_contig: %d\n", n_long_contig);

	pthread_t *thr = (pthread_t *)calloc(opt->n_threads, sizeof(pthread_t));
	pthread_attr_t attr = init_thread_attr();

	struct contig_edge *list_candidate_edges = NULL;
	struct params_build_candidate_edges *params_candidate = calloc(1,
			sizeof(struct params_build_candidate_edges));
	assert(list_long_contig != NULL);
	params_candidate->g = g;
	params_candidate->i = 0;
	params_candidate->big_table = build_big_table(g, opt);
	params_candidate->n_long_contig = n_long_contig;
	params_candidate->list_long_contig = list_long_contig;
	params_candidate->n_candidate_edges = 0;
	params_candidate->list_candidate_edges = NULL;

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(&thr[i], &attr, process_build_candidate_edges, params_candidate);
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(thr[i], NULL);

	struct params_check_edge *para =  calloc(1, sizeof(struct params_check_edge));
	para->g = g;
	para->i = 0;
	para->n_bucks = global_n_buck;
	para->avg_bin_hash = avg_bin_hash;
	para->thres_score = global_thres_bucks_score;
	para->out_file = out_file;
	para->list_candidate_edges = params_candidate->list_candidate_edges;
	para->n_candidate_edges = params_candidate->n_candidate_edges;
	VERBOSE_FLAG(1, "n candidate edges %d", para->n_candidate_edges);
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(&thr[i], &attr, process_check_edge, para);
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(thr[i], NULL);
	free(thr);
	pthread_attr_destroy(&attr);
	free(list_long_contig);
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
	free(E);
	free(list_one_dir_E);
	free(res);
}

void connect_contig(FILE *fp, FILE *out_file, FILE *out_graph, struct asm_graph_t *g)
{
	init_global_params(g);
	check_global_params(g);
	int n_v, *listV = NULL;
	fscanf(fp, "n_v: %d\n", &n_v);
	listV = realloc(listV , n_v * sizeof(int)); for (int i = 0; i < n_v; i++) {
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

void scaffolding_test(struct asm_graph_t *g, struct opt_proc_t *opt)
{
	init_global_params(g);
	check_global_params(g);
	build_big_table(g, opt);
}

