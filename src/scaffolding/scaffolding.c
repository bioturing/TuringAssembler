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
static pthread_mutex_t lock_put_table = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t lock_append_edges = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t lock_write_file = PTHREAD_MUTEX_INITIALIZER;


struct params{
	struct opt_proc_t *opt;
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
				float tmp  = get_score_multiple_buck(pa->g, e, &e->bucks[j], &e->bucks[j1], pa->opt);
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
	float avg_bin_hash, thres_score ;
	float *list_candidate_scores;
	FILE *out_file;
	struct opt_proc_t *opt;
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
		int i_candidate;
		pthread_mutex_lock(&lock_id);
		if (pa_check_edge->i == n_candidate_edges) {
			pthread_mutex_unlock(&lock_id);
			break;
		} else {
			i_candidate = pa_check_edge->i++;
		}
		pthread_mutex_unlock(&lock_id);
		int e0 = list_candidate_edges[i_candidate].src;
		int e1 = list_candidate_edges[i_candidate].des;
		assert(e1 < g->n_e && e0 < g->n_e);
		struct bucks_score score = get_score_edges_res(e0, e1, g, n_bucks, 
						avg_bin_hash, pa_check_edge->opt);
		assert(!isnan(score.score));
		int check = 0;
		if (e0 == e1){
			// todo @huu what about metagenomics
			float cvr = global_genome_coverage;
			if (__get_edge_cov(&g->edges[e0], g->ksize)/cvr > 1.8) {
				check = 1;
			}
		} else if (e0 == g->edges[e1].rc_id) {
		} else {
			check = 1;
		}
		if (check) {
			pa_check_edge->list_candidate_scores[i_candidate] = score.score;
		} else {
			pa_check_edge->list_candidate_scores[i_candidate] = -1;
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
	pthread_mutex_t lock_entry ;
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

void count_pos(int *count, struct list_position *pos)
{
	assert(pos != NULL && count != NULL);
	VERBOSE_FLAG(3, "n pos %d\n", pos->n_pos);
	for (int i = 0; i < pos->n_pos; i++){
		count[pos->i_contig[i]] += pos->count[i];
	}
}

void find_local_nearby_contig(int i_edge, struct params_build_candidate_edges *params, int *n_local_edges, 
		struct candidate_edge **list_local_edges)
{
	int *count = calloc(params->g->n_e, sizeof(int));
	struct asm_edge_t *e = &params->g->edges[i_edge];
	khash_t(big_table) *big_table = params->big_table;
	int n_bucks = (get_edge_len(e) + 1000-1) / 1000;
	for (int i = n_bucks - 3; i <  n_bucks -1; i++) {
		struct barcode_hash_t *buck = &e->bucks[i];
		for (int j = 0; j < buck->n_item; j++){
			if (buck->cnts[j] > 20) {
				uint64_t barcode = buck->keys[j];
				khint_t k = kh_get(big_table, big_table, barcode);
				if (k == kh_end(big_table))
					continue;
				struct list_position *pos = kh_value(big_table, k);
				count_pos(count, pos);
			}
		} 
	}

	for (int i = 0; i < params->n_long_contig; i++){
		int long_contig = params->list_long_contig[i];
		float edge_cov  = __get_edge_cov(&params->g->edges[long_contig], params->g->ksize);
		float value = count[long_contig] ;
		*list_local_edges = realloc(*list_local_edges, (*n_local_edges+1) *
				sizeof(struct candidate_edge));
		struct candidate_edge *new_candidate_edge = calloc(1, sizeof(struct candidate_edge));
		new_candidate_edge->src = i_edge;
		new_candidate_edge->des = long_contig;
		new_candidate_edge->score = value;

		(*list_local_edges)[*n_local_edges] = *new_candidate_edge;
		++*n_local_edges;
	}

	qsort(*list_local_edges, *n_local_edges, sizeof(struct candidate_edge), decending_candidate_edge);
	*n_local_edges = MIN(*n_local_edges, 50);
//	for (int i = 0; i < *n_local_edges; i++){
//		struct candidate_edge e = (*list_local_edges)[i];
//		if ((*list_local_edges)[i].score < global_filter_constant / 2 )
//		{
//			*n_local_edges = i;
//			break;
//		}
//	}
	
}

struct params_build_big_table {
	int i,j;
	struct asm_graph_t *g;
	khash_t(big_table) *big_table;
	int n_long_contig, *list_long_contig;
};

void *process_build_big_table(void *data)
{
	void next_index(struct params_build_big_table *params, struct asm_graph_t *g)
	{
		int len = 0;
		struct asm_edge_t *e = &g->edges[params->list_long_contig[params->i]];
		int n_bucks = (get_edge_len(e) + 1000-1) / 1000;
		(params->j)++;
		if (params->j == n_bucks) {
			params->i++;
			params->j = 0;
		}
		if (params->j == 3)
			params->j = n_bucks - 3;
	}
	
	// ________________________________________BEGIN___________________________________
	struct params_build_big_table *params = (struct params_build_big_table *) data;
	struct asm_graph_t *g = params->g;
	khash_t(big_table) *big_table = params->big_table;
	int *list_long_contig = params->list_long_contig;

	do {
		pthread_mutex_lock(&lock_id);
		if (params->i == params->n_long_contig){
			pthread_mutex_unlock(&lock_id);
			break;
		}
		int i = params->i;
		int j = params->j;
		int new_i_contig = params->list_long_contig[i];
		int new_i_bin = j;
		int new_count;
		next_index(params, g);
		pthread_mutex_unlock(&lock_id);
		
		struct asm_edge_t *e = &g->edges[list_long_contig[i]];
		assert(e != NULL);
		struct barcode_hash_t *buck = &e->bucks[j];
		for (int l = 0; l < buck->size; l++){
			if (buck->cnts[l] > 0){
				new_count = buck->cnts[l];
				uint64_t barcode  = buck->keys[l];
				pthread_mutex_lock(&lock_put_table);
				khint_t k = kh_get(big_table, big_table, barcode);
				if (k == kh_end(big_table)) {
					int tmp=1;
					k = kh_put(big_table, big_table, barcode, &tmp);
					VERBOSE_FLAG(3, "tmp is %d\n", tmp);
					assert(tmp == 1);
				}
				pthread_mutex_unlock(&lock_put_table);
				assert(k != kh_end(big_table));
				pthread_mutex_lock(&lock_put_table);
				struct list_position *pos = kh_value(big_table, k);
				if (pos == NULL){
					pos = calloc(1, sizeof(struct list_position));
					kh_value(big_table,k) = pos;
					pos->lock_entry = (pthread_mutex_t)PTHREAD_MUTEX_INITIALIZER ;
				} 
				pthread_mutex_unlock(&lock_put_table);
				pthread_mutex_lock(&pos->lock_entry);
				int v = pos->n_pos;
				pos->n_pos++;
				if ((v&(v-1)) == 0) {
					pos->i_contig = realloc(pos->i_contig, (v*2+1) *
							sizeof(int));
					pos->i_bin = realloc(pos->i_bin, (v*2+1) *
							sizeof(int));
					pos->count = realloc(pos->count, (v*2+1) * sizeof(int));
				}
					
				pos->i_contig[v] = new_i_contig;
				pos->i_bin[v] = new_i_bin;
				pos->count[v] = new_count;
				pthread_mutex_unlock(&pos->lock_entry);
			}
   		}
		// todo @huu lock at only this entry pthread_mutex_lock(
	} while (1);
}

khash_t(big_table) *build_big_table(struct asm_graph_t *g, struct opt_proc_t *opt, int n_long_contig, int *list_long_contig)
{
	pthread_t *thr = (pthread_t *)calloc(opt->n_threads, sizeof(pthread_t));
	pthread_attr_t attr = init_thread_attr();

	struct params_build_big_table *params_build_table = calloc(1, sizeof(struct
				params_build_big_table));
	params_build_table->i = 0;
	params_build_table->j = 0;
	params_build_table->g = g;
	params_build_table->big_table = kh_init(big_table);
	params_build_table->n_long_contig = n_long_contig;
	params_build_table->list_long_contig = list_long_contig;
	// todo @huu auto resize
	kh_resize(big_table, params_build_table->big_table, 100000000);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(&thr[i], &attr, process_build_big_table, params_build_table);
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(thr[i], NULL);

	VERBOSE_FLAG(1, "build done\n");
	khash_t(big_table) *big_table = params_build_table->big_table ;

	//test
	//end test
	return big_table;
}

void *process_build_candidate_edges(void *data)
{
	VERBOSE_FLAG(1, "build edgeee \n");
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
		VERBOSE_FLAG(1, "build edge from %d\n", i);
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
		pthread_mutex_unlock(&lock_append_edges);
	} while(1);
	pthread_exit(NULL);
}

void init_params_build_candidate_edges(struct params_build_candidate_edges **params_candidate, 
	struct asm_graph_t *g, int n_long_contig, int *list_long_contig, struct opt_proc_t *opt)
{
	(*params_candidate) = calloc(1, sizeof(struct params_build_candidate_edges));
	(*params_candidate)->g = g;
	(*params_candidate)->i = 0;
	(*params_candidate)->big_table = build_big_table(g, opt, n_long_contig, list_long_contig);
	(*params_candidate)->n_long_contig = n_long_contig;
	(*params_candidate)->list_long_contig = list_long_contig;
	(*params_candidate)->n_candidate_edges = 0;
	(*params_candidate)->list_candidate_edges = NULL;
}

void run_parallel_build_candidate_edges(struct params_build_candidate_edges *params_candidate, int n_threads)
{
	assert(params_candidate != NULL);
	pthread_t *thr = (pthread_t *)calloc(n_threads, sizeof(pthread_t));
	pthread_attr_t attr = init_thread_attr();
	for (int i = 0; i < n_threads; ++i)
		pthread_create(&thr[i], &attr, process_build_candidate_edges, params_candidate);
	for (int i = 0; i < n_threads; ++i)
		pthread_join(thr[i], NULL);
	free(thr);
	pthread_attr_destroy(&attr);
}

void init_params_check_edge(struct params_check_edge **params, struct asm_graph_t *g, 
	float avg_bin_hash, FILE *out_file, struct params_build_candidate_edges *params_candidate, 
	struct opt_proc_t *opt)
{
	*params = calloc(1, sizeof(struct params_check_edge));
	(*params)->g = g;
	(*params)->i = 0;
	(*params)->n_bucks = global_n_buck;
	(*params)->avg_bin_hash = avg_bin_hash;
	(*params)->thres_score = global_thres_bucks_score;
	(*params)->out_file = out_file;
	(*params)->n_candidate_edges = params_candidate->n_candidate_edges;
	(*params)->list_candidate_edges = params_candidate->list_candidate_edges;
	(*params)->list_candidate_scores = calloc(params_candidate->n_candidate_edges, sizeof(float));
	(*params)->opt = opt;
}

void parallel_build_edge_score(struct params_check_edge *para, int n_threads)
{
	pthread_t *thr = (pthread_t *)calloc(n_threads, sizeof(pthread_t));
	pthread_attr_t attr = init_thread_attr();
	for (int i = 0; i < n_threads; ++i)
		pthread_create(&thr[i], &attr, process_check_edge, para);
	for (int i = 0; i < n_threads; ++i)
		pthread_join(thr[i], NULL);
	free(thr);
	pthread_attr_destroy(&attr);
}

void remove_lov_cov(struct asm_graph_t *g)
{
	float cvr = global_genome_coverage;
	int count = 0;
	for (int i_e = 0; i_e < g->n_e; i_e++) {
		float edge_cov = __get_edge_cov(&g->edges[i_e], g->ksize)/cvr;
		VERBOSE_FLAG(0, "edge %d len:%d cov: %f\n", i_e , get_edge_len(&g->edges[i_e]), edge_cov);
		if (edge_cov > 0.5){
			int rc_id = g->edges[i_e].rc_id;
			g->edges[rc_id].rc_id = count;
			g->edges[count] =  g->edges[i_e];
			count++;
		}
	}
	g->n_e = count;
}

void build_list_edges(struct asm_graph_t *g, FILE *out_file, struct opt_proc_t *opt) 
{
	init_global_params(g);
	check_global_params(g);
	// todo @huu: uncomment the following line. Comment for debug purpose only
//	if (!opt->metagenomics)
//		remove_lov_cov(g);
	for (int i = 0 ; i < g->n_e; i++) {
		printf("this is i %d rc_id %d huy %d\n", i, g->edges[i].rc_id, g->edges[i].source);
	}

	float cvr = global_genome_coverage;
	for (int i_e = 0; i_e < g->n_e; i_e++) {
		float edge_cov = __get_edge_cov(&g->edges[i_e], g->ksize)/cvr;
		VERBOSE_FLAG(0, "edge %d len:%d cov: %f\n", i_e , get_edge_len(&g->edges[i_e]), edge_cov);
	}
	int n_long_contig = 0, *list_long_contig = NULL;
	get_long_contig(g, global_thres_length, &n_long_contig, &list_long_contig);
	assert(list_long_contig != NULL);
	print_list_long_contig(out_file, n_long_contig, list_long_contig);
	float avg_bin_hash = get_avg_unique_bin_hash(g);

	VERBOSE_FLAG(1, "avg_bin_hash %f\n", avg_bin_hash);
	VERBOSE_FLAG(1, "n_long_contig: %d\n", n_long_contig);

	struct params_build_candidate_edges *params_candidate ;
	init_params_build_candidate_edges(&params_candidate, g, n_long_contig, list_long_contig, opt);
	run_parallel_build_candidate_edges(params_candidate, opt->n_threads);
	VERBOSE_FLAG(1, "n candidate edges %d", params_candidate->n_candidate_edges);
	VERBOSE_FLAG(1, "done build candidate edge");

	struct params_check_edge *para;
	init_params_check_edge(&para, g, avg_bin_hash, out_file, params_candidate, opt);
	parallel_build_edge_score(para, opt->n_threads);

	VERBOSE_FLAG(1, "find real edge");
	int n_bucks = para->n_bucks;
	for (int i = 0; i < n_long_contig; i++) {
		int i_edge =  list_long_contig[i];
		struct contig_edge *list_E = NULL;
		int n_contig_edge = 0;
		for (int j = 0; j < para->n_candidate_edges; j++) 
						if (para->list_candidate_edges[j].src == i_edge) {
			int i2_edge = para->list_candidate_edges[j].des; 
			float score = para->list_candidate_scores[j];
			struct contig_edge new_edge;
			assert(!isnan(score));
			new_edge.src = i_edge;
			new_edge.des = i2_edge;
			new_edge.score0 = score;
			n_contig_edge++;
			list_E = realloc(list_E, n_contig_edge * (sizeof(struct contig_edge)));
			list_E[n_contig_edge-1] = new_edge;
		}
		qsort(list_E, n_contig_edge, sizeof(struct contig_edge), decending_edge_score);
		for (int j = 0; j < MIN(global_number_degree, n_contig_edge) ; j++) {
			int j_edge = list_E[j].des;
			float score = list_E[j].score0;
			if ((score < para->thres_score) || score <=0) {
				break;
			}
			// todo @huu should I used this function? It seem not good now.
//			if (!check_replicate_contig_edge(g, i_edge, j_edge, n_bucks, para->thres_score*(n_bucks * n_bucks - 1), avg_bin_hash)) {
				fprintf(out_file, "score: %f edge: %d %d\n", list_E[j].score0, list_E[j].src, list_E[j].des); 
//			}
		}
		free(list_E);
	}

	free(list_long_contig);
}

void find_hamiltonian_contig_edge(FILE *out_file, struct asm_graph_t *g, 
		struct contig_edge *listE_ori, int n_e, int n_v, int *listV, float avg_bin_hash)
{
	struct contig_edge *list_one_dir_E = calloc(n_e, sizeof(struct contig_edge));
	for (int i = 0; i < n_e; i++) {
		list_one_dir_E[i] = listE_ori[i];
		normalize_one_dir(g, list_one_dir_E+i);
	}
	for (int i = 0; i < n_e; i++) {
		VERBOSE_FLAG(3, "unique4 %d %d\n", list_one_dir_E[i].src, list_one_dir_E[i].des);
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
		uint32_t v_rc_id = g->edges[v].rc_id;
		uint32_t u_rc_id = g->edges[u].rc_id;
		E[v_rc_id * n_v + u_rc_id] = list_one_dir_E[i].score0;
	}
	int *res = calloc(n_v, sizeof(int)), n_res=0;
	algo_find_hamiltonian(out_file, g, E, n_v, res, &n_res, listV, avg_bin_hash);
	free(E);
	free(list_one_dir_E);
	free(res);
}

void connect_contig(FILE *fp, FILE *out_file, FILE *out_graph, struct asm_graph_t *g, struct
		opt_proc_t *opt)
{
	init_global_params(g);
	check_global_params(g);
	int n_v, *listV = NULL;
	fscanf(fp, "n_long_contig: %d\n", &n_v);
	listV = realloc(listV , n_v * sizeof(int)); for (int i = 0; i < n_v; i++) {
		fscanf(fp, "vertex:%d\n", &listV[i]);
	}
	VERBOSE_FLAG(1, "load n_long_contig %d done", n_v);
	float score;
	int src, des;
	struct contig_edge *listE = NULL;
	int n_e=0;
	float avg_bin_hash = get_avg_unique_bin_hash(g);
	while (fscanf(fp,"score: %f edge: %d %d\n", &score, &src, &des) !=EOF) {
		n_e += 1;
		listE = realloc(listE, n_e*sizeof(struct contig_edge));
		add_contig_edge(g, listE, n_e-1, src, des, score, opt);
	}
	VERBOSE_FLAG(1, "done add contig edge\n");
	qsort(listE, n_e, sizeof(struct contig_edge), less_contig_edge);
	unique_edge(listE, &n_e);
	VERBOSE_FLAG(1, "n_e after unique %d \n", n_e);

	char *tmp = str_concate(opt->out_dir, "/list_unique_contig");
	FILE *f_unique = fopen(tmp, "w");
	for (int i = 0; i < n_e; i++) {
		normalize_one_dir(g, &listE[i]);
		fprintf(f_unique, "%d %d\n", listE[i].src, listE[i].des);
	}
	fclose(f_unique);


//	print_gfa_from_E(g, n_e, listE, n_v, listV, out_graph);

	find_hamiltonian_contig_edge(out_file, g, listE, n_e, n_v, listV, avg_bin_hash);
	free(listE);
	free(listV);
}

void scaffolding_test(struct asm_graph_t *g, struct opt_proc_t *opt)
{
	init_global_params(g);
	__VERBOSE("n_e: %ld\n", g->n_e);
	for (int i = 0; i < g->n_e; i++)
		__VERBOSE("edge %d length %d\n", i, g->edges[i].seq_len);

	check_global_params(g);
	//build_big_table(g, opt);
}

