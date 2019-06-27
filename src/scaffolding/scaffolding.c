#include <stdlib.h>
#include "assembly_graph.h"
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
#include "scaffolding/scaffold.h"

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

struct params_check_edge {
	struct asm_graph_t *g;
	int i, j;
	int *list_contig, n_list_contig;
	int n_candidate_edges;
	struct candidate_edge *list_candidate_edges;
	float avg_bin_hash, thres_score ;
	FILE *out_file;
	struct opt_proc_t *opt;
};

void destroy_params_check_edge(struct params_check_edge *para)
{
	free(para->list_contig);
}

void *process_build_edge_score(void *data)
{
	struct params_check_edge *pa_check_edge = (struct params_check_edge *) data;
	struct candidate_edge *list_candidate_edges = pa_check_edge->list_candidate_edges;
	struct asm_graph_t *g = pa_check_edge->g;
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
		struct pair_contigs_score *score = get_score_edges_res(e0, e1, g, 
						avg_bin_hash, pa_check_edge->opt);
		int check = 0;
		if (e0 == e1){
			// todo @huu what about metagenomics
			check = 1;
		} else if (e0 == g->edges[e1].rc_id) {
		} else {
			check = 1;
		}
		if (check) {
			pa_check_edge->list_candidate_edges[i_candidate].score = *score;
		} else {
			struct pair_contigs_score *sc = calloc(1, sizeof(struct pair_contigs_score));
			sc->bc_score =-1;
			pa_check_edge->list_candidate_edges[i_candidate].score = *sc;
		}
		VERBOSE_FLAG(3, "score %f\n", score); 
	} while (1); 
	pthread_exit(NULL);
	return NULL;
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
	int *i_contig, *i_bin;
	pthread_mutex_t lock_entry ;
} ;

KHASH_MAP_INIT_INT64(big_table, struct list_position*)

struct params_build_candidate_edges{
	struct asm_graph_t *g;
	int i;
	int n_candidate_edges;
	struct candidate_edge *list_candidate_edges;
	khash_t(big_table) *big_table;
};

void destroy_params_build_candidate(struct params_build_candidate_edges *para)
{
	free(para->list_candidate_edges);
	free(para->big_table);
	free(para);
}

void count_pos(int *count, struct list_position *pos)
{
	assert(pos != NULL && count != NULL);
	VERBOSE_FLAG(3, "n pos %d\n", pos->n_pos);
	for (int i = 0; i < pos->n_pos; i++){
		count[pos->i_contig[i]]++ ;
	}
}

void find_local_nearby_contig(int i_edge, struct params_build_candidate_edges *params, int *n_local_edges, 
		struct candidate_edge **list_local_edges)
{
	struct asm_graph_t *g = params->g;
	int i_rev_edge = g->edges[i_edge].rc_id;
	struct asm_edge_t *e = &params->g->edges[i_edge];
	struct asm_edge_t *rev_e= &params->g->edges[i_rev_edge];
	int *count = calloc(params->g->n_e, sizeof(int));

	khash_t(big_table) *big_table = params->big_table;
	if (get_edge_len(e) < global_thres_short_len)
		return;

	struct barcode_hash_t *buck = &rev_e->barcodes;
	for (int j = 0; j < buck->size; j++){
		if (buck->keys[j] != (uint64_t)(-1)) {
			uint64_t barcode = buck->keys[j];
			khint_t k = kh_get(big_table, big_table, barcode);
			if (k == kh_end(big_table))
				continue;
			struct list_position *pos = kh_value(big_table, k);
			count_pos(count, pos);
		}
	} 

	for (int i_contig = 0; i_contig < g->n_e; i_contig++) {
		if (is_very_short_contig(&g->edges[i_contig]))
			continue;
		float edge_cov  = __get_edge_cov(&params->g->edges[i_contig], params->g->ksize);
		int value = count[i_contig] ;
		*list_local_edges = realloc(*list_local_edges, (*n_local_edges+1) *
				sizeof(struct candidate_edge));
		struct candidate_edge *new_candidate_edge = calloc(1, sizeof(struct candidate_edge));
		new_candidate_edge->src = i_edge;
		new_candidate_edge->des = i_contig;
		new_candidate_edge->score.bc_score = 1.0*value;

		(*list_local_edges)[*n_local_edges] = *new_candidate_edge;
		++*n_local_edges;
	}

	qsort(*list_local_edges, *n_local_edges, sizeof(struct candidate_edge), decending_candidate_edge);
	*n_local_edges = MIN(global_n_candidate, *n_local_edges);
	for (int i = 0; i < *n_local_edges; i++){
		struct candidate_edge e = (*list_local_edges)[i];
		if ((*list_local_edges)[i].score.bc_score < 0.001)
		{
			*n_local_edges = i;
			break;
		}
	}
}

struct params_build_big_table {
	int i;
	struct asm_graph_t *g;
	khash_t(big_table) *big_table;
};

void *process_build_big_table(void *data)
{
	void next_index(struct params_build_big_table *params, struct asm_graph_t *g)
	{
		int len = 0;
		struct asm_edge_t *e = &g->edges[params->i];
		params->i++;
		if (params->i == g->n_e) 
			return;
	}
	
	// ________________________________________BEGIN___________________________________
	struct params_build_big_table *params = (struct params_build_big_table *) data;
	struct asm_graph_t *g = params->g;
	khash_t(big_table) *big_table = params->big_table;

	do {
		pthread_mutex_lock(&lock_id);
		if (params->i == g->n_e){
			pthread_mutex_unlock(&lock_id);
			break;
		}
		int i_contig = params->i;
		int new_i_contig = i_contig;
		int new_count;
		next_index(params, g);
		pthread_mutex_unlock(&lock_id);
		
		struct asm_edge_t *e = &g->edges[i_contig];
		struct barcode_hash_t *buck = &e->barcodes;
		for (int l = 0; l < buck->size; l++){
			if (buck->keys[l] != (uint64_t)(-1)){
				uint64_t barcode = buck->keys[l];
				pthread_mutex_lock(&lock_put_table);
				khint_t k = kh_get(big_table, big_table, barcode);
				if (k == kh_end(big_table)) {
					int tmp=1;
					k = kh_put(big_table, big_table, barcode, &tmp);
					VERBOSE_FLAG(3, "tmp is %d\n", tmp);
					assert(k != kh_end(big_table) && tmp == 1);
				}
				pthread_mutex_unlock(&lock_put_table);
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
				}
				pos->i_contig[v] = new_i_contig;
				pthread_mutex_unlock(&pos->lock_entry);
			}
   		}
	} while (1);
}

khash_t(big_table) *build_big_table(struct asm_graph_t *g, struct opt_proc_t *opt)
{
	pthread_t *thr = (pthread_t *)calloc(opt->n_threads, sizeof(pthread_t));
	pthread_attr_t attr = init_thread_attr();

	struct params_build_big_table *params_build_table = calloc(1, sizeof(struct
				params_build_big_table));
	params_build_table->i = 0;
	params_build_table->g = g;
	params_build_table->big_table = kh_init(big_table);
	// todo @huu auto resize
	kh_resize(big_table, params_build_table->big_table, 100000000);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(&thr[i], &attr, process_build_big_table, params_build_table);
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(thr[i], NULL);

	VERBOSE_FLAG(1, "build done\n");
	khash_t(big_table) *big_table = params_build_table->big_table ;

	return big_table;
}

void *process_build_candidate_edges(void *data)
{
	VERBOSE_FLAG(1, "build candidate edgeee \n");
	struct params_build_candidate_edges *params_candidate;
	params_candidate = (struct params_build_candidate_edges*) data;
	struct asm_graph_t *g = params_candidate->g;
	do {
		int i_contig;
		pthread_mutex_lock(&lock_id);
		if (params_candidate->i == g->n_e) {
			pthread_mutex_unlock(&lock_id);
			break;
		} else {
			i_contig = params_candidate->i++;
		}
		VERBOSE_FLAG(1, "build edge from %d\n", i_contig);
		pthread_mutex_unlock(&lock_id);
		int n_local_edges = 0; 
		struct candidate_edge *list_local_edges = NULL; 
		if (is_very_short_contig(&g->edges[i_contig]))
			continue;
		find_local_nearby_contig(i_contig, params_candidate,
				&n_local_edges,
				&list_local_edges);
		pthread_mutex_lock(&lock_append_edges);
		int n = params_candidate->n_candidate_edges;
		params_candidate->list_candidate_edges =
			realloc(params_candidate->list_candidate_edges, (n + n_local_edges) *
					sizeof(struct candidate_edge));
		COPY_ARR(list_local_edges, params_candidate->list_candidate_edges+n,
				n_local_edges);
		VERBOSE_FLAG(0, "nlocaledge %d \n", n_local_edges);
		free(list_local_edges);
		params_candidate->n_candidate_edges += n_local_edges; 
		pthread_mutex_unlock(&lock_append_edges);
	} while(1);
	pthread_exit(NULL);
}

struct params_build_candidate_edges* new_params_build_candidate_edges(
	struct asm_graph_t *g, struct opt_proc_t *opt)
{
	struct params_build_candidate_edges *params_candidate = 
		calloc(1, sizeof(struct params_build_candidate_edges));
	params_candidate->g = g;
	params_candidate->i = 0;
	params_candidate->big_table = build_big_table(g, opt);
	params_candidate->n_candidate_edges = 0;
	params_candidate->list_candidate_edges = NULL;
	return params_candidate;
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

struct params_check_edge *new_params_check_edge(struct asm_graph_t *g, struct opt_proc_t *opt,
	struct params_build_candidate_edges *params_candidate) 
{
	struct params_check_edge *params = calloc(1, sizeof(struct params_check_edge));
	params->g = g;
	params->i = 0;
	params->avg_bin_hash = get_avg_barcode(g);
	params->thres_score = global_thres_bucks_score;
	params->n_candidate_edges = params_candidate->n_candidate_edges;
	params->list_candidate_edges = params_candidate->list_candidate_edges;
	params->opt = opt;
	return params;
}

void parallel_build_edge_score(struct params_check_edge *para, int n_threads)
{
	pthread_t *thr = (pthread_t *)calloc(n_threads, sizeof(pthread_t));
	pthread_attr_t attr = init_thread_attr();
	for (int i = 0; i < n_threads; ++i)
		pthread_create(&thr[i], &attr, process_build_edge_score, para);
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
		if (edge_cov > 0.5){
			int rc_id = g->edges[i_e].rc_id;
			g->edges[rc_id].rc_id = count;
			g->edges[count] =  g->edges[i_e];
			count++;
		}
	}
	g->n_e = count;
}

void pre_calc_score(struct asm_graph_t *g,struct opt_proc_t* opt, struct edges_score_type *edges_score)
{ 
	struct params_build_candidate_edges *params_candidate = 
		new_params_build_candidate_edges(g, opt);
	run_parallel_build_candidate_edges(params_candidate, opt->n_threads);
	// sorted candidate edge
	VERBOSE_FLAG(1, "n candidate edges %d", params_candidate->n_candidate_edges);
	VERBOSE_FLAG(1, "done build candidate edge");

	struct params_check_edge *para = new_params_check_edge(g, opt, 
			params_candidate);
	parallel_build_edge_score(para, opt->n_threads);

	VERBOSE_FLAG(1, "find real edge");
	int i_candidate_edge = 0;
	qsort(para->list_candidate_edges, para->n_candidate_edges, sizeof(struct candidate_edge),
			ascending_index_edge);
	for (int i_contig = 0; i_contig < g->n_e; i_contig++){
		VERBOSE_FLAG(0, "i contig %d\n", i_contig);
		struct scaffold_edge *list_edges = NULL;
		int size_list_edge = 0; 
		int n_edges = 0;
		while (i_candidate_edge < para->n_candidate_edges && 
			para->list_candidate_edges[i_candidate_edge].src <= i_contig) {
			int i1_contig = para->list_candidate_edges[i_candidate_edge].des; 
			struct pair_contigs_score score = para->list_candidate_edges[i_candidate_edge].score;
			VERBOSE_FLAG(0, "i candidate %d src %d des %d score %f\n", i_candidate_edge,
					para->list_candidate_edges[i_candidate_edge].src, i1_contig, score.bc_score);
			struct scaffold_edge *new_edge = new_scaffold_edge(i_contig, i1_contig, &score);
			n_edges++;
			if (n_edges > size_list_edge) {
				size_list_edge = get_new_size(size_list_edge);
				list_edges = realloc(list_edges, size_list_edge * (sizeof(struct scaffold_edge)));
			}
			list_edges[n_edges-1] = *new_edge;
			i_candidate_edge++;
		}
		qsort(list_edges, n_edges, sizeof(struct scaffold_edge), decending_edge_score);
		for (int i_edge = 0; i_edge < MIN(global_number_degree, n_edges) ; i_edge++) {
			struct pair_contigs_score *score = &list_edges[i_edge].score;
			append_edge_score(edges_score, &(list_edges[i_edge]));
			//todo huu  not hardcode
			if ((score->bc_score == 0) ) {
				break;
			}
		}
		free(list_edges);
	}
	destroy_params_build_candidate(params_candidate);
}

int get_highest_cov_contig(struct asm_graph_t *g, int *mark, int start)
{
//	float max_edge_cov = -1;
// 	int max_edge_index = -1;
//	for (int i = 0; i < g->n_e; i++) if (is_long_contig(&g->edges[i]) && mark[i] == 0) {
//		float edge_cov = __get_edge_cov(&g->edges[i], g->ksize);
//		if (edge_cov > max_edge_cov) {
//			max_edge_cov = edge_cov;
//			max_edge_index = i;
//		}
//	}
//	return max_edge_index;
	for (int i = start; i < g->n_e; i++) if (is_long_contig(&g->edges[i]) && mark[i]) {
		return i;
	}
	return -1;
}

struct pair_contigs_score *get_score_edge(struct edges_score_type *edges_score, int src, int des)
{
	// todo @huu use binary search or matrix to speed up
	for (int i = 0; i < edges_score->n_edge; i++){
		struct scaffold_edge *edge = &edges_score->list_edge[i];
		if (edge->src == src && edge->des == des)
			return &edge->score;
	}
	struct pair_contigs_score *sc = calloc(1, sizeof(struct pair_contigs_score));
	sc->bc_score = 0;
	return sc;
}

struct pair_contigs_score *get_score(struct asm_graph_t *g, struct scaffold_path *path, int start_contig, 
	struct scaffold_edge *edge, struct edges_score_type *edges_score)
{
	struct pair_contigs_score *score = calloc(1, sizeof(struct pair_contigs_score));
	int des = edge->des;
	int distance = 0;
	//todo using distance instead of /2 
	struct pair_contigs_score *i_score = get_score_edge(edges_score, path->list_i_contig[path->n_contig-1], des);
	score->m_score = i_score->m_score;
	score->m2_score = i_score->m2_score;
	score->bc_score = i_score->bc_score;
	i_score = get_score_edge(edges_score, path->list_i_contig[path->n_contig-1], g->edges[des].rc_id);
	score->bc_score += i_score->bc_score/2;
	//todo not return here
	for (int i = path->n_contig-2; i>=0; i--){
		int src = path->list_i_contig[i];
		i_score = get_score_edge(edges_score, src, des);
		score->bc_score += i_score->bc_score;
		distance += get_edge_len(&g->edges[src]);
		if (distance > 100000)
			break;
	}
	VERBOSE_FLAG(0, "scascore %d %d %f %d %d\n", path->list_i_contig[path->n_contig-1], des, score->bc_score, score->m_score, score->m2_score);
	return score;
}

int better_edge(struct pair_contigs_score *sc0, struct pair_contigs_score *sc1)
{
	//todo wtd
	if (sc0->bc_score + sc0->m_score > 1.1 * sc1->bc_score + sc1->m_score) return 1;
	if (sc0->bc_score + sc0->m_score * 1.1 < sc1->bc_score + sc1->m_score) return 0;
	if (sc0->m2_score != sc1->m2_score)
		return (sc0->m2_score > sc1->m2_score);
	return (sc0->bc_score + sc0->m_score > sc1->bc_score + sc1->m_score);
}

void find_best_edge(struct asm_graph_t *g, struct edges_score_type *edges_score, int *start_contig,
		struct scaffold_path *path, int *found, int *mark)
{
	int n_edge_adj = 0;
	struct scaffold_edge *list_edge_adj = NULL;
	find_edge_from(edges_score, *start_contig, &n_edge_adj, &list_edge_adj);
	if (n_edge_adj == 0){
		*found = 0;
		return;
	}
	struct pair_contigs_score *max_score = calloc(1, sizeof(struct pair_contigs_score));
	max_score->bc_score = 0.3;
	int best_edge = -1;
	for (int i = 0; i < n_edge_adj; i++) {
		int des = list_edge_adj[i].des;
		if (des == *start_contig)
			continue;
		if (!mark[des]) continue;
		struct pair_contigs_score *score = get_score(g, path, *start_contig, &list_edge_adj[i], edges_score);
		if (better_edge(score, max_score)) {
			max_score = score;
			best_edge = i;
		}
	}
	if (best_edge != -1) {
		int i_contig = list_edge_adj[best_edge].des;
		add_i_contig(path, i_contig);
		(*found)++;
		*start_contig = i_contig;
		mark[i_contig]--;
		int rc = get_rc_id(g, i_contig); 
		mark[rc]--;
		return;
	}
	*found = 0;
	//todo @huu: mark by cov
}

void unmark(struct asm_graph_t *g, struct scaffold_path *path, int *mark)
{
	for (int i = 0; i < path->n_contig; i++) {
		int c = path->list_i_contig[i], rc = get_rc_id(g, c);
		mark[c]++;
		mark[rc]++;
		VERBOSE_FLAG(0, "mark %d %d\n", mark[c], mark[rc]);
	}
}

void domark(struct asm_graph_t *g, struct scaffold_path *path, int *mark)
{
	for (int i = 0; i < path->n_contig; i++) {
		int c = path->list_i_contig[i], rc = get_rc_id(g, c);
		mark[c]--;
		mark[rc]--;
		VERBOSE_FLAG(0, "mark %d %d\n", mark[c], mark[rc]);
	}
}

void find_scaffolds(struct asm_graph_t *g,struct opt_proc_t *opt, struct edges_score_type *edges_score,
 		struct scaffold_type *scaffold)
{
	VERBOSE_FLAG(0, "start find scaffold\n");
	int *mark = calloc(g->n_e, sizeof(int));
	for (int i = 0; i < g->n_e; i++) {
		float edge_cov = __get_edge_cov(&g->edges[i], g->ksize)/global_genome_coverage;
		mark[i] = lround(edge_cov);
	}
	int found = 0;
	int start = 0;
	struct scaffold_path *best_path;
	do{
		start = 0;
		best_path = calloc(1, sizeof(struct scaffold_path));
		while (1) {
			int start_contig = get_highest_cov_contig(g, mark, start);
			start = start_contig + 1;
			if (start_contig == -1)
				break;
			struct scaffold_path *path = calloc(1, sizeof(struct scaffold_path));
			mark[start_contig]--;
			mark[get_rc_id(g, start_contig)]--;
			add_i_contig(path, start_contig);
			VERBOSE_FLAG(1, "start find scaffolds from %d\n", start_contig);
			int len_component = 0;
			do {
				find_best_edge(g, edges_score, &start_contig, path, &found, mark);
				VERBOSE_FLAG(1, "find path %d ", start_contig);
			} while (found);
			unmark(g, path, mark);
			if (path->n_contig > best_path->n_contig) {
				struct scaffold_path *tmp = best_path;
				best_path = path;
				destroy_path(tmp);
			}
		}
		add_path(scaffold, best_path);
		domark(g, best_path, mark);
	} while (best_path->n_contig > 0);
	for (int i = 0; i < g->n_e; i++) if (!is_very_short_contig(&g->edges[i]) && mark[i]) {
		struct scaffold_path *path = calloc(1, sizeof(struct scaffold_path));
		add_path(scaffold, path);
	}
	print_scaffold_contig(scaffold);
}

void insert_short_contig()
{
	//todo @huu
}

void scaffolding(FILE *out_file, struct asm_graph_t *g,
		struct opt_proc_t *opt) 
{
	init_global_params(g);
	VERBOSE_FLAG(0, "init global params done\n");
	check_global_params(g);
//	if (!opt->metagenomics) todo @huu uncomment
//		remove_lov_cov(g);

	float cvr = global_genome_coverage;
	for (int i_e = 0; i_e < g->n_e; i_e++) {
		float edge_cov = __get_edge_cov(&g->edges[i_e], g->ksize)/cvr;
		VERBOSE_FLAG(0, "edge %d len:%d cov: %f\n", i_e , get_edge_len(&g->edges[i_e]), edge_cov);
	}

	struct edges_score_type *edges_score = calloc(1, sizeof(struct edges_score_type));
	struct scaffold_type *scaffold = new_scaffold_type();
	pre_calc_score(g, opt, edges_score);
	print_edge_score(edges_score);

//	normalize_min_index(g, edges_score);
//	qsort(edges_score->list_edge, edges_score->n_edge, sizeof(struct scaffold_edge), 
//			ascending_scaffold_edge_index);
//	unique_edge(edges_score->list_edge, &edges_score->n_edge);
//	normalize_one_dir(g, edges_score);

	sort_edges_score(edges_score);
	find_scaffolds(g, opt, edges_score, scaffold);
	print_scaffold(g, out_file, scaffold);
	insert_short_contig();
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

