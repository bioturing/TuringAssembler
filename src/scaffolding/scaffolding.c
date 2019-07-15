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

static pthread_mutex_t lock_id = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t lock_put_table = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t lock_append_edges = PTHREAD_MUTEX_INITIALIZER;


struct params {
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
	float avg_bin_hash;
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
		int check = 0;
		if (e0 == e1) {
			// todo @huu what about metagenomics
			check = 1;
		} else if (e0 == g->edges[e1].rc_id) {
		} else {
			check = 1;
		}
		if (check) {
		} else {
			struct pair_contigs_score *sc = calloc(1, sizeof(struct pair_contigs_score));
			sc->bc_score =-1;
			pa_check_edge->list_candidate_edges[i_candidate].score = *sc;
		}
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
	float avg_bin_hash;
	khash_t(big_table) *big_table;
};

void destroy_params_build_candidate(struct params_build_candidate_edges *para)
{
	free(para->list_candidate_edges);
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
		int value = count[i_contig] ;
		*list_local_edges = realloc(*list_local_edges, (*n_local_edges+1) *
				sizeof(struct candidate_edge));
		struct candidate_edge *new_candidate_edge = calloc(1, sizeof(struct candidate_edge));
		new_candidate_edge->src = i_edge;
		new_candidate_edge->des = i_contig;
		if (value != 0) {
			int cnt0 = g->edges[get_rc_id(g, i_edge)].barcodes.n_item ;
			int cnt1 = g->edges[i_contig].barcodes.n_item;
			new_candidate_edge->score.bc_score = get_bc_score(value, cnt0, cnt1, params->avg_bin_hash); 
			(*list_local_edges)[*n_local_edges] = *new_candidate_edge;
			++*n_local_edges;
		}
		free(new_candidate_edge);
	}

	qsort(*list_local_edges, *n_local_edges, sizeof(struct candidate_edge), decending_candidate_edge);
	*n_local_edges = MIN(global_n_candidate, *n_local_edges);
	for (int i = 0; i < *n_local_edges; i++){
		struct candidate_edge e = (*list_local_edges)[i];
		if ((*list_local_edges)[i].score.bc_score == 0)
		{
			for (int j = 0; j < *n_local_edges; j++) {
				struct candidate_edge e_j = (*list_local_edges)[j];
				VERBOSE_FLAG(0, "ffff %d\n", e_j.score.bc_score);
			}
			*n_local_edges = i;
			break;
		}
	}
	free(count);
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
				struct list_position *pos = NULL;
				pos = kh_value(big_table, k);
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
	free(params_build_table);

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
			params_candidate->n_candidate_edges += n_local_edges; 
		COPY_ARR(list_local_edges, params_candidate->list_candidate_edges+n,
				n_local_edges);
		pthread_mutex_unlock(&lock_append_edges);
		VERBOSE_FLAG(0, "nlocaledge from %d %d \n", i_contig, n_local_edges);
		free(list_local_edges);
	} while(1);
	pthread_exit(NULL);
}

struct params_build_candidate_edges* new_params_build_candidate_edges(
	struct asm_graph_t *g, struct opt_proc_t *opt, khash_t(big_table) *big_table)
{
	struct params_build_candidate_edges *params_candidate = 
		calloc(1, sizeof(struct params_build_candidate_edges));
	params_candidate->g = g;
	params_candidate->i = 0;
	params_candidate->big_table = big_table;
	params_candidate->n_candidate_edges = 0;
	params_candidate->list_candidate_edges = NULL;
	params_candidate->avg_bin_hash = get_avg_barcode(g);
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

void remove_lov_high_cov(struct asm_graph_t *g)
{
	float cvr = global_genome_coverage;
	for (int i_e = 0; i_e < g->n_e; i_e++) {
		float edge_cov = __get_edge_cov(&g->edges[i_e], g->ksize)/cvr;
		if (edge_cov >= 3 || edge_cov <= 0.25) {
			g->edges[i_e].seq_len = 0;
		}
//		if (edge_cov <3 ){
//			int rc_id = g->edges[i_e].rc_id;
//			g->edges[rc_id].rc_id = count;
//			g->edges[count] =  g->edges[i_e];
//			count++;
//		}
	}
}

void pre_calc_score(struct asm_graph_t *g,struct opt_proc_t* opt, struct edges_score_type *edges_score)
{ 
	VERBOSE_FLAG(1, "start build big table");
	khash_t(big_table) *big_table = build_big_table(g, opt);
	struct params_build_candidate_edges *params_candidate = 
		new_params_build_candidate_edges(g, opt, big_table);
	run_parallel_build_candidate_edges(params_candidate, opt->n_threads);
	// sorted candidate edge
	VERBOSE_FLAG(1, "n candidate edges %d", params_candidate->n_candidate_edges);
	VERBOSE_FLAG(1, "done build candidate edge");

	struct params_check_edge *para = new_params_check_edge(g, opt, 
			params_candidate);
	VERBOSE_FLAG(1, "avg bin hash %f\n", para->avg_bin_hash);
	parallel_build_edge_score(para, opt->n_threads);

	VERBOSE_FLAG(1, "find real edge");
	int i_candidate_edge = 0;
	qsort(para->list_candidate_edges, para->n_candidate_edges, sizeof(struct candidate_edge),
			ascending_index_edge);
	for (int i_contig = 0; i_contig < g->n_e; i_contig++) if (!is_very_short_contig(&g->edges[i_contig])) {
		VERBOSE_FLAG(0, "i contig %d\n", i_contig);
		struct scaffold_edge *list_edges = NULL;
		int size_list_edge = 0; 
		int n_edges = 0;
		while (i_candidate_edge < para->n_candidate_edges && 
			para->list_candidate_edges[i_candidate_edge].src <= i_contig) {
			int src = para->list_candidate_edges[i_candidate_edge].src;
			assert(src == i_contig);
			int des = para->list_candidate_edges[i_candidate_edge].des;
			struct pair_contigs_score score = para->list_candidate_edges[i_candidate_edge].score;
			score.m_score = get_share_mate(g, src, des);
			VERBOSE_FLAG(0, "i candidate %d src %d des %d score %f\n", i_candidate_edge,
					para->list_candidate_edges[i_candidate_edge].src, des, score.bc_score);
			struct scaffold_edge *new_edge = new_scaffold_edge(i_contig, des, &score);
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
			if (score->bc_score == 0) {
				break;
			}
		}
		free(list_edges);
	}
	free(para);
	destroy_params_build_candidate(params_candidate);
	kh_destroy(big_table, big_table);
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
	struct pair_contigs_score *sc = calloc(1, sizeof(struct pair_contigs_score));
	int l = 0, r = edges_score->n_edge - 1;
	while (l != r) {
		int mid = (l + r)/2;
		if (edges_score->list_edge[mid].src < src) 
			l = mid+1;
		else
			r = mid;
	}
	while (l < edges_score->n_edge && edges_score->list_edge[l].src == src) {
		struct scaffold_edge *edge = &edges_score->list_edge[l];
		if (edge->des == des) {
			*sc = edge->score;
			return sc;
		}
		l++;
	}
	sc->bc_score = 0;
	return sc;
}

struct pair_contigs_score *divn(struct pair_contigs_score *score, int n)
{
	struct pair_contigs_score *res = calloc(1, sizeof(struct pair_contigs_score));
	if (n == 0) {
		return res;
	}
	res->bc_score = score->bc_score/n;
	res->m_score = score->m_score/n;
	res->m2_score = score->m2_score/n;
	return res;
}

struct pair_contigs_score *get_score(struct asm_graph_t *g, struct scaffold_path *path, int start_contig, 
	struct scaffold_edge *edge, struct edges_score_type *edges_score, int is_left)
{
	struct pair_contigs_score *score = calloc(1, sizeof(struct pair_contigs_score));
	struct pair_contigs_score *second_score = calloc(1, sizeof(struct pair_contigs_score));
	int des = edge->des;
	struct pair_contigs_score *i_score; 
	int last = get_last_n(path, is_left, 0);
	if (is_left) last = get_rc_id(g, last);
	i_score = get_score_edge(edges_score,last , des);
	*score = *i_score;
	i_score = get_score_edge(edges_score, last, g->edges[des].rc_id);
	score->bc_score += i_score->bc_score/2;
	int i = 0;
	int distance = get_edge_len(&g->edges[last]);
	VERBOSE_FLAG(0, "scascore %d %d %f\n", last, des, score->bc_score);
	while (1) {
		if (distance > global_distance)
			break;
		int src = get_last_n(path, is_left, i);
		if (src == -1)
			break;
		if (is_left) 
			src = get_rc_id(g, src);
		i_score = get_score_edge(edges_score, src, des);
		second_score->bc_score += i_score->bc_score;
		VERBOSE_FLAG(0, "more %d %f ", src, second_score->bc_score);
		distance += get_edge_len(&g->edges[src]);
	}
	if (i != 0)
		score->bc_score += second_score->bc_score/(i*3);
	VERBOSE_FLAG(0, "\ndonescascore %f %d %d\n", score->bc_score, score->m_score, score->m2_score);
	free(second_score);
	//todo @huu MAX(i/2, 1) because far contig have less score
	return score;
}

int better_edge(struct pair_contigs_score *sc0, struct pair_contigs_score *sc1)
{
	//todo wtd
	if (sc0->bc_score > 1.1 * sc1->bc_score) return 1;
	if (sc0->bc_score * 1.1 < sc1->bc_score) return 0;
	if (sc0->m_score + sc0->m2_score != sc1->m_score + sc1->m2_score)
		return (sc0->m_score + sc0->m2_score > sc1->m_score +sc1->m2_score);
	return (sc0->bc_score > sc1->bc_score);
}

struct best_next_contig {
	int i_contig;
	struct pair_contigs_score score;
};

struct best_next_contig *find_best_edge(struct asm_graph_t *g, struct edges_score_type *edges_score, int start_contig,
		struct scaffold_path *path, int *mark, int is_left, struct pair_contigs_score *thres_score)
{
	int n_edge_adj = 0;
	struct scaffold_edge *list_edge_adj = NULL;
	find_edge_from(edges_score, start_contig, &n_edge_adj, &list_edge_adj);
	struct pair_contigs_score *max_score = calloc(1, sizeof(struct pair_contigs_score));
	//todo @huu dynamic this thres
//	max_score->bc_score = 0.05;
	VERBOSE_FLAG(0, "find best edge from %d\n", start_contig);
	int best_edge = -1;
	for (int i = 0; i < n_edge_adj; i++) {
		int des = list_edge_adj[i].des;
		VERBOSE_FLAG(0, "des %d\n", des);
		if (des == start_contig)
			continue;
		if (!mark[des]) continue;
		struct pair_contigs_score *score = get_score(g, path, start_contig, &list_edge_adj[i], edges_score, is_left);
		struct pair_contigs_score *tmp = score;
		if (better_edge(score, max_score)) {
			tmp = max_score;
			max_score = score;
			best_edge = des;
		}
		free(tmp);
	}
	struct best_next_contig *res = calloc(1, sizeof(struct best_next_contig));
	res->i_contig = best_edge;
	res->score = *max_score;
	if (!better_edge(&res->score, thres_score))
		res->i_contig = -1;
	VERBOSE_FLAG(0,"res icontig %d thres score %f\n\n", res->i_contig, thres_score->bc_score);
	return res;
}

//void unmark(struct asm_graph_t *g, struct scaffold_path *path, int *mark)
//{
//	for (int i = 0; i < path->n_contig; i++) {
//		int c = path->list_i_contig[i], rc = get_rc_id(g, c);
//		mark[c]++;
//		mark[rc]++;
//		VERBOSE_FLAG(0, "mark %d %d\n", mark[c], mark[rc]);
//	}
//}
//
//void domark(struct asm_graph_t *g, struct scaffold_path *path, int *mark)
//{
//	for (int i = 0; i < path->n_contig; i++) {
//		int c = path->list_i_contig[i], rc = get_rc_id(g, c);
//		mark[c]--;
//		mark[rc]--;
//		VERBOSE_FLAG(0, "mark %d %d\n", mark[c], mark[rc]);
//	}
//}

void mark_contig(struct asm_graph_t *g, int *mark, int i_contig)
{
	mark[i_contig]--;
	mark[get_rc_id(g, i_contig)]--;
}

struct pair_contigs_score *get_score_tripple(struct edges_score_type *edges_score, int l, int m, int r)
{
	struct pair_contigs_score *score0 = get_score_edge(edges_score, l, m);
	struct pair_contigs_score *score1 = get_score_edge(edges_score, m, r);
	struct pair_contigs_score *res = calloc(1, sizeof(struct pair_contigs_score));
	res->bc_score = score0->bc_score + score1->bc_score;
	res->m_score = score0->m_score + score1->m_score;
	res->m2_score = score0->m2_score + score1->m2_score;
	free(score0);
	free(score1);
	return res;
}

void refine_path(struct asm_graph_t *g, struct edges_score_type *edges_score, struct scaffold_path *path) 
{
	//todo @huu refine base on nearby contig in distance not by 2 next to contig
	int n = path->n_left_half + path->n_right_half;
	for (int j = 1; j < n-1; j++) {
		int left = get_last_n(path, 1, j-1);
		int mid = get_last_n(path, 1, j);
		int right = get_last_n(path, 1, j+1);
		struct pair_contigs_score *normal_score = get_score_tripple(edges_score, left, mid, right);
		struct pair_contigs_score *reverse_score = get_score_tripple(edges_score, left, get_rc_id(g, mid), right);
		if (better_edge(reverse_score, normal_score)) {
			reverse_n_th(g, path, 1, j);
		}
		free(normal_score);
		free(reverse_score);
	}
}

void refine_scaffold(struct asm_graph_t *g, struct edges_score_type *edges_score, struct scaffold_type *scaffold) 
{
	for (int i = 0; i < scaffold->n_path; i++) {
		struct scaffold_path *path = &scaffold->path[i];
		refine_path(g, edges_score, path);
	}
}

struct scaffold_path *find_path(struct opt_proc_t *opt, struct asm_graph_t *g, 
		struct edges_score_type *edges_score, int *mark, int start_contig,
		struct pair_contigs_score *thres_score, int *count)
{
	struct scaffold_path *path = calloc(1, sizeof(struct scaffold_path));
	mark_contig(g, mark, start_contig);
//		add_i_contig(path, start_contig);
	append_i_contig(path, start_contig);
	int i_r_contig = start_contig, i_l_contig = get_rc_id(g, start_contig);
	if (opt->metagenomics) {
		thres_score->bc_score = 0;
		thres_score->m_score = 0;
		thres_score->m2_score = 0;
		*count = 0;
	}
	while (1) {
		struct best_next_contig *next_l_contig, *next_r_contig, *next_contig;
		next_l_contig = find_best_edge(g, edges_score, i_l_contig, path, mark, 1, divn(thres_score, 5*(*count)));
		next_r_contig = find_best_edge(g, edges_score, i_r_contig, path, mark, 0, divn(thres_score, 5*(*count)));

		VERBOSE_FLAG(0, "next l contig %d next r contig %d\n", next_l_contig->i_contig, next_r_contig->i_contig);
		if (next_r_contig->i_contig == -1 && next_l_contig->i_contig == -1) {
			break;
		}
		if (next_r_contig->i_contig == -1 || better_edge(&next_l_contig->score, &next_r_contig->score)) {
			prepend_i_contig(path, get_rc_id(g, next_l_contig->i_contig));
			i_l_contig = next_l_contig->i_contig;
			next_contig = next_l_contig;
		} else {
			append_i_contig(path, next_r_contig->i_contig);
			i_r_contig = next_r_contig->i_contig;
			next_contig = next_r_contig;
		}
		mark_contig(g, mark, next_contig->i_contig);
		thres_score->bc_score = (thres_score->bc_score + next_contig->score.bc_score);
		thres_score->m_score = (thres_score->m_score + next_contig->score.m_score);
		thres_score->m2_score = (thres_score->m2_score + next_contig->score.m2_score);
		(*count)++;
	}
	return path;
}

void init_mark(struct asm_graph_t *g, struct opt_proc_t *opt, int *mark)
{
	if (!opt->metagenomics) {
		for (int i = 0; i < g->n_e; i++) {
			float edge_cov = __get_edge_cov(&g->edges[i], g->ksize)/global_genome_coverage;
			mark[i] = MIN(lround(edge_cov), 3);
		}
	} else {
		for (int i = 0; i < g->n_e; i++) {
			mark[i] = 1;
		}
	}
}

void find_scaffolds(struct asm_graph_t *g,struct opt_proc_t *opt, struct edges_score_type *edges_score,
 		struct scaffold_type *scaffold)
{
	VERBOSE_FLAG(0, "start find scaffold\n");
	int *mark = calloc(g->n_e, sizeof(int));
	init_mark(g, opt, mark);
	int count = 0;
	struct pair_contigs_score *thres_score = calloc(1, sizeof(struct pair_contigs_score));
	for (int i = 0; i < g->n_e; i++) if (mark[i] && is_long_contig(&g->edges[i])){
		int start_contig = i;
		VERBOSE_FLAG(1, "start find scaffolds from %d\n", start_contig);
		struct scaffold_path *path = find_path(opt, g, edges_score, mark, start_contig, 
							thres_score, &count);
		add_path(scaffold, path);
		free(path);
	}
	for (int i = 0; i < g->n_e; i++) if (is_short_contig(&g->edges[i]) && mark[i]) {
		struct scaffold_path *path = calloc(1, sizeof(struct scaffold_path));
		append_i_contig(path, i);
		add_path(scaffold, path);
	}
	free(mark);
	free(thres_score);
	print_scaffold_contig(scaffold);
}

void insert_short_contig()
{
	//todo @huu
}

int count_bc(struct asm_edge_t *e)
{
	struct barcode_hash_t *buck = &e->barcodes;
	int count = 0;
	for (int j = 0; j < buck->size; j++){
		if (buck->keys[j] != (uint64_t)(-1)) {
			count++;
		}
	} 
	return count;
}

void scaffolding(FILE *out_file, struct asm_graph_t *g,
		struct opt_proc_t *opt) 
{
	init_global_params(g);
	VERBOSE_FLAG(0, "init global params done\n");
	check_global_params(g);
	if (!opt->metagenomics) {
		remove_lov_high_cov(g);
	}

	float cvr = global_genome_coverage;
	for (int i_e = 0; i_e < g->n_e; i_e++) {
		struct asm_edge_t *edge = &g->edges[i_e];
		float edge_cov = __get_edge_cov(edge, g->ksize)/cvr;
		VERBOSE_FLAG(0, "edge %d len:%d cov: %f count_bc %d\n", 
				i_e , get_edge_len(&g->edges[i_e]), edge_cov, g->edges[i_e].barcodes.n_item);
	}

	struct edges_score_type *edges_score = calloc(1, sizeof(struct edges_score_type));
	struct scaffold_type *scaffold = new_scaffold_type();
	pre_calc_score(g, opt, edges_score);

	sort_edges_score(edges_score);
	print_edge_score(edges_score);
	find_scaffolds(g, opt, edges_score, scaffold);
	insert_short_contig();
	refine_scaffold(g, edges_score, scaffold);
	print_scaffold_contig(scaffold);
	print_scaffold(g, out_file, scaffold);
}

void scaffolding_test(struct asm_graph_t *g, struct opt_proc_t *opt)
{
	init_global_params(g);
	__VERBOSE("n_e: %ld\n", g->n_e);
	for (int i = 0; i < g->n_e; i++)
		__VERBOSE("edge %d length %d\n", i, g->edges[i].seq_len);

	check_global_params(g);
}

