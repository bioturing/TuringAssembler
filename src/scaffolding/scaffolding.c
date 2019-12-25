#include <stdlib.h>
#include "assembly_graph.h"
#include "pthread.h"
#include "verbose.h"
#include <string.h>
#include <assert.h>
#include <barcode_builder.h>
#include <barcode_resolve2.h>
#include <minimizers/minimizers.h>
#include <cluster_molecules.h>
#include <barcode_graph.h>
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
#include "scaffolding.h"
#include "smart_load.h"
#include "log.h"
#include "yeast_analyze_utils.h"

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
    struct scaffold_edge *list_candidate_edges;
    float avg_bin_hash;
    FILE *out_file;
    struct opt_proc_t *opt;
};

void destroy_params_check_edge(struct params_check_edge *para)
{
	free(para->list_contig);
}

int too_different(float a, float b)
{
	if ((a < 1.0 / 3 * b) || (a > 3 * b))
		return 1;
	return 0;
}

pthread_attr_t init_thread_attr()
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	return attr;
}

struct list_position {
    int n_pos;
    int *i_contig, *i_bin;
    pthread_mutex_t lock_entry;
};

KHASH_MAP_INIT_INT64(btable_sig, struct list_position*);

struct params_build_candidate_edges {
    struct asm_graph_t *g;
    int i;
    struct edges_score_type *list_candidate_edges;
    float avg_bin_hash;
    int metagenomics;
    khash_t(btable_sig) *big_table;
};

void destroy_params_build_candidate(struct params_build_candidate_edges *para)
{
//	free(para->list_candidate_edges); this array is keep
	free(para);
}

int count_pos(int *count, struct list_position *pos)
{
	assert(pos != NULL && count != NULL);
	int res = 0;
	for (int i = 0; i < pos->n_pos; i++) {
		res++;
		count[pos->i_contig[i]]++;
	}
	return res;
}

void find_local_nearby_contig(int i_edge, struct params_build_candidate_edges *params, int *n_local_edges,
                              struct scaffold_edge **list_local_edges)
{
	struct asm_graph_t *g = params->g;
	int i_rev_edge = g->edges[i_edge].rc_id;
	struct asm_edge_t *e = &params->g->edges[i_edge];
	struct asm_edge_t *rev_e = &params->g->edges[i_rev_edge];
	//todo calloc n long contigs to reduce RAMs when scale to big genome
	int *count = calloc(params->g->n_e, sizeof(int));

	khash_t(btable_sig) *big_table = params->big_table;

	struct barcode_hash_t *buck = &rev_e->barcodes_scaf;
	for (int j = 0; (uint32_t) j < buck->size; j++) {
		if (buck->keys[j] != (uint64_t) (-1)) {
			uint64_t barcode = buck->keys[j];
			khint_t k = kh_get(btable_sig, big_table, barcode);
			if (k == kh_end(big_table))
				continue;
			struct list_position *pos = kh_value(big_table, k);
			count_pos(count, pos);
		}
	}

	for (int i_contig = 0; i_contig < g->n_e; i_contig++) {
		if (is_very_short_contig(&g->edges[i_contig]))
			continue;
		*list_local_edges = realloc(*list_local_edges, (*n_local_edges + 1) *
		                                               sizeof(struct scaffold_edge));
		struct scaffold_edge *new_candidate_edge = calloc(1, sizeof(struct scaffold_edge));
		if ((i_contig == i_edge && g->edges[i_contig].seq_len < 50000) || i_contig == get_rc_id(g, i_edge)) {
			free(new_candidate_edge);
			continue;
		}
		new_candidate_edge->src = i_edge;
		new_candidate_edge->des = i_contig;
		float e1_cov = __get_edge_cov(&g->edges[i_edge], g->ksize);
		float e2_cov = __get_edge_cov(&g->edges[i_contig], g->ksize);
		int value = count[i_contig];
		if (too_different(e1_cov, e2_cov))
			value = 0;
		if (value != 0) {
			int cnt0 = g->edges[get_rc_id(g, i_edge)].barcodes_scaf.n_item;
			int cnt1 = g->edges[i_contig].barcodes_scaf.n_item;
			new_candidate_edge->score.bc_score = get_bc_score(value, cnt0, cnt1, params->avg_bin_hash, i_edge, i_contig);
			(*list_local_edges)[*n_local_edges] = *new_candidate_edge;
			++*n_local_edges;
		}
		free(new_candidate_edge);
	}

	qsort(*list_local_edges, *n_local_edges, sizeof(struct scaffold_edge), decending_scaffold_edge);
	*n_local_edges = MIN(global_n_candidate, *n_local_edges);
	for (int i = 0; i < *n_local_edges; i++) {
		struct scaffold_edge e = (*list_local_edges)[i];
		if ((*list_local_edges)[i].score.bc_score == 0
		     || (i > 0 && (*list_local_edges)[i].score.bc_score < 0.5 * (*list_local_edges)[i-1].score.bc_score)) {
			for (int j = 0; j < *n_local_edges; j++) {
				struct scaffold_edge e_j = (*list_local_edges)[j];
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
    khash_t(btable_sig) *big_table;
};

void *process_build_big_table(void *data)
{
	struct params_build_big_table *params = (struct params_build_big_table *) data;
	struct asm_graph_t *g = params->g;
	khash_t(btable_sig) *big_table = params->big_table;

	do {
		pthread_mutex_lock(&lock_id);
		if (params->i == g->n_e) {
			pthread_mutex_unlock(&lock_id);
			break;
		}
		int i_contig = params->i;
		params->i++;
		pthread_mutex_unlock(&lock_id);
		int new_i_contig = i_contig;

		struct asm_edge_t *e = &g->edges[i_contig];
		if (!is_long_contig(e))
			continue;
		struct barcode_hash_t *buck = &e->barcodes_scaf;
		for (int l = 0; (uint32_t) l < buck->size; l++) {
			if (buck->keys[l] != (uint64_t) (-1)) {
				uint64_t barcode = buck->keys[l];
				pthread_mutex_lock(&lock_put_table);
				khint_t k = kh_get(btable_sig, big_table, barcode);
				if (k == kh_end(big_table)) {
					int tmp = 1;
					k = kh_put(btable_sig, big_table, barcode, &tmp);
					kh_value(big_table, k) = NULL;
					assert(tmp == 1);
				}
				struct list_position *pos = kh_value(big_table, k);
				if (pos == NULL) {
					pos = calloc(1, sizeof(struct list_position));
					kh_value(big_table, k) = pos;
					pos->lock_entry = (pthread_mutex_t) PTHREAD_MUTEX_INITIALIZER;
				}
				pthread_mutex_unlock(&lock_put_table);
				pthread_mutex_lock(&pos->lock_entry);
				int v = pos->n_pos;
				pos->n_pos++;
				if ((v & (v - 1)) == 0) {
					pos->i_contig = realloc(pos->i_contig, (v * 2 + 1) *
					                                       sizeof(int));
					pos->i_bin = realloc(pos->i_bin, (v * 2 + 1) *
					                                 sizeof(int));
				}
				pos->i_contig[v] = new_i_contig;
				pthread_mutex_unlock(&pos->lock_entry);
			}
		}
	} while (1);
	return NULL;
}

khash_t(btable_sig) *build_big_table(struct asm_graph_t *g, struct opt_proc_t *opt)
{
	log_info("----- Start build big table ------");
	pthread_t *thr = (pthread_t *) calloc(opt->n_threads, sizeof(pthread_t));
	pthread_attr_t attr = init_thread_attr();

	struct params_build_big_table *params_build_table = calloc(1, sizeof(struct
		params_build_big_table));
	params_build_table->i = 0;
	params_build_table->g = g;
	params_build_table->big_table = kh_init(btable_sig);
	// todo @huu auto resize
	kh_resize(btable_sig, params_build_table->big_table, 100000000);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(&thr[i], &attr, process_build_big_table, params_build_table);
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(thr[i], NULL);

	khash_t(btable_sig) *big_table = params_build_table->big_table;
	free(params_build_table);
	free(thr);

	log_info("----- Build done -------");
	return big_table;
}

void *process_build_candidate_edges(void *data)
{
	struct params_build_candidate_edges *params_candidate;
	params_candidate = (struct params_build_candidate_edges *) data;
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
		struct scaffold_edge *list_local_edges = NULL;
		if (!is_long_contig(&g->edges[i_contig]))
			continue;
		if (!params_candidate->metagenomics) {
			float edge_cov = __get_edge_cov(&g->edges[i_contig], g->ksize) / global_genome_coverage;
			if (edge_cov < MIN_EDGE_COV_SCAFFOLD) {
				continue;
			}
		}
		find_local_nearby_contig(i_contig, params_candidate,
		                         &n_local_edges,
		                         &list_local_edges);
		pthread_mutex_lock(&lock_append_edges);

		for (int i = 0; i < n_local_edges; i++) {
			append_edge_score(params_candidate->list_candidate_edges, &list_local_edges[i]);
		}
		pthread_mutex_unlock(&lock_append_edges);
		free(list_local_edges);
	} while (1);
	pthread_exit(NULL);
}

struct params_build_candidate_edges *new_params_build_candidate_edges(
	struct asm_graph_t *g, struct opt_proc_t *opt, khash_t(btable_sig) *big_table)
{
	struct params_build_candidate_edges *params_candidate =
		calloc(1, sizeof(struct params_build_candidate_edges));
	params_candidate->g = g;
	params_candidate->i = 0;
	params_candidate->big_table = big_table;
	params_candidate->list_candidate_edges = calloc(1, sizeof(struct edges_score_type));
	params_candidate->avg_bin_hash = get_avg_barcode(g);
	params_candidate->metagenomics = opt->metagenomics;
	return params_candidate;
}

void run_parallel_build_candidate_edges(struct params_build_candidate_edges *params_candidate, int n_threads)
{
	assert(params_candidate != NULL);
	pthread_t *thr = (pthread_t *) calloc(n_threads, sizeof(pthread_t));
	pthread_attr_t attr = init_thread_attr();
	for (int i = 0; i < n_threads; ++i)
		pthread_create(&thr[i], &attr, process_build_candidate_edges, params_candidate);
	for (int i = 0; i < n_threads; ++i)
		pthread_join(thr[i], NULL);
	free(thr);
	pthread_attr_destroy(&attr);
}

void remove_lov_high_cov(struct asm_graph_t *g)
{
	float cvr = global_genome_coverage;
	int64_t total_len = 0;
	for (int i_e = 0; i_e < g->n_e; i_e++) {
		float edge_cov = __get_edge_cov(&g->edges[i_e], g->ksize) / cvr;
		if (edge_cov < MIN_EDGE_COV_SCAFFOLD) {
			total_len += g->edges[i_e].seq_len;
			g->edges[i_e].seq_len = 0;
			log_debug("Remove edge: %d, coverage: %.2f, threshold %.2f", i_e, edge_cov, MIN_EDGE_COV_SCAFFOLD);
		}
	}
	log_info("remove %ld bp because have lower than 0.25 cov", total_len);
}

void deep_kh_destroy(khash_t(btable_sig) *big_table)
{
	for (khiter_t it = kh_begin(big_table); it != kh_end(big_table); it++) {
		if (!(kh_exist(big_table, it))) {
			continue;
		}
		struct list_position *pos = kh_value(big_table, it);
		free(pos->i_bin);
		free(pos->i_contig);
		free(pos);
	}
	kh_destroy(btable_sig, big_table);
}

void calc_score_pairwise(struct asm_graph_t *g, struct opt_proc_t *opt, struct edges_score_type **edges_score)
{
	log_info("----- Build candidate edge -----");
	khash_t(btable_sig) *big_table = build_big_table(g, opt);
	struct params_build_candidate_edges *params_candidate =
		new_params_build_candidate_edges(g, opt, big_table);
	run_parallel_build_candidate_edges(params_candidate, opt->n_threads);
	// sorted candidate edge
	log_info("----- Done build candidate edge -----");

	*edges_score = params_candidate->list_candidate_edges;

	destroy_params_build_candidate(params_candidate);
	deep_kh_destroy(big_table);
}

struct pair_contigs_score *get_score_edge(struct edges_score_type *edges_score, int src, int des)
{
	struct pair_contigs_score *sc = calloc(1, sizeof(struct pair_contigs_score));
	int l = 0, r = edges_score->n_edge - 1;
	while (l != r) {
		int mid = (l + r) / 2;
		if (edges_score->list_edge[mid].src < src)
			l = mid + 1;
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
	return res;
}

struct pair_contigs_score *get_score(struct asm_graph_t *g, struct scaffold_path *path,
                                     struct scaffold_edge *edge, struct edges_score_type *edges_score, int is_left)
{
	struct pair_contigs_score *score = calloc(1, sizeof(struct pair_contigs_score));
	struct pair_contigs_score *second_score = calloc(1, sizeof(struct pair_contigs_score));
	int des = edge->des;
	int last = get_last_n(path, is_left, 0);
	if (is_left) last = get_rc_id(g, last);
	score = get_score_edge(edges_score, last, des);
	struct pair_contigs_score *tmp_score = get_score_edge(edges_score, last, g->edges[des].rc_id);
	score->bc_score += tmp_score->bc_score / 2;
	int i = 0;
	int distance = get_edge_len(&g->edges[last]);
	log_trace("scascore %d %d %f", last, des, score->bc_score);
	while (1) {
		i++;
		int src = get_last_n(path, is_left, i);
		if (src == -1)
			break;
		if (is_left)
			src = get_rc_id(g, src);
		struct pair_contigs_score *tmp_score = get_score_edge(edges_score, src, des);
		second_score->bc_score += tmp_score->bc_score;
		log_trace("more %d %f ", src, second_score->bc_score);
//		if (tmp_score->bc_score ==0) {
//			log_trace("TERMINATE");
//			score->bc_score = 0;
//			return score;
//		}

		distance += get_edge_len(&g->edges[src]);
		free(tmp_score);
		if (distance > global_distance)
			break;
	}
	if (i != 0)
		score->bc_score += second_score->bc_score/(i*3);
	log_trace("donescascore %f", score->bc_score);
	free(second_score);
	free(tmp_score);
	//todo @huu MAX(i/2, 1) because far contig have less score
	return score;
}

int better_edge(struct pair_contigs_score *sc0, struct pair_contigs_score *sc1)
{
	//todo wtd
//	if (sc0->bc_score > 1.1 * sc1->bc_score) return 1;
//	if (sc0->bc_score * 1.1 < sc1->bc_score) return 0;
//	if (sc0->m_score + sc0->m2_score != sc1->m_score + sc1->m2_score)
//		return (sc0->m_score + sc0->m2_score > sc1->m_score +sc1->m2_score);
	return (sc0->bc_score > sc1->bc_score);
}

struct best_next_contig {
    int i_contig;
    struct pair_contigs_score score;
};

struct best_next_contig *find_best_edge(struct asm_graph_t *g, struct edges_score_type *edges_score, int start_contig,
                                        struct scaffold_path *path, int *mark, int is_left,
                                        struct pair_contigs_score *thres_score)
{
	int n_edge_adj = 0;
	struct scaffold_edge *list_edge_adj = NULL;
	find_edge_from(edges_score, start_contig, &n_edge_adj, &list_edge_adj);
	struct pair_contigs_score *max_score = calloc(1, sizeof(struct pair_contigs_score));
	//todo @huu dynamic this thres
//	max_score->bc_score = 0.05;
	log_trace("Find best edge from %d", start_contig);
	int best_edge = -1;
	for (int i = 0; i < n_edge_adj; i++) {
		int des = list_edge_adj[i].des;
		log_trace("des %d", des);
		if (des == start_contig)
			continue;
		if (!mark[des]) continue;
		struct pair_contigs_score *score = get_score(g, path, &list_edge_adj[i], edges_score, is_left);
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
	if (!better_edge(&res->score, thres_score)) {
		res->i_contig = -1;
		log_trace("below thres score %f", thres_score->bc_score);
	}
	log_trace("res icontig %d thres score %f", res->i_contig, thres_score->bc_score);
	free(max_score);
	return res;
}

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
	free(score0);
	free(score1);
	return res;
}

void refine_path(struct asm_graph_t *g, struct edges_score_type *edges_score, struct scaffold_path *path)
{
	//todo @huu refine base on nearby contig in distance not by 2 next to contig
	int n = path->n_left_half + path->n_right_half;
	for (int j = 1; j < n - 1; j++) {
		int left = get_last_n(path, 1, j - 1);
		int mid = get_last_n(path, 1, j);
		int right = get_last_n(path, 1, j + 1);
		struct pair_contigs_score *normal_score = get_score_tripple(edges_score, left, mid, right);
		struct pair_contigs_score *reverse_score = get_score_tripple(edges_score, left, get_rc_id(g, mid),
		                                                             right);
		log_trace("get score tripple %d %d %d normal %f reverse %f", left, mid, right, normal_score->bc_score,
		          reverse_score->bc_score);
		if (better_edge(reverse_score, normal_score)) {
			reverse_n_th(g, path, 1, j);
			j++;
		}
		free(normal_score);
		free(reverse_score);
	}
}

void refine_scaffold(struct asm_graph_t *g, struct edges_score_type *edges_score, struct scaffold_type *scaffold)
{
	log_info("----- Start refine scaffold------");
	for (int i = 0; i < scaffold->n_path; i++) {
		struct scaffold_path *path = &scaffold->path[i];
		refine_path(g, edges_score, path);
	}
	log_info("----- Refine scaffold done------");
}

struct scaffold_path *find_path(struct opt_proc_t *opt, struct asm_graph_t *g,
                                struct edges_score_type *edges_score, int *mark, int start_contig,
                                struct pair_contigs_score *thres_score, int *count)
{
	struct scaffold_path *path = calloc(1, sizeof(struct scaffold_path));
	mark_contig(g, mark, start_contig);
	append_i_contig(path, start_contig);
	int i_r_contig = start_contig, i_l_contig = get_rc_id(g, start_contig);
	if (opt->metagenomics) {
		thres_score->bc_score = 0;
		*count = 0;
	}
	while (1) {
		struct best_next_contig *next_l_contig, *next_r_contig, *next_contig;
		struct pair_contigs_score *thres_score_d5 = divn(thres_score, 5 * (*count));
		next_l_contig = find_best_edge(g, edges_score, i_l_contig, path, mark, 1, thres_score_d5);
		next_r_contig = find_best_edge(g, edges_score, i_r_contig, path, mark, 0, thres_score_d5);

		log_trace("next l contig %d next r contig %d", next_l_contig->i_contig, next_r_contig->i_contig);
		if (next_r_contig->i_contig == -1 && next_l_contig->i_contig == -1) {
			free(next_l_contig);
			free(next_r_contig);
			free(thres_score_d5);
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
		assert(g->edges[next_contig->i_contig].seq_len > 0);
		thres_score->bc_score = (thres_score->bc_score + next_contig->score.bc_score);
		free(next_l_contig);
		free(next_r_contig);
		free(thres_score_d5);
		(*count)++;
	}
	return path;
}

void init_mark(struct asm_graph_t *g, struct opt_proc_t *opt, int *mark)
{
	if (!opt->metagenomics) {
		for (int i = 0; i < g->n_e; i++) {
			float edge_cov =
				MAX(__get_edge_cov(&g->edges[i], g->ksize) / global_genome_coverage, 1);
			mark[i] = MIN(lround(edge_cov), 3);
		}
	} else {
		for (int i = 0; i < g->n_e; i++) {
			mark[i] = 1;
		}
	}
	for(int i = 0 ; i < g->n_e; i++) {
		assert(mark[i] == mark[g->edges[i].rc_id]);
	}
}

void find_scaffolds(struct asm_graph_t *g, struct opt_proc_t *opt, struct edges_score_type *edges_score,
                    struct scaffold_type *scaffold)
{
	log_info("----- Start find scaffold------");
	int *mark = calloc(g->n_e, sizeof(int));
	init_mark(g, opt, mark);
	int count = 0;
	struct pair_contigs_score *thres_score = calloc(1, sizeof(struct pair_contigs_score));
	//TODO: Start from the longest edge
	for (int i = 0; i < g->n_e; i++)
		if (mark[i] && is_long_contig(&g->edges[i])) {
			int start_contig = i;
			log_trace("Start find scaffolds from %d", start_contig);
			struct scaffold_path *path = find_path(opt, g, edges_score, mark, start_contig,
			                                       thres_score, &count);
			add_path(scaffold, path);
			free(path);
		}
	int64_t total_very_short = 0;
	for (int i = 0; i < g->n_e; i++) {
		if (is_short_contig(&g->edges[i]) && mark[i]) {
			struct scaffold_path *path = calloc(1, sizeof(struct scaffold_path));
			mark_contig(g, mark, i);
			append_i_contig(path, i);
			assert(g->edges[i].seq_len > 0);
			add_path(scaffold, path);
			free(path);
		} else if (is_very_short_contig(&g->edges[i])) {
			int len = g->edges[i].seq_len;
			total_very_short += len;
		}
	}
	log_info("contig shorter than %d have total length: %d", global_thres_short_len, total_very_short);
	free(mark);
	free(thres_score);
	log_info("----- Find scaffold done ------");
}

void print_contig_info(struct asm_graph_t *g)
{
	float cvr = global_genome_coverage;
	for (int i_e = 0; i_e < g->n_e; i_e++) {
		struct asm_edge_t *edge = &g->edges[i_e];
		float edge_cov = __get_edge_cov(edge, g->ksize) / cvr;
		log_debug("edge %d len:%d cov: %f count_bc %d",
		          i_e, get_edge_len(&g->edges[i_e]), edge_cov, g->edges[i_e].barcodes_scaf.n_item);
	}
}

int check_should_local_assembly(struct scaffold_type *scaffold)
{
	int lc = 0;
	for (int i = 0; i < scaffold->n_path; i++) {
		if (scaffold->path[i].n_left_half + scaffold->path[i].n_right_half > 1) {
			lc = 1;
			break;
		}
	}
	if (lc == 0) {
		log_warn("Too many low quality contigs for scaffolding. Stop program.");
		return 0;
	}
	return 1;
}

void copyfile(char *in_name, char *out_name)
{
	FILE *in_file = fopen(in_name, "rb");
	if (in_file == NULL) {
		log_error("file open for reading");
	}
	FILE *out_file = fopen(out_name, "wb");
	if (out_file == NULL) {
		log_error("file open for writing");
	}
	size_t n, m;
	unsigned char buff[8192];
	do {
		n = fread(buff, 1, sizeof buff, in_file);
		if (n) m = fwrite(buff, 1, n, out_file);
		else m = 0;
	} while ((n > 0) && (n == m));
	if (m) {
		log_error("copy");
	}
	if (fclose(out_file)) perror("close output file");
	if (fclose(in_file)) perror("close input file");
}

void scaffolding(FILE *out_file, struct asm_graph_t *g,
                 struct opt_proc_t *opt)
{
	init_global_params(g);

//	if (!opt->metagenomics) {
//		remove_lov_high_cov(g);
//	}

	print_contig_info(g);

	struct edges_score_type *edges_score = NULL;
	calc_score_pairwise(g, opt, &edges_score);

	struct scaffold_type *scaffold = new_scaffold_type();
	sort_edges_score(edges_score);
	print_edge_score(edges_score);
	find_scaffolds(g, opt, edges_score, scaffold);

	refine_scaffold(g, edges_score, scaffold);

	print_scaffold_contig(opt, scaffold); /* print scaffold path into local_assembly_scaffold_path.txt */
	print_scaffold(g, out_file, scaffold);
	int check_should_do_local = check_should_local_assembly(scaffold);
	destroy_scaffold_type(scaffold);
	destroy_edges_score_type(edges_score);
	if (!check_should_do_local) {
		char *in_name = str_concate(opt->out_dir, "/scaffolds.fasta");
		char *out_name = str_concate(opt->out_dir, "/scaffold.full.fasta");
		copyfile(in_name, out_name);
		exit(0);
	}
}

struct asm_graph_t *create_and_load_graph(struct opt_proc_t *opt)
{
	struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_file);
	test_asm_graph(g0);
	return g0;
}

void load_list_barcode(int *n_barcodes, char ***barcodes, int **freq)
{
	FILE *f = fopen("barcode_frequencies.txt", "r");
	int n_bc = 0;
	char **list_barcodes = NULL;
	char *bc = calloc(19, 1);
	int fre;
	int *arr_fre = NULL;


	while (fscanf(f, "%s", bc) != EOF) {
		fscanf(f, "%d", &fre);
//		printf("%s %d\n", bc, fre);
		list_barcodes = realloc(list_barcodes, (n_bc+1) * sizeof(char*));
		arr_fre = realloc(arr_fre, (n_bc+1) * sizeof(int));
		arr_fre[n_bc] = fre;
		list_barcodes[n_bc] = bc;
		n_bc++;
		bc = calloc(19, 1);
	}
	*n_barcodes = n_bc;
	*barcodes = list_barcodes;
	*freq = arr_fre;
}

inline int get_nu(const uint32_t *seq, int pos)
{
	return (seq[pos >> 4] >> ((pos & 15) << 1)) & 3;
}

inline void set_nu(uint32_t *seq, int pos, int val)
{
	seq[pos >> 4] |= val << ((pos & 15) << 1);
}

int concate_edge(struct asm_graph_t *g, int n_path, int *path, int ksize,
		uint32_t **res_seq)
{
	int total_len = ksize;
	for (int i = 0; i < n_path; ++i)
		total_len += g->edges[path[i]].seq_len - ksize;
	uint32_t *res = calloc((total_len + 15) / 16, 4);
	int cur_len = g->edges[path[0]].seq_len;
	for (int i = 0 ;  i < cur_len; i++) {
		set_nu(res, i, get_nu(g->edges[path[0]].seq, i));
	}

	for(int i = 1; i < n_path; i++) {
		struct asm_edge_t *a = &g->edges[path[i]];
		for (int j = ksize; j < a->seq_len; j++) {
			set_nu(res, cur_len + j-ksize, get_nu(a->seq, j));
		}
		cur_len += a->seq_len - ksize;
	}
	*res_seq = res;
	return 0;
}

//struct asm_graph_t* huu_create_new_graph(struct asm_graph_t *g, char *path)
//{
//	FILE *f = fopen(path, "r");
//	int a, b;
//	khash_t(long_spath) *spath = kh_init(long_spath);
//	get_all_shortest_paths_dp(g, spath);
//
//	struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
//	int count = 0 ;
//	FILE *out = fopen("have_path.txt","w");
//	while (fscanf(f, "%d %d\n", &a, &b) != EOF){
//		printf("get path from %d %d\n", a,b);
//		count++;
//		int *path = NULL;
//		int n_path;
//		int res = extract_shortest_path(g, spath, a, b, &path, &n_path);
//		if (res == -1) {
//			int t_b = g->edges[a].rc_id;
//			int t_a = g->edges[b].rc_id;
//			b = t_b;
//			a = t_a;
//			res = extract_shortest_path(g, spath, a, b, &path, &n_path);
//		}
//		if (res == -1)
//			continue;
//		fprintf(out, "%d %d\n", a,b);
////		assert(res != -1);
////		for(int i = 0 ; i < n_path; i++) {
////			printf("%d\n", path[i]);
////		}
////		printf("%d\n", count);
//		g0->edges = realloc(g0->edges, (g0->n_e+1) * sizeof(struct asm_edge_t));
//		int total_len;
//		g0->edges[g0->n_e].seq = concate_seq(g, n_path, path, g->ksize, &total_len);
//		g0->edges[g0->n_e].seq_len = total_len;
//		g0->edges[g0->n_e].count = g->edges[a].count + g->edges[b].count;
//		g0->edges[g0->n_e].n_holes = 0;
//		g0->n_e++;
//		free(path);
//	}
//	fclose(out);
//	for (int i = 0; i < g->n_e; i++) {
//		g0->edges = realloc(g0->edges, (g0->n_e+1) * sizeof(struct asm_edge_t));
//		g0->edges[g0->n_e].seq = g->edges[i].seq;
//		g0->edges[g0->n_e].seq_len = g->edges[i].seq_len;
//		g0->edges[g0->n_e].count = g->edges[i].count;
//		g0->edges[g0->n_e].n_holes = 0;
//		g0->n_e++;
//	}
//	return g0;
//}

int compare(const void *a, const void *b)
{
	int *x = *(int (*)[3]) a;
	int *y= *(int (*)[3]) b;
	if (x[0] < y[0] || (x[0] == y[0]  && x[1] < y[1])) return -1;
	if (x[0] == y[0]  && x[1] == y[1]) return 0;
	return 1;
}

void filter_xxx(struct asm_graph_t *g)
{
	FILE *in = fopen("/home/che/bioturing/data/yeast/metadata/all_shortest_paths.txt", "r");
	FILE *out = fopen("res.txt" ,"w");
	int a, b;
	char *tmp = calloc(1000,1), n_tmp;
	fprintf(out, "a,b\n");
	while (fscanf(in, "%d to %d:\n", &a, &b) != EOF) {
		fgets(tmp, 10000, in);
		if (g->edges[a].seq_len > 1000 && g->edges[b].seq_len > 1000) {
			fprintf(out, "%d,%d\n",a,b);
		}
		fflush(out);
	}
	fclose(out);
	fclose(in);
}

void dirty(struct asm_graph_t *g, struct opt_proc_t *opt)
{
	char log_path[1024];
	sprintf(log_path, "%s/get_long_contig.log", opt->out_dir);
	init_logger(opt->log_level, log_path);
	log_info("Version: %s", GIT_SHA);
	set_log_stage("Get long contig");
	get_list_contig(opt, g);
}

void test_sort_read(struct read_path_t *read_sorted_path, struct asm_graph_t *g)
{
	khash_t(bcpos) *dict = kh_init(bcpos);
	construct_read_index(read_sorted_path, dict);
	uint64_t *bc = calloc(5000, sizeof(uint64_t));
	int n_bc = 0;
	for (int i = 0; i < g->n_e; i++)
		if (g->edges[i].barcodes[2].n_item > 3000) {
			struct barcode_hash_t *buck = &g->edges[i].barcodes[2];
			for (int l = 0; (uint32_t) l < buck->size; l++) {
				if (buck->keys[l] != (uint64_t) (-1)) {
					uint64_t barcode = buck->keys[l];
					bc[n_bc++] = barcode;
					if (n_bc > 3000) {
						break;
					}
				}
			}
			break;
		}
	test_same_barcode(read_sorted_path, dict,
	                  bc, n_bc);
}
