#include <errno.h>
#include <sys/stat.h>
#include "build_bridge.h"
#include "helper.h"
#include "sort_read.h"
#include "barcode_resolve2.h"
#include "kmer_hash.h"
#include "resolve.h"
#include "utils.h"
#include "log.h"
#include "io_utils.h"
#include "unit_test.h"
#include "barcode_builder.h"
#define MIN_PROCESS_COV 500
#define SYNC_KEEP_GLOBAL 0
#define SYNC_KEEP_LOCAL 1
#define SYNC_MAX_GLOBAL 2
#define SYNC_MAX_LOCAL 3

void combine_edges(struct asm_graph_t lg, int *path, int path_len, char **seq)
{
	*seq = calloc(1, sizeof(char));
	for (int i = 0; i < path_len; ++i){
		char *cur_seq;
		decode_seq(&cur_seq, lg.edges[i].seq, lg.edges[i].seq_len);
		join_seq(seq, cur_seq + lg.ksize);
		free(cur_seq);
	}
}

/**
 * Brief: finds the id of the original contig e1 in the local graph
 * @param g: the global graph
 * @param lg: the local graph
 * @param e_id: the original contig id
 * @param emap: the mapping between the original edge and its counterpart
 * 	in the local graph
 */
int get_local_edge_head(struct asm_graph_t g, struct asm_graph_t lg,
		int e_id, struct edge_map_info_t *emap)
{
	int res;
	emap->gl_e = e_id;
	int *edge_id = &(emap->lc_e);
	struct subseq_pos_t *gpos = &(emap->gpos);
	struct subseq_pos_t *lpos = &(emap->lpos);
	struct asm_edge_t e = g.edges[e_id];
	struct map_contig_t mct;
	init_map_contig(&mct, g.edges[e.rc_id], lg);
	*edge_id = find_match(&mct);
	if (*edge_id == -1 || lg.edges[*edge_id].seq_len < WINDOW_SIZE)
		goto no_local_edge_found;
	get_match_pos(&mct, gpos, lpos);
	if (gpos->start > gpos->end || lpos->start > lpos->end)
		goto no_local_edge_found;
	*edge_id = lg.edges[*edge_id].rc_id;
	gpos->start = e.seq_len - gpos->start - WINDOW_SIZE;
	gpos->end = e.seq_len - gpos->end - WINDOW_SIZE;
	swap(&gpos->start, &gpos->end, sizeof(khint32_t));

	lpos->start = lg.edges[*edge_id].seq_len - lpos->start - WINDOW_SIZE;
	lpos->end = lg.edges[*edge_id].seq_len - lpos->end - WINDOW_SIZE;
	swap(&lpos->start, &lpos->end, sizeof(khint32_t));
	if (gpos->start > gpos->end || lpos->start > lpos->end)
		goto no_local_edge_found;
	if (gpos->start < 0 || gpos->end >= g.edges[emap->gl_e].seq_len)
		goto no_local_edge_found;
	if (lpos->start < 0 || lpos->end >= lg.edges[emap->lc_e].seq_len)
		goto no_local_edge_found;
	res = 1;
	goto end_function;
no_local_edge_found:
	log_debug("Mapping failed");
	res = 0;
	emap->lc_e = -1;
end_function:
	map_contig_destroy(&mct);
	return res;
}

/**
 * Brief: finds the id of the original contig e2 in the local graph
 * @param g: the global graph
 * @param lg: the local graph
 * @param e_id: the original contig id
 * @param emap: the mapping between the original edge and its counterpart
 * 	in the local graph
 */
int get_local_edge_tail(struct asm_graph_t g, struct asm_graph_t lg,
		int e_id, struct edge_map_info_t *emap)
{
	int res;
	emap->gl_e = e_id;
	int *edge_id = &(emap->lc_e);
	struct subseq_pos_t *gpos = &(emap->gpos);
	struct subseq_pos_t *lpos = &(emap->lpos);
	struct asm_edge_t e = g.edges[e_id];
	struct map_contig_t mct;
	init_map_contig(&mct, e, lg);
	*edge_id = find_match(&mct);
	if (*edge_id == -1 || lg.edges[*edge_id].seq_len < WINDOW_SIZE)
		goto no_local_edge_found;
	get_match_pos(&mct, gpos, lpos);
	if (gpos->start > gpos->end || lpos->start > lpos->end)
		goto no_local_edge_found;
	if (gpos->start < 0 || gpos->end >= g.edges[emap->gl_e].seq_len)
		goto no_local_edge_found;
	if (lpos->start < 0 || lpos->end >= lg.edges[emap->lc_e].seq_len)
		goto no_local_edge_found;
	res = 1;
	goto end_function;
no_local_edge_found:
	log_debug("Mapping failed");
	res = 0;
	emap->lc_e = -1;
end_function:
	map_contig_destroy(&mct);
	return res;
}

void sync_global_local_edge(struct asm_edge_t global, struct asm_edge_t local,
		struct subseq_pos_t global_pos, struct subseq_pos_t local_pos,
		int sync_type, char **res_seq)
{
	char *global_seq, *local_seq;
	decode_seq(&global_seq, global.seq, global.seq_len);
	decode_seq(&local_seq, local.seq, local.seq_len);
	if (sync_type == SYNC_KEEP_GLOBAL){
		int len = global_pos.start + local.seq_len - local_pos.start;
		*res_seq = (char *) calloc(len + 1, sizeof(char));
		strncpy(*res_seq, global_seq, global_pos.start);
		strcpy(*res_seq + global_pos.start, local_seq + local_pos.start);
	} else if (sync_type == SYNC_KEEP_LOCAL){
		int len = local_pos.end + global.seq_len - global_pos.end;
		*res_seq = (char *) calloc(len + 1, sizeof(char));
		strncpy(*res_seq, local_seq, local_pos.end);
		strcpy(*res_seq + local_pos.end, global_seq + global_pos.end);
	} else if (sync_type == SYNC_MAX_GLOBAL){
		int len = global_pos.start + local_pos.end - local_pos.start
			+ max(local.seq_len - local_pos.end,
					global.seq_len - global_pos.end);
		*res_seq = calloc(len + 1, sizeof(char));
		strncat(*res_seq, global_seq, global_pos.start);
		strncat(*res_seq, local_seq + local_pos.start,
				local_pos.end - local_pos.start);
		if (global.seq_len - global_pos.end > local.seq_len - local_pos.end)
			strcat(*res_seq, global_seq + global_pos.end);
		else
			strcat(*res_seq, local_seq + local_pos.end);
	} else {
		int len = global.seq_len - global_pos.end + local_pos.end - local_pos.start
			+ max(local_pos.start, global_pos.start);
		*res_seq = calloc(len + 1, sizeof(char));
		if (global_pos.start > local_pos.start)
			strncat(*res_seq, global_seq, global_pos.start);
		else
			strncat(*res_seq, local_seq, local_pos.start);
		strncat(*res_seq, local_seq + local_pos.start,
				local_pos.end - local_pos.start);
		strcat(*res_seq, global_seq + global_pos.end);
	}
	free(global_seq);
	free(local_seq);
}

void unrelated_filter(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2,
		int *scaffolds, int n_scaff, struct asm_graph_t *lg)
{
	log_debug("Filter irrelevant edges");
	log_debug("Before filter: %d edges", lg->n_e);
	int e1 = emap1->gl_e;
	int e2 = emap2->gl_e;
	int *bad = (int *) calloc(lg->n_e, sizeof(int));
	for (int i = 0; i < n_scaff; ++i){
		struct map_contig_t mct;
		init_map_contig(&mct, g->edges[scaffolds[i]], *lg);
		find_match(&mct);
		for (int j = 0; j < lg->n_e; ++j){
			int rc_id = lg->edges[j].rc_id;
			bad[j] |= mct.is_match[j] || mct.is_match[rc_id];
		}
		map_contig_destroy(&mct);
	}
	struct map_contig_t mct_1;
	init_map_contig(&mct_1, g->edges[g->edges[e1].rc_id], *lg);
	int lc_e1 = lg->edges[find_match(&mct_1)].rc_id;

	struct map_contig_t mct_2;
	init_map_contig(&mct_2, g->edges[e2], *lg);
	int lc_e2 = find_match(&mct_2);

	bad[lc_e1] = bad[lc_e2] = bad[lg->edges[lc_e1].rc_id]
		= bad[lg->edges[lc_e2].rc_id] = 0;
	for (int i = 0; i < lg->n_e; ++i){
		if (bad[i])
			asm_remove_edge(lg, i);
	}
	char tmp_name[1024];
	sprintf(tmp_name, "%s/tmp_graph_%d_%d.bin", opt->out_dir, emap1->gl_e,
			emap2->gl_e);
	struct asm_graph_t lg1;
	struct asm_graph_t g_bak;
	asm_clone_graph(lg, &g_bak, tmp_name);
	asm_condense(lg, &lg1);
	if (check_degenerate_graph(g, &lg1, emap1->gl_e, emap2->gl_e)){
		log_debug("Condensed graph degenerated, aborting filtering!");
		asm_graph_destroy(lg);
		asm_graph_destroy(&lg1);
		*lg = g_bak;
	} else {
		asm_graph_destroy(lg);
		asm_graph_destroy(&g_bak);
		*lg = lg1;
		log_debug("After filter: %d edges", lg1.n_e);
		get_local_edge_head(*g, lg1, emap1->gl_e, emap1);
		get_local_edge_tail(*g, lg1, emap2->gl_e, emap2);
		print_log_edge_map(emap1, emap2);
	}
	free(bad);
	map_contig_destroy(&mct_1);
	map_contig_destroy(&mct_2);
}

/**
 * @brief: trys to fill the gap between the two bridge contigs by doing
 * 	local assembly
 * @param opt: application options
 * @param g: the global graph
 * @param lg: the local graph
 * @param e1, e2: the two bridge contigs
 * @param scaffolds: list of contigs that are in the same scaffold path of e1 and e2
 * @param n_scaff: size of scaffolds
 * @param res_seq, seq_len: the result sequence and its length
 * @return: the result of the local assembly (BRIDGE_LOCAL_NOT_FOUND,
 * 	BRIDGE_TRIVIAL_BRIDGE, BRIDGE_PATH_NOT_FOUND, BRIDGE_MULTIPLE_PATHS_FOUND)
 * @description:
 * 	Given the global graph, the local graph, the bridge edges, the function
 * 	needs to find a sequence that bridge from one edge to the other
 * 		+ First, it needs to map the bridge edges in the global graph
 * 			to the local graph
 * 		+ Then, it tries to find a path between the mapped edges
 * 			in the local graph
 *
 * @potential bugs:
 * 	Uncorrect mapping might cause errors
 *
 * @post-conditions:
 * 	+ the function needs to always output a sequence
 * 	+ res must be either BRIDGE_LOCAL_NOT_FOUND, BRIDGE_TRIVIAL_BRIDGE, BRIDGE_PATH_NOT_FOUND
 * 	or BRIDGE_MULTIPLE_PATHS_FOUND
 */
int get_bridge(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, int e1, int e2, int *scaffolds,
		int n_scaff, char **res_seq, int *seq_len)
{
	struct edge_map_info_t emap1;
	int res_head = get_local_edge_head(*g, *lg, e1, &emap1);


	struct edge_map_info_t emap2;
	int res_tail = get_local_edge_tail(*g, *lg, e2, &emap2);

	print_log_edge_map(&emap1, &emap2);
	int res = try_bridging(opt, g, lg, scaffolds, n_scaff, &emap1, &emap2,
			res_seq, seq_len);
	post_test(get_bridge, *res_seq, *seq_len, res);
	return res;
}

void print_log_edge_map(struct edge_map_info_t *emap1, struct edge_map_info_t *emap2)
{
	log_trace("Local edge 1: %d", emap1->lc_e);
	log_trace("Global edge starts from: %d, ends at: %d",
			emap1->gpos.start, emap1->gpos.end);
	log_trace("Local edge starts from: %d, ends at: %d",
			emap1->lpos.start, emap1->lpos.end);
	log_trace("Local edge 2: %d", emap2->lc_e);
	log_trace("Global edge starts from: %d, ends at: %d",
			emap2->gpos.start, emap2->gpos.end);
	log_trace("Local edge starts from: %d, ends at: %d",
			emap2->lpos.start, emap2->lpos.end);
}

void join_complex_path(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_edge_t lc_e1, struct asm_edge_t lc_e2,
		struct subseq_pos_t gpos1, struct subseq_pos_t lpos1,
		struct subseq_pos_t gpos2, struct subseq_pos_t lpos2,
		char **res_seq)
{
	char *first, *second;
	sync_global_local_edge(e1, lc_e1, gpos1, lpos1, SYNC_KEEP_GLOBAL,
			&first);
	sync_global_local_edge(e2, lc_e2, gpos2, lpos2, SYNC_KEEP_LOCAL,
			&second);
	char *tmp;
	get_dump_N(&tmp);
	*res_seq = (char *) calloc(1, sizeof(char));
	join_seq(res_seq, first);
	join_seq(res_seq, tmp);
	join_seq(res_seq, second);
	free(first);
	free(second);
	free(tmp);
}

void join_middle_edge(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_edge_t lc_e1, struct asm_edge_t lc_e2,
		struct subseq_pos_t gpos1, struct subseq_pos_t lpos1,
		struct subseq_pos_t gpos2, struct subseq_pos_t lpos2,
		struct asm_edge_t middle, char **res_seq)
{
	char *first, *second, *mid;
	sync_global_local_edge(e1, lc_e1, gpos1, lpos1, SYNC_KEEP_GLOBAL,
			&first);
	sync_global_local_edge(e2, lc_e2, gpos2, lpos2, SYNC_KEEP_LOCAL,
			&second);
	decode_seq(&mid, middle.seq, middle.seq_len);
	char *dump_N;
	get_dump_N(&dump_N);
	*res_seq = (char *) calloc(1, sizeof(char));
	join_seq(res_seq, first);
	join_seq(res_seq, dump_N);
	join_seq(res_seq, mid);
	join_seq(res_seq, dump_N);
	join_seq(res_seq, second);
	free(first);
	free(second);
	free(dump_N);
	free(mid);

}

int try_bridging(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, int *scaffolds, int n_scaff,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2,
		char **res_seq, int *seq_len)
{
	int bridge_type;
	char *bridge_seq;
	int e1 = emap1->gl_e;
	int e2 = emap2->gl_e;
	int lc_e1 = emap1->lc_e;
	int lc_e2 = emap2->lc_e;
	struct subseq_pos_t gpos1 = emap1->gpos;
	struct subseq_pos_t gpos2 = emap2->gpos;
	struct subseq_pos_t lpos1 = emap1->lpos;
	struct subseq_pos_t lpos2 = emap2->lpos;
	if (lc_e1 == -1 || lc_e2 == -1){
		bridge_type = BRIDGE_LOCAL_NOT_FOUND;
		join_bridge_dump(g->edges[e1], g->edges[e2], &bridge_seq);
		goto end_function;
	} else if (lc_e1 == lc_e2){
		bridge_type = BRIDGE_TRIVIAL_BRIDGE;
		join_trivial_bridge(g->edges[e1], g->edges[e2], *lg, emap1,
				emap2, &bridge_seq);
		goto end_function;
	} else {
		int *path;
		int path_len;
		get_best_path(opt, g, lg, emap1, emap2, scaffolds, n_scaff,
				&path, &path_len);
		if (path_len == 0){
			bridge_type = BRIDGE_PATH_NOT_FOUND;
			// Graph is changed so the local edge id might not stay
			// the same
			join_bridge_no_path(g, lg, emap1, emap2, &bridge_seq);
			goto end_function;
		} else {
			bridge_type = BRIDGE_MULTIPLE_PATH;
			join_bridge_by_path(g->edges[emap1->gl_e], g->edges[emap2->gl_e],
					*lg, path, path_len, emap1, emap2, &bridge_seq);
			free(path);
			goto end_function;
		}
	}
end_function:
	*seq_len = strlen(bridge_seq);
	*res_seq = bridge_seq;
	return bridge_type;
}

/**
 * Finds the best path from the first local edge to the second one
 * @param opt: application options
 * @param g, lg: the global and local graph
 * @param emap1, emap2: mapping betwwen the global edges and local edges
 * @param scaffolds: the list of contigs that are on the same scaffold
 * 	path of e1 and e2
 * @param n_scaff: the size of scaffolds
 * @param path: the path that is found
 * @param path_len: the size of path
 */
void get_best_path(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, struct edge_map_info_t *emap1,
		struct edge_map_info_t *emap2, int *scaffolds, int n_scaff,
		int **path, int *path_len)
{
	//		 STAGE 1 getting local reads 			//
	*path = NULL;
	*path_len = 0;
	int e1 = emap1->gl_e;
	int e2 = emap2->gl_e;
	int lc_e1 = emap1->lc_e;
	int lc_e2 = emap2->lc_e;
	struct asm_edge_t edge_pre_e1;
	struct asm_edge_t edge_next_e2;

	struct read_path_t local_read_path;
	int ret = get_reads_kmer_check(opt, g, e1, e2, &local_read_path);
	if (!ret){
		log_warn("Something supicious happens, probably that read files are empty, please check at %s and %s",
				local_read_path.R1_path, local_read_path.R2_path);
		goto ignore_stage_1;
	}

	//		 STAGE 2 filtering local graph 			//
	unrelated_filter(opt, g, emap1, emap2, scaffolds, n_scaff, lg);
	connection_filter(opt, g, lg, emap1, emap2);


	// BUGGY CODE
	/*__VERBOSE("resolve local loop\n");
	asm_resolve_local_loop(lg);
	get_local_edge_head(*g, *lg, emap1->gl_e, emap1);
	get_local_edge_tail(*g, *lg, emap2->gl_e, emap2);
	print_log_edge_map(emap1, emap2);
	__VERBOSE("DONE\n");*/
	//link_filter(opt, g, lg, emap1, emap2);
	print_graph(opt, lg, emap1->gl_e, emap2->gl_e);

	// 		 STAGE 3 finding paths 			//
	log_info("Start finding paths");
	struct path_info_t pinfo;
	path_info_init(&pinfo);
	khash_t(kmer_int) *kmer_count = get_kmer_hash(local_read_path.R1_path,
			local_read_path.R2_path, KSIZE_CHECK);
	get_all_paths_kmer_check(lg, emap1, emap2, &pinfo, KSIZE_CHECK,
			kmer_count);
	kh_destroy(kmer_int, kmer_count);

	if (pinfo.n_paths == 0)
		goto ignore_stage_2;
	log_debug("Found %d paths, finding the best one", pinfo.n_paths);
	float *scores;
	float *error;
	get_path_scores(opt, &local_read_path, g, lg, &pinfo, e1, e2, &scores,
			&error);
	float best_score = 0;
	int best_path = 0;
	int min_score = 1e9;
	int max_err = 0;
	for (int i = 0; i < pinfo.n_paths; ++i){
		min_score = min(min_score, scores[i]);
		max_err = max(max_err, error[i]);
	}
	for (int i = 0; i < pinfo.n_paths; ++i){
		if (scores[i] - min_score + max_err - error[i]  > best_score){
			best_path = i;
			best_score = scores[i] - min_score + max_err - error[i];
		}
	}
	log_debug("Found best path id: %d, scores: %.3f\n",
			best_path, best_score);
	*path_len = pinfo.path_lens[best_path];
	*path = (int *) calloc(*path_len, sizeof(int));
	memcpy(*path, pinfo.paths[best_path], sizeof(int) * *path_len);
	free(scores);
	free(error);
ignore_stage_2:
	path_info_destroy(&pinfo);
ignore_stage_1:
	destroy_read_path(&local_read_path);
}

void get_path_scores(struct opt_proc_t *opt, struct read_path_t *local_read_path,
		struct asm_graph_t *g, struct asm_graph_t *lg,
		struct path_info_t *pinfo, int e1, int e2, float **scores,
		float **error)
{
	char cand_path[1024];
	sprintf(cand_path, "%s/%d_%d_all.fasta", opt->out_dir, e1, e2);
	FILE *f = fopen(cand_path, "w");
	for (int i = 0; i < pinfo->n_paths; ++i){
		fprintf(f, ">%d\n", i);
		char *seq;
		join_bridge_center_by_path(lg, pinfo->paths[i],
				pinfo->path_lens[i], &seq);
		fprintf(f, "%s\n", seq);
		free(seq);
	}
	fclose(f);
	khash_t(contig_count) *ctg_cnt = kh_init(contig_count);
	khash_t(contig_count) *count_err = kh_init(contig_count);
	for (int i = 0; i < pinfo->n_paths; ++i){
		int ret;
		khiter_t it = kh_put(contig_count, ctg_cnt, i, &ret);
		kh_val(ctg_cnt, it) = 0;

		it = kh_put(contig_count, count_err, i, &ret);
		kh_val(count_err, it) = 0;
	}

	count_readpair_err_path(opt->n_threads, local_read_path, cand_path,
			ctg_cnt, count_err);
	*scores = (float *) calloc(pinfo->n_paths, sizeof(float));
	for (khiter_t it = kh_begin(ctg_cnt); it != kh_end(ctg_cnt); ++it){
		if (!kh_exist(ctg_cnt, it))
			continue;
		int key = kh_key(ctg_cnt, it);
		int val = kh_val(ctg_cnt, it);
		(*scores)[key] = val;
	}

	*error = (float *) calloc(pinfo->n_paths, sizeof(float));
	for (khiter_t it = kh_begin(count_err); it != kh_end(count_err); ++it){
		if (!kh_exist(count_err, it))
			continue;
		int key = kh_key(count_err, it);
		int val = kh_val(count_err, it);
		(*error)[key] = val;
	}
	kh_destroy(contig_count, ctg_cnt);
	kh_destroy(contig_count, count_err);
}

void join_seq(char **dest, char *source)
{
	int old_len = strlen(*dest);
	int new_len = old_len + strlen(source);
	*dest = (char *) realloc(*dest, (new_len + 1) * sizeof(char));
	strcpy(*dest + old_len, source);
}

/**
 * @brief: joins the two bridge contigs by the local assembly result
 * @param e1, e2: the two bridge contigs
 * @param lg: the local graph
 * @param emap1, emap2: mapping between the original edges and the local edges
 * @param res_seq: the result seq
 */
void join_trivial_bridge(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_graph_t lg, struct edge_map_info_t *emap1,
		struct edge_map_info_t *emap2, char **res_seq)
{
	int local_edge = emap1->lc_e;
	struct subseq_pos_t gpos1 = emap1->gpos;
	struct subseq_pos_t gpos2 = emap2->gpos;
	struct subseq_pos_t lpos1 = emap1->lpos;
	struct subseq_pos_t lpos2 = emap2->lpos;
	char *edge_seq1, *edge_seq2, *local_seq;
	decode_seq(&edge_seq1, e1.seq, e1.seq_len);
	decode_seq(&edge_seq2, e2.seq, e2.seq_len);
	decode_seq(&local_seq, lg.edges[local_edge].seq,
			lg.edges[local_edge].seq_len);

	if (lpos2.start < lpos1.end){
		int diff = lpos2.start - lpos1.end;
		lpos2.start = lpos1.end;
		gpos2.start += diff;
	}
	int len = gpos1.end + lpos2.start - lpos1.end + e2.seq_len
			- gpos2.start;
	*res_seq = (char *) calloc(len + 1, sizeof(char));
	len = 0;
	strncpy(*res_seq + len, edge_seq1, gpos1.end);
	len += gpos1.end;
	strncpy(*res_seq + len, local_seq + lpos1.end, lpos2.start - lpos1.end);
	len += lpos2.start - lpos1.end;
	strncpy(*res_seq + len, edge_seq2, e2.seq_len - gpos2.start);
	len += e2.seq_len - gpos2.start;
	free(edge_seq1);
	free(edge_seq2);
	free(local_seq);
}

/**
 * @brief: join the two bridge edges by the path that is found from one edge
 * 	to another
 * @param e1, e2: the two bridge edges
 * @param lg: the local graph
 * @param path, path_len: the found path and its length
 * @param emap1, emap2: the mapping between the global and local edges
 * @param res-seq: the result sequence
 */
void join_bridge_by_path(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_graph_t lg, int *path, int path_len,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2,
		char **res_seq)
{
	struct subseq_pos_t gpos1 = emap1->gpos;
	struct subseq_pos_t gpos2 = emap2->gpos;
	struct subseq_pos_t lpos1 = emap1->lpos;
	struct subseq_pos_t lpos2 = emap2->lpos;
	char *head_seq, *tail_seq;
	int lc_e1 = path[0];
	int lc_e2 = path[path_len - 1];
	log_debug("Joining from %d to %d", lc_e1, lc_e2);
	sync_global_local_edge(e1, lg.edges[lc_e1], gpos1, lpos1,
			SYNC_KEEP_GLOBAL, &head_seq);
	sync_global_local_edge(e2, lg.edges[lc_e2], gpos2, lpos2,
			SYNC_KEEP_LOCAL, &tail_seq);
	*res_seq = (char *) calloc(1, sizeof(char));
	join_seq(res_seq, head_seq);
	for (int i = 1; i < path_len - 1; ++i){
		char *tmp;
		struct asm_edge_t local_edge = lg.edges[path[i]];
		decode_seq(&tmp, local_edge.seq, local_edge.seq_len);
		join_seq(res_seq, tmp + lg.ksize);
		free(tmp);
	}
	join_seq(res_seq, tail_seq + lg.ksize);
	free(head_seq);
	free(tail_seq);
}

void join_bridge_center_by_path(struct asm_graph_t *lg, int *path, int path_len,
		char **seq)
{
	int len = lg->edges[path[0]].seq_len;
	for (int i = 1; i < path_len; ++i)
		len += lg->edges[path[i]].seq_len - lg->ksize;
	*seq = (char *) calloc(len + 1, sizeof(char));
	for (int i = 0; i < path_len; ++i){
		char *tmp;
		decode_seq(&tmp, lg->edges[path[i]].seq,
				lg->edges[path[i]].seq_len);
		int pos = i == 0 ? 0 : lg->ksize;
		strcat(*seq, tmp + pos);
		free(tmp);
	}
}

void join_bridge_no_path(struct asm_graph_t *g, struct asm_graph_t *lg,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2,
		char **res_seq)
{
	*res_seq = (char *) calloc(1, sizeof(char));
	char *first, *second;
	sync_global_local_edge(g->edges[emap1->gl_e], lg->edges[emap1->lc_e],
			emap1->gpos, emap1->lpos, SYNC_MAX_GLOBAL, &first);
	sync_global_local_edge(g->edges[emap2->gl_e], lg->edges[emap2->lc_e],
			emap2->gpos, emap2->lpos, SYNC_MAX_LOCAL, &second);
	char *dump_N;
	get_dump_N(&dump_N);
	join_seq(res_seq, first);
	join_seq(res_seq, dump_N);
	join_seq(res_seq, second);
	free(first);
	free(second);
	free(dump_N);
}

/**
 * @brief: simply joins the two contigs by filling Ns in the gap
 * @param e1, e2: the two contigs in the order that will be joined
 * @param res_seq: the result sequence
 */
void join_bridge_dump(struct asm_edge_t e1, struct asm_edge_t e2,
		char **res_seq)
{
	char *first, *second;
	decode_seq(&first, e1.seq, e1.seq_len);
	decode_seq(&second, e2.seq, e2.seq_len);
	char *dump_N;
	get_dump_N(&dump_N);
	*res_seq = (char *) calloc(1, sizeof(char));
	join_seq(res_seq, first);
	join_seq(res_seq, dump_N);
	join_seq(res_seq, second);
}

void cov_filter(struct asm_graph_t *g, struct asm_graph_t *lg,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2)
{
	log_info("Filter by coverage");
	log_debug("Before filter: %ld edges\n", lg->n_e);
	int thresh = (int) (MIN_DEPTH_RATIO *
			min(__get_edge_cov(lg->edges + emap1->lc_e, lg->ksize),
			__get_edge_cov(lg->edges + emap2->lc_e, lg->ksize)));
	for (int i = 0; i < lg->n_e; ++i){
		if (__get_edge_cov(lg->edges + i, lg->ksize) < thresh)
			asm_remove_edge(lg, i);
	}
	struct asm_graph_t lg1;
	struct asm_graph_t g_bak;
	char tmp_name[1024];
	sprintf(tmp_name, "tmp_graph_%d_%d.bin", emap1->gl_e, emap2->gl_e);
	asm_clone_graph(lg, &g_bak, tmp_name);
	asm_condense(lg, &lg1);
	if (check_degenerate_graph(g, &lg1, emap1->gl_e, emap2->gl_e)){
		log_debug("Condensed graph degenerated, aborting filtering!\n");
		asm_graph_destroy(lg);
		asm_graph_destroy(&lg1);
		*lg = g_bak;
	} else {
		asm_graph_destroy(lg);
		asm_graph_destroy(&g_bak);
		*lg = lg1;
		log_debug("After filter: %d edges\n", lg1.n_e);
		get_local_edge_head(*g, lg1, emap1->gl_e, emap1);
		get_local_edge_tail(*g, lg1, emap2->gl_e, emap2);
		print_log_edge_map(emap1, emap2);
	}
}

void connection_filter(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, struct edge_map_info_t *emap1,
		struct edge_map_info_t *emap2)
{
	log_info("Filter by connections");
	log_debug("Before filter: %d edges", lg->n_e);
	int *forward_len, *backward_len;
	struct graph_info_t ginfo;
	graph_info_init(lg, &ginfo, emap1->lc_e, emap2->lc_e);
	bfs(lg, &ginfo, emap1->lc_e, &forward_len);
	bfs(lg, &ginfo, lg->edges[emap2->lc_e].rc_id, &backward_len);
	graph_info_destroy(&ginfo);

	int *bad = (int *) calloc(lg->n_e, sizeof(int));
	for (int i = 0; i < lg->n_e; ++i){
		int l1 = forward_len[i];
		int l2 = backward_len[lg->edges[i].rc_id];
		if (l1 == -1 || l2 == -1 || l1 + l2 > MIN_PATH_LENGTH)
			bad[i] = 1;
	}
	for (int i = 0; i < lg->n_e; ++i)
		if (!bad[lg->edges[i].rc_id])
			bad[i] = 0;
	bad[emap1->lc_e] = bad[lg->edges[emap1->lc_e].rc_id]
		= bad[emap2->lc_e] = bad[lg->edges[emap2->lc_e].rc_id] = 0;
	for (int i = 0; i < lg->n_e; ++i){
		if (bad[i])
			asm_remove_edge(lg, i);
	}
	free(bad);
	struct asm_graph_t lg1;
	struct asm_graph_t g_bak;
	char tmp_name[1024];
	sprintf(tmp_name, "%s/tmp_graph_%d_%d.bin", opt->out_dir, emap1->gl_e,
			emap2->gl_e);
	asm_clone_graph(lg, &g_bak, tmp_name);
	asm_condense(lg, &lg1);
	if (check_degenerate_graph(g, &lg1, emap1->gl_e, emap2->gl_e)){
		log_debug("Condensed graph degenerated, aborting filtering!");
		asm_graph_destroy(lg);
		asm_graph_destroy(&lg1);
		*lg = g_bak;
	} else {
		asm_graph_destroy(lg);
		asm_graph_destroy(&g_bak);
		*lg = lg1;
		log_debug("After filter: %d edges", lg1.n_e);
		get_local_edge_head(*g, lg1, emap1->gl_e, emap1);
		get_local_edge_tail(*g, lg1, emap2->gl_e, emap2);
		print_log_edge_map(emap1, emap2);
	}
	free(forward_len);
	free(backward_len);
}

void link_filter(struct opt_proc_t *opt, struct asm_graph_t *g, struct asm_graph_t *lg,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2)
{
	char r1_path[1024], r2_path[1024];
	sprintf(r1_path, "%s/local_assembly_shared_%d_%d/R1.sub.fq", opt->out_dir,
			emap1->gl_e, emap2->gl_e);
	sprintf(r2_path, "%s/local_assembly_shared_%d_%d/R2.sub.fq", opt->out_dir,
			emap1->gl_e, emap2->gl_e);
	khash_t(kmer_int) *kmer_count = get_kmer_hash(r1_path, r2_path, KSIZE_CHECK);
	int watch = 32;
	for (int e1 = 0; e1 < lg->n_e; ++e1){
		int tg = lg->edges[e1].target;
		float max_score = 0;
		printf("edge %d\n", e1);
		for (int i = 0; i < lg->nodes[tg].deg; ++i){
			int e2 = lg->nodes[tg].adj[i];
			char *first, *second;
			decode_seq(&first, lg->edges[e1].seq, lg->edges[e1].seq_len);
			decode_seq(&second, lg->edges[e2].seq, lg->edges[e2].seq_len);
			int len1 = strlen(first);
			int len2 = strlen(second);
			int m = min(len1, KSIZE_CHECK);
			int pos = len1 - m;
			char *join = (char *) calloc(m + len2 - lg->ksize + 1, sizeof(char));
			strcpy(join, first + pos);
			strncpy(join + m, second + lg->ksize,
					min(len2, KSIZE_CHECK) - lg->ksize);
			float score = count_kmer_on_seq(kmer_count, join, KSIZE_CHECK);

			max_score = max(max_score, score);
			printf("e2 %d m %d score %.3f\n", e2, m, score);
			print_kmer_count_on_seq(kmer_count, join, KSIZE_CHECK);
			free(first);
			free(second);
			free(join);
		}
	}
	kh_destroy(kmer_int, kmer_count);
}

int check_degenerate_graph(struct asm_graph_t *g, struct asm_graph_t *lg,
		int e1, int e2)
{
	struct edge_map_info_t emap1 = {0};
	get_local_edge_head(*g, *lg, e1, &emap1);
	struct edge_map_info_t emap2 = {0};
	get_local_edge_head(*g, *lg, e2, &emap2);

	if (emap1.lc_e == emap2.lc_e || emap1.lc_e == -1 || emap2.lc_e == -1)
		return 1;
	return 0;
}

/**
 * @brief The main function for local assembly process
 * @param opt main options struct
 * @param f the final scaffolds file with gap closed.
 */
void build_bridge(struct opt_proc_t *opt)
{
	FILE *f = fopen(opt->lc, "w");
	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_file); /* The global assembly graph */
	test_asm_graph(g0);

	struct scaffold_record_t scaffolds;
	struct query_record_t query_record;
	get_scaffolds_info(opt, &scaffolds);
	get_local_assembly_query(&scaffolds, &query_record);

	int *mark = (int *) calloc(g0->n_e, sizeof(int));
	for (int i = 0; i < scaffolds.n_paths; ++i){
		int path_len = scaffolds.path_lens[i];
		for (int j = 0; j < path_len; ++j){
			mark[scaffolds.paths[i][j]] = 1;
			mark[g0->edges[scaffolds.paths[i][j]].rc_id] = 1;
		}
	}

	log_info("Done initializing scaffold paths\n");

	log_info("Getting all local graphs");
	get_all_local_graphs(opt, g0, &query_record); /* Iteratively build the local assembly graph */
	log_info("Done getting all local graphs");

	char **bridges = calloc(query_record.n_process, sizeof(char *));

	pthread_mutex_t query_lock;
	pthread_mutex_init(&query_lock, NULL);
	pthread_mutex_t bridge_lock;
	pthread_mutex_init(&bridge_lock, NULL);
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	pthread_mutex_t *bridge_process_locks = calloc(query_record.n_process,
			sizeof(pthread_mutex_t));
	for (int i = 0; i < query_record.n_process; ++i)
		pthread_mutex_init(bridge_process_locks + i, NULL);
	//opt->n_threads = 1; /* Edit to 1 thread */
	struct build_bridge_bundle_t *worker_bundles = calloc(opt->n_threads,
			sizeof(struct build_bridge_bundle_t));
	for (int i = 0; i < opt->n_threads; ++i){
		worker_bundles[i].opt = opt;
		worker_bundles[i].g = g0;
		worker_bundles[i].query_record = &query_record;
		worker_bundles[i].query_lock = &query_lock;
		worker_bundles[i].bridge_lock = &bridge_lock;
		worker_bundles[i].bridges = bridges;
		worker_bundles[i].scaffold_record = &scaffolds;
		worker_bundles[i].bridge_process_locks = bridge_process_locks;
	}
	log_info("Building bridges on scaffold");
	pthread_t *worker_threads = calloc(opt->n_threads, sizeof(pthread_t));
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, build_bridge_iterator,
				worker_bundles + i);
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);
	free(query_record.e1);
	free(query_record.e2);
	free(query_record.path_id);
	pthread_mutex_destroy(&query_lock);
	pthread_mutex_destroy(&bridge_lock);
	for (int i = 0; i < query_record.n_process; ++i)
		pthread_mutex_destroy(bridge_process_locks + i);
	free(bridge_process_locks);
	free(worker_bundles);
	free(worker_threads);

	log_info("Done local assembly. Now print all the bridged sequences");
	print_bridges(f, g0, &scaffolds, bridges);
	for (int i = 0; i < query_record.n_process; ++i)
		free(bridges[i]);
	free(bridges);

	log_info("Print remain sequences");
	for (int i = 0; i < g0->n_e; ++i){
		if (g0->edges[i].seq_len < MIN_OUTPUT_CONTIG_LEN)
			continue;
		if (mark[i] == 0){
			int rc = g0->edges[i].rc_id;
			char *tmp;
			decode_seq(&tmp, g0->edges[i].seq, g0->edges[i].seq_len);
			fprintf(f, ">%d_%d\n", i, rc);
			fprintf(f, "%s\n", tmp);
			free(tmp);
			mark[rc] = 1;
		}
	}
	free(mark);
	free(scaffolds.path_lens);
	for (int i = 0; i < scaffolds.n_paths; ++i)
		free(scaffolds.paths[i]);
	free(scaffolds.paths);

	asm_graph_destroy(g0);
	free(g0);

	cleanup(opt);
	fclose(f);
}

/**
 * @brief: a worker funtion for local assembly
 * @param data: a bundle to store all the stuffs needed for local assembly
 */
void *build_bridge_iterator(void *data)
{
	int local_asm_res;
	struct build_bridge_bundle_t *bundle = (struct build_bridge_bundle_t *)
			data;
	while (1){
		pthread_mutex_lock(bundle->query_lock);
		int process_pos = bundle->query_record->process_pos;
		++bundle->query_record->process_pos;
		pthread_mutex_unlock(bundle->query_lock);
		if (process_pos >= bundle->query_record->n_process)
			break;

		int e1 = bundle->query_record->e1[process_pos];
		int e2 = bundle->query_record->e2[process_pos];
		int scaf_id = -1;
		for (int i = 0; i < bundle->query_record->n_process; ++i){
			if (bundle->query_record->e1[i] == e1 &&
				bundle->query_record->e2[i] == e2){
				scaf_id = i;
				break;
			}
		}
		if (scaf_id == -1)
			log_error("Cannot find the index of the bridge in the scaffold path, probably something went wrong");
		pthread_mutex_lock(bundle->bridge_process_locks + scaf_id);

		int path_id = bundle->query_record->path_id[process_pos];
		int *scaffolds = bundle->scaffold_record->paths[path_id];
		int n_scaff = bundle->scaffold_record->path_lens[path_id];
		struct asm_graph_t *g = bundle->g;
		char *seq;
		int seq_len;
		if (__get_edge_cov(g->edges + e1, g->ksize) > MIN_PROCESS_COV
			|| __get_edge_cov(g->edges + e2, g->ksize) > MIN_PROCESS_COV){
			log_info("Graph is too complex, filling Ns");
			join_bridge_dump(g->edges[e1], g->edges[e2], &seq);
		} else {
			char graph_bin_path[1024];
			sprintf(graph_bin_path, "%s/local_assembly_%d_%d/graph_k_%d_local_lvl_1.bin",
					bundle->opt->out_dir,
					bundle->g->edges[e1].rc_id,
					e2, bundle->opt->lk);
			if (access(graph_bin_path, F_OK) == -1){
				log_info("Local graph not found, filling Ns");
				join_bridge_dump(g->edges[e1], g->edges[e2],
						&seq);
			} else {
				struct asm_graph_t lg;
				load_asm_graph(&lg, graph_bin_path);
				local_asm_res = get_bridge(bundle->opt, bundle->g, &lg,
						e1, e2, scaffolds, n_scaff, &seq,
						&seq_len);
				asm_graph_destroy(&lg);
			}
		}
		pthread_mutex_lock(bundle->bridge_lock);
		log_debug("Local assembly status for edges %d and %d: %s", e1, e2,
				local_asm_result[local_asm_res]);
		bundle->bridges[process_pos] = seq;
		pthread_mutex_unlock(bundle->bridge_lock);

		pthread_mutex_unlock(bundle->bridge_process_locks + scaf_id);
	}
}

/**
 * @brief Get all the local graphs for local assembly
 * @param opt: options
 * @param g: the original graph (global graph)
 * @param query: a list of bridges pairs
 * @description:
 * 	Given a list of bridges pairs, the function needs to build all the local graphs
 * 	on those pairs.
 * 	For each pairs, it needs to:
 * 		+ Get the reads in that region
 * 		+ Use KMC to get kmer table
 * 		+ Build graph level 0 from the kmer table
 * 		+ Resolve 0-1
 *
 * @preconditions:
 * 	Original read files must be sorted
 */
void get_all_local_graphs(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct query_record_t *query)
{
	pre_test(get_all_local_graphs, opt);
	struct read_path_t read_sorted_path = parse_read_path_from_opt(opt);
	khash_t(bcpos) *dict = kh_init(bcpos);
	construct_read_index(&read_sorted_path, dict);

	for (int i = 0; i < query->n_process; ++i){
		log_info("Processing %d on %d local graphs", i, query->n_process);
		int e1 = query->e1[i];
		int e2 = query->e2[i];
		log_debug("Bridge contigs: %d %d", e1, e2);
		/*
		 * Build the local assembly graph for two edges: e1.rev and e2
		 */
		if (__get_edge_cov(g->edges + e1, g->ksize) > MIN_PROCESS_COV
			|| __get_edge_cov(g->edges + e2, g->ksize) > MIN_PROCESS_COV){
			log_debug("Too complex region, continue");
			continue;
		}
		get_local_assembly(opt, g, g->edges[e1].rc_id, e2, dict);
	}
	log_info("All of the local assembly graph are constructed. Now trying to bridging each pair of edges");
	kh_destroy(bcpos, dict);
}

void cleanup(struct opt_proc_t *opt)
{
	if (opt->log_level <= LOG_DEBUG_TECH){
		log_info("Currently in technical debug mode, nothing to be cleaned");
		return;
	}
	log_info("Cleaning up all files in %s", opt->out_dir);
	int flag = 0;
	flag = recursive_delete(opt->out_dir);
	if (flag)
		log_info("Some files are not deleted, please check in %s",
				opt->out_dir);
	log_info("Done cleaning up");
}

/**
 * @brief: gets the scaffolds paths info
 * @param opt: application options
 * @param scaffolds: a structure to store the scaffolds paths
 */
void get_scaffolds_info(struct opt_proc_t *opt, struct scaffold_record_t *scaffolds)
{
	FILE *fp = xfopen(opt->in_fasta, "r");
	int n_paths;
	fscanf(fp, "%d\n", &n_paths); /* Total number of paths in the scaffolds */
	scaffolds->n_paths = n_paths;
	scaffolds->path_lens = calloc(n_paths, sizeof(int));
	scaffolds->paths = calloc(n_paths, sizeof(int *));
	/*
	 * Prepare the bundle structs for local assembly
	 * One local assembly process needs an e1, e2, pre e1, next e2
	 */
	for (int i = 0; i < n_paths; ++i){
		int path_len;
		fscanf(fp, "%d\n", &path_len);
		scaffolds->path_lens[i] = path_len;
		scaffolds->paths[i] = calloc(path_len, sizeof(int));
		for (int j = 0; j < path_len; ++j)
			fscanf(fp, "%d", scaffolds->paths[i] + j);
	}
	fclose(fp);
}

void get_local_assembly_query(struct scaffold_record_t *scaffolds,
		struct query_record_t *query)
{
	query->e1 = NULL;
	query->e2 = NULL;
	query->path_id = NULL;
	query->process_pos = 0;
	query->n_process = 0;
	for (int i = 0; i < scaffolds->n_paths; ++i){
		int path_len = scaffolds->path_lens[i];
		query->n_process += path_len - 1;
		query->e1 = realloc(query->e1, query->n_process
				* sizeof(int));
		query->e2 = realloc(query->e2, query->n_process
				* sizeof(int));
		query->path_id = realloc(query->path_id, query->n_process
				* sizeof(int));
		for (int j = 1; j < path_len; ++j){
			query->e1[query->process_pos] = scaffolds->paths[i][j - 1];
			query->e2[query->process_pos] = scaffolds->paths[i][j];
			query->path_id[query->process_pos] = i;
			++query->process_pos;
		}
	}
	query->process_pos = 0;
}

void print_bridges(FILE *f, struct asm_graph_t *g, struct scaffold_record_t *scaffolds,
		char **bridges)
{
	int p = 0;
	int **paths = scaffolds->paths;
	int *path_lens = scaffolds->path_lens;
	for (int i = 0; i < scaffolds->n_paths; ++i){
		fprintf(f, ">contig_%d\n", i);
		char *seq;
		decode_seq(&seq, g->edges[paths[i][0]].seq,
				g->edges[paths[i][0]].seq_len);
		fprintf(f, "%s", seq);
		for (int j = 1; j < path_lens[i]; ++j){
			fprintf(f, "%s", bridges[p] +
					g->edges[paths[i][j - 1]].seq_len);
			++p;
		}
		free(seq);
		fprintf(f, "\n");
	}
}
