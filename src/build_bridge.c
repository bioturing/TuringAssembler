#include "build_bridge.h"
#include "helper.h"

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

void get_local_edge_head(struct asm_graph_t g, struct asm_graph_t lg,
		int e_id, struct edge_map_info_t *emap)
{
	emap->gl_e = e_id;
	int *edge_id = &(emap->lc_e);
	struct subseq_pos_t *gpos = &(emap->gpos);
	struct subseq_pos_t *lpos = &(emap->lpos);
	struct asm_edge_t e = g.edges[e_id];
	struct map_contig_t mct;
	init_map_contig(&mct, g.edges[e.rc_id], lg);
	*edge_id = find_match(&mct);
	if (*edge_id == -1)
		goto no_local_edge_found;
	get_match_pos(&mct, gpos, lpos);
	*edge_id = lg.edges[*edge_id].rc_id;
	gpos->start = e.seq_len - gpos->start - WINDOW_SIZE;
	gpos->end = e.seq_len - gpos->end - WINDOW_SIZE;
	swap(&gpos->start, &gpos->end, sizeof(khint32_t));

	lpos->start = lg.edges[*edge_id].seq_len - lpos->start - WINDOW_SIZE;
	lpos->end = lg.edges[*edge_id].seq_len - lpos->end - WINDOW_SIZE;
	swap(&lpos->start, &lpos->end, sizeof(khint32_t));
no_local_edge_found:
	map_contig_destroy(&mct);
}

void get_local_edge_tail(struct asm_graph_t g, struct asm_graph_t lg,
		int e_id, struct edge_map_info_t *emap)
{
	emap->gl_e = e_id;
	int *edge_id = &(emap->lc_e);
	struct subseq_pos_t *gpos = &(emap->gpos);
	struct subseq_pos_t *lpos = &(emap->lpos);
	struct asm_edge_t e = g.edges[e_id];
	struct map_contig_t mct;
	init_map_contig(&mct, e, lg);
	*edge_id = find_match(&mct);
	if (*edge_id == -1)
		goto no_local_edge_found;
	get_match_pos(&mct, gpos, lpos);
no_local_edge_found:
	map_contig_destroy(&mct);
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
	} else {
		int len = local_pos.end + global.seq_len - global_pos.end;
		*res_seq = (char *) calloc(len + 1, sizeof(char));
		strncpy(*res_seq, local_seq, local_pos.end);
		strcpy(*res_seq + local_pos.end, global_seq + global_pos.end);
	}
	free(global_seq);
	free(local_seq);
}

void unrelated_filter(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_edge_t pre_e1, struct asm_edge_t next_e2,
		struct asm_graph_t *lg, struct graph_info_t *ginfo)
{
	__VERBOSE_LOG("UNRELATED FILTER", "+++++++++++++++++++++++++++\n");
	__VERBOSE_LOG("", "Before filter: %d edges\n", lg->n_e);
	int *bad = (int *) calloc(lg->n_e, sizeof(int));
	struct map_contig_t mct_1;
	init_map_contig(&mct_1, e1, *lg);
	int lc_e1 = find_match(&mct_1);

	struct map_contig_t mct_2;
	init_map_contig(&mct_2, e2, *lg);
	int lc_e2 = find_match(&mct_2);

	for (int i = 0; i < lg->n_e; ++i){
		int rc_id = lg->edges[i].rc_id;
		bad[i] |= mct_1.is_match[i] || mct_1.is_match[rc_id]
			|| mct_2.is_match[i] || mct_2.is_match[rc_id];
	}
	if (pre_e1.source != -1){
		struct map_contig_t mct_3;
		init_map_contig(&mct_3, pre_e1, *lg);
		find_match(&mct_3);
		
		for (int i = 0; i < lg->n_e; ++i){
			int rc_id = lg->edges[i].rc_id;
			bad[i] |= mct_3.is_match[i] || mct_3.is_match[rc_id];
		}
		map_contig_destroy(&mct_3);
	}
	if (next_e2.source != -1){
		struct map_contig_t mct_4;
		init_map_contig(&mct_4, next_e2, *lg);
		find_match(&mct_4);
		
		for (int i = 0; i < lg->n_e; ++i){
			int rc_id = lg->edges[i].rc_id;
			bad[i] |= mct_4.is_match[i] || mct_4.is_match[rc_id];
		}
		map_contig_destroy(&mct_4);
	}
	bad[lc_e1] = bad[lc_e2] = bad[lg->edges[lc_e1].rc_id]
		= bad[lg->edges[lc_e2].rc_id] = 0;
	int is_disabled = 0;
	for (int i = 0; i < lg->n_e; ++i){
		if (bad[i]){
			mark_edge_trash(ginfo, i);
			++is_disabled;
		}
	}
	__VERBOSE_LOG("", "After filter: %d edges\n", lg->n_e - is_disabled);
	free(bad);
	map_contig_destroy(&mct_1);
	map_contig_destroy(&mct_2);
}

int get_bridge(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, int e1, int e2, int pre_e1, int next_e2,
		char **res_seq, int *seq_len)
{
	__VERBOSE("Matching edges...\n");
	struct edge_map_info_t emap1;
	get_local_edge_head(*g, *lg, e1, &emap1);

	__VERBOSE_LOG("", "Local edge 1: %d\n", emap1.lc_e);
	__VERBOSE_LOG("", "Global edge starts from: %d, ends at: %d\n",
			emap1.gpos.start, emap1.gpos.end);
	__VERBOSE_LOG("", "Local edge starts from: %d, ends at: %d\n",
			emap1.lpos.start, emap1.lpos.end);

	struct edge_map_info_t emap2;
	get_local_edge_tail(*g, *lg, e2, &emap2);
	__VERBOSE_LOG("", "Local edge 2: %d\n", emap2.lc_e);
	__VERBOSE_LOG("", "Global edge starts from: %d, ends at: %d\n",
			emap2.gpos.start, emap2.gpos.end);
	__VERBOSE_LOG("", "Local edge starts from: %d, ends at: %d\n",
			emap2.lpos.start, emap2.lpos.end);

	int res = try_bridging(opt, g, lg, pre_e1, next_e2, &emap1, &emap2,
			res_seq, seq_len);
	return res;
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
		struct asm_graph_t *lg, int pre_e1, int next_e2,
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
		get_best_path(opt, g, lg, emap1, emap2, pre_e1, next_e2,
				&path, &path_len);
		if (path_len == 0){
			bridge_type = BRIDGE_PATH_NOT_FOUND;
			join_bridge_no_path(g->edges[e1], g->edges[e2],
					g->edges[lc_e1], g->edges[lc_e2],
					emap1, emap2, &bridge_seq);
			goto end_function;
		} else {
			bridge_type = BRIDGE_MULTIPLE_PATH;
			join_bridge_by_path(g->edges[e1], g->edges[e2], *lg,
					path, path_len, emap1, emap2,
					&bridge_seq);
			free(path);
			goto end_function;
		}
	}
end_function:
	*seq_len = strlen(bridge_seq);
	*res_seq = bridge_seq;
	return bridge_type;
}

void get_best_path(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, struct edge_map_info_t *emap1,
		struct edge_map_info_t *emap2, int pre_e1, int next_e2,
		int **path, int *path_len)
{
	*path = NULL;
	*path_len = 0;
	int e1 = emap1->gl_e;
	int e2 = emap2->gl_e;
	int lc_e1 = emap1->lc_e;
	int lc_e2 = emap2->lc_e;
	struct graph_info_t ginfo;
	graph_info_init(lg, &ginfo, lc_e1, lc_e2);
	struct asm_edge_t edge_pre_e1;
	struct asm_edge_t edge_next_e2;
	if (pre_e1 == -1)
		edge_pre_e1.source = -1;
	else
		edge_pre_e1 = g->edges[pre_e1];
	if (next_e2 == -1)
		edge_next_e2.source = -1;
	else
		edge_next_e2 = g->edges[next_e2];
	unrelated_filter(g->edges[e1], g->edges[e2], edge_pre_e1,
			edge_next_e2, lg, &ginfo);
	struct path_info_t pinfo;
	path_info_init(&pinfo);
	get_all_paths(lg, &ginfo, &pinfo);
	if (pinfo.n_paths == 0)
		goto end_function;
	__VERBOSE_LOG("PATH", "Found %d paths, finding the best one\n",
			pinfo.n_paths);
	float *scores;
	get_path_scores(opt, g, lg, &pinfo, e1, e2, &scores);
	float best_score = 0;
	int best_path = 0;
	for (int i = 0; i < pinfo.n_paths; ++i){
		if (scores[i] > best_score){
			best_path = i;
			best_score = scores[i];
		}
	}
	__VERBOSE_LOG("", "Found best path id: %d, scores: %.3f\n",
			best_path, best_score);
	for (int i = 0; i < pinfo.path_lens[best_path]; ++i)
		__VERBOSE("%d ", pinfo.paths[best_path][i]);
	__VERBOSE("\n");
	__VERBOSE_LOG("", "-----------------------------------------\n");
	*path_len = pinfo.path_lens[best_path];
	*path = (int *) calloc(*path_len, sizeof(int));
	memcpy(*path, pinfo.paths[best_path], sizeof(int) * *path_len);
	free(scores);
end_function:
	graph_info_destroy(&ginfo);
	path_info_destroy(&pinfo);
}

void get_path_scores(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, struct path_info_t *pinfo,
		int e1, int e2, float **scores)
{
	int *seq_lens = (int *) calloc(pinfo->n_paths, sizeof(int));
	char cand_path[1024];
	sprintf(cand_path, "%s/%d_%d_all.fasta", opt->out_dir, e1, e2);
	FILE *f = fopen(cand_path, "w");
	for (int i = 0; i < pinfo->n_paths; ++i){
		fprintf(f, ">%d\n", i);
		char *seq;
		join_bridge_center_by_path(lg, pinfo->paths[i],
				pinfo->path_lens[i], &seq);
		fprintf(f, "%s\n", seq);
		seq_lens[i] = strlen(seq);
		free(seq);
	}
	fclose(f);
	khash_t(contig_count) *ctg_cnt = kh_init(contig_count);
	for (int i = 0; i < pinfo->n_paths; ++i){
		int ret;
		khiter_t it = kh_put(contig_count, ctg_cnt, i, &ret);
		kh_val(ctg_cnt, it) = 0;
	}

	struct read_path_t read_path;
	read_path.R1_path = (char *) calloc(1024, sizeof(char));
	sprintf(read_path.R1_path, "%s/local_assembly_%d_%d/R1.sub.fq",opt->out_dir,
			g->edges[e1].rc_id, e2);
	read_path.R2_path = (char *) calloc(1024, sizeof(char));
	sprintf(read_path.R2_path, "%s/local_assembly_%d_%d/R2.sub.fq",opt->out_dir,
			g->edges[e1].rc_id, e2);
	count_readpair_path(opt->n_threads, &read_path, cand_path, ctg_cnt);
	*scores = (float *) calloc(pinfo->n_paths, sizeof(float));
	for (khiter_t it = kh_begin(ctg_cnt); it != kh_end(ctg_cnt); ++it){
		if (!kh_exist(ctg_cnt, it))
			continue;
		int key = kh_key(ctg_cnt, it);
		int val = kh_val(ctg_cnt, it);
		(*scores)[key] = 1.0f * val / seq_lens[key];
	}
}

void join_seq(char **dest, char *source)
{
	int old_len = strlen(*dest);
	int new_len = old_len + strlen(source);
	char *old_dest = *dest;
	*dest = (char *) realloc(*dest, (new_len + 1) * sizeof(char));
	strcpy(*dest + old_len, source);
}

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
	__VERBOSE("Joining from %d to %d\n", lc_e1, lc_e2);
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

void get_contig_from_scaffold_path(struct opt_proc_t *opt, struct asm_graph_t *g,
		int *path, int path_len, char **contig)
{
	*contig = (char *) calloc(1, sizeof(char));
	{
		char *tmp;
		decode_seq(&tmp, g->edges[path[0]].seq, g->edges[path[0]].seq_len);
		join_seq(contig, tmp);
		free(tmp);
	}
	int bridge_types[N_BRIDGE_TYPE] = {};
	char gap_path[1024];
	sprintf(gap_path, "%s/gap.txt", opt->out_dir);
	FILE *f = fopen(gap_path, "w");
	for (int i = 1; i < path_len; ++i){
		int pre_u = i == 1 ? -1 : path[i - 2];
		int u = path[i - 1];
		int v = path[i];
		int next_v = i + 1 == path_len ? -1 : path[i + 1];
		struct asm_graph_t lg = get_local_assembly(opt, g,
				g->edges[u].rc_id, v);
		__VERBOSE("\n+------------------------------------------------------------------------------+\n");
		__VERBOSE_LOG("INFO", "Processing %d on %d bridges\n", i, path_len - 1);
		__VERBOSE_LOG("PATH", "Building bridge from %d to %d\n", u, v);
		char *seq;
		int leng;
		int res = get_bridge(opt, g, &lg, u, v, pre_u, next_v,
				&seq, &leng);
		join_seq(contig, seq + g->edges[path[i - 1]].seq_len);
		int closed_gap = max(1, leng - g->edges[path[i - 1]].seq_len
			- g->edges[path[i]].seq_len);

		/*if (res == LOCAL_NOT_FOUND){
			char *dump_N;
			get_dump_N(&dump_N);
			join_seq(contig, dump_N);
			free(dump_N);
			char *tmp;
			decode_seq(&tmp, 
			join_seq(contig, 
		} else if (res == TRIVIAL_BRIDGE){
			__VERBOSE_LOG("", "Trivial bridge found\n");
			join_seq(contig, seq + g->edges[path[i - 1]].seq_len);
			closed_gap = max(1, leng - g->edges[path[i - 1]].seq_len
				- g->edges[path[i]].seq_len);
		} else if (res == MULTIPLE_PATH){
			__VERBOSE_LOG("", "Multiple paths found\n");
			closed_gap = max(1, leng - g->edges[path[i - 1]].seq_len
				- g->edges[path[i]].seq_len);
			join_seq(contig, seq + g->edges[path[i - 1]].seq_len);
		} else {
			__VERBOSE_LOG("", "No path found\n");
			char *dump_N;
			get_dump_N(&dump_N);
			join_seq(contig, dump_N);

			char *tmp;
			decode_seq(&tmp, g->edges[path[i]].seq,
					g->edges[path[i]].seq_len);
			join_seq(contig, tmp);

			free(tmp);
			free(dump_N);
		}*/
		++bridge_types[res];
		__VERBOSE_LOG("GAP", "Closed gap: %d\n", closed_gap);
		fprintf(f, "Gap from %d to %d: %d\n", path[i],
				path[i - 1], closed_gap);
		fprintf(f, "%d %d %d\n", g->edges[path[i - 1]].seq_len,
				g->edges[path[i]].seq_len, leng);
		free(seq);
	}
	fclose(f);
	__VERBOSE_LOG("INFO", "Path summary:\n");
	__VERBOSE_LOG("", "Number of trivial bridges: %d\n",
			bridge_types[BRIDGE_TRIVIAL_BRIDGE]);
	__VERBOSE_LOG("", "Number of cases where local edges not found: %d\n",
			bridge_types[BRIDGE_LOCAL_NOT_FOUND]);
	__VERBOSE_LOG("", "Number of multiple paths: %d\n",
			bridge_types[BRIDGE_MULTIPLE_PATH]);
	__VERBOSE_LOG("", "Number of disconnected region: %d\n",
			bridge_types[BRIDGE_PATH_NOT_FOUND]);
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

void join_bridge_no_path(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_edge_t lc_e1, struct asm_edge_t lc_e2,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2,
		char **res_seq)
{
	*res_seq = (char *) calloc(1, sizeof(char));
	char *first, *second;
	sync_global_local_edge(e1, lc_e1, emap1->gpos, emap1->lpos,
			SYNC_KEEP_GLOBAL, &first);
	sync_global_local_edge(e2, lc_e2, emap2->gpos, emap1->lpos,
			SYNC_KEEP_LOCAL, &second);
	char *dump_N;
	get_dump_N(&dump_N);
	join_seq(res_seq, first);
	join_seq(res_seq, dump_N);
	join_seq(res_seq, second);
	free(first);
	free(second);
	free(dump_N);
}

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
