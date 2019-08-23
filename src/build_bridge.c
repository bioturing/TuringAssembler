#include "build_bridge.h"
#include "helper.h"
#include "graph_search.h"

void print_graph(struct asm_graph_t g, int *is_disable, char *path)
{
	FILE *f = fopen(path, "w");
	int *mark = (int *) calloc(g.n_e, sizeof(int));
	for (int u = 0; u < (int) g.n_e; ++u){
		if (is_disable[u])
			continue;
		if (mark[u] || mark[g.edges[u].rc_id])
			continue;
		mark[u] = mark[g.edges[u].rc_id] = 1;
		fprintf(f, "S\t%d_%ld_%.3f\t", u, g.edges[u].rc_id,
				get_cov(g, u));
		for (int i = 0; i < (int)g.edges[u].seq_len; ++i)
			fprintf(f, "%c", int_to_base(__binseq_get(g.edges[u].seq, i)));
		fprintf(f, "\tKC:i:%ld\n", g.edges[u].count);
	}
	for (int u = 0; u < (int) g.n_e; ++u)
		mark[u] = 0;
	for (int u = 0; u < (int) g.n_e; ++u){
		if (is_disable[u])
			continue;
		if (mark[u] || mark[g.edges[u].rc_id])
			continue;
		mark[u] = mark[g.edges[u].rc_id] = 1;
		int tg = g.edges[u].target;
		for (int i = 0; i < (int) g.nodes[tg].deg; ++i){
			int v = g.nodes[tg].adj[i];
			if (is_disable[v])
				continue;
			fprintf(f, "L\t%d_%ld_%.3f\t+\t%d_%ld_%.3f\t+\t%dM\n", u,
					g.edges[u].rc_id, get_cov(g, u),
					v, g.edges[v].rc_id, get_cov(g, v),
					g.ksize);
		}
	}
	fclose(f);
}

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
		struct asm_edge_t e, int *edge_id, struct subseq_pos_t *gpos,
		struct subseq_pos_t *lpos)
{
	struct map_contig_t mct;
	init_map_contig(&mct, g.edges[e.rc_id], lg);
	*edge_id = find_match(&mct);
	__VERBOSE("%d\n", *edge_id);
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
		struct asm_edge_t e, int *edge_id, struct subseq_pos_t *gpos,
		struct subseq_pos_t *lpos)
{
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

int get_bridge(struct asm_graph_t *g, struct asm_graph_t *lg, int e1, int e2,
		uint32_t **ret_seq, uint32_t *seq_len)
{
	__VERBOSE("Matching edges...\n");
	int lc_e1;
	struct subseq_pos_t gpos1, lpos1;
	get_local_edge_head(*g, *lg, g->edges[e1], &lc_e1, &gpos1, &lpos1);

	__VERBOSE_LOG("", "Local edge 1: %d\n", lc_e1);
	__VERBOSE_LOG("", "Global edge starts from: %d, ends at: %d\n", gpos1.start,
			gpos1.end);
	__VERBOSE_LOG("", "Local edge starts from: %d, ends at: %d\n", lpos1.start,
			lpos1.end);

	int lc_e2;
	struct subseq_pos_t gpos2, lpos2;
	get_local_edge_tail(*g, *lg, g->edges[e2], &lc_e2, &gpos2, &lpos2);
	__VERBOSE_LOG("", "Local edge 2: %d\n", lc_e2);
	__VERBOSE_LOG("", "Global edge starts from: %d, ends at: %d\n", gpos2.start,
			gpos2.end);
	__VERBOSE_LOG("", "Local edge starts from: %d, ends at: %d\n", lpos2.start,
			lpos2.end);

	int res = try_bridging(g, lg, e1, e2, ret_seq, seq_len, lc_e1, lc_e2,
			gpos1, lpos1, gpos2, lpos2);
	/*if (res == PATH_NOT_FOUND){
		__VERBOSE_LOG("", "Checkiing if edge is reversed\n");
		get_local_edge_head(*g, *lg, g->edges[e2], &lc_e2, &gpos2,
				&lpos2);
		__VERBOSE_LOG("", "Local edge 2 reversed: %d\n", lc_e2);
		__VERBOSE_LOG("", "Global edge starts from: %d, ends at: %d\n",
				gpos2.start, gpos2.end);
		__VERBOSE_LOG("", "Local edge starts from: %d, ends at: %d\n",
				lpos2.start, lpos2.end);
		int res = try_bridging(g, lg, e1, e2, ret_seq, seq_len, lc_e1,
				lc_e2, gpos1, lpos1, gpos2, lpos2);
		if (res == TRIVIAL_CASE)
			res = MIS_SCAFFOLD;
	}*/
	return res;
}

int try_bridging(struct asm_graph_t *g, struct asm_graph_t *lg, int e1, int e2,
		uint32_t **ret_seq, uint32_t *seq_len, int lc_e1, int lc_e2,
		struct subseq_pos_t gpos1, struct subseq_pos_t lpos1,
		struct subseq_pos_t gpos2, struct subseq_pos_t lpos2)
{
	int bridge_type;
	char *bridge_seq;
	if (lc_e1 == -1 || lc_e2 == -1){
		goto path_not_found;
	} else if (lc_e1 == lc_e2){
		bridge_type = TRIVIAL_CASE;
		join_trivial_bridge(g->edges[e1], g->edges[e2], *lg, lc_e1,
				gpos1, lpos1, gpos2, lpos2, &bridge_seq);
		goto path_found;
	} else {
		int *path;
		int path_len;
		int mid_edge;
		int path_type = get_path(lg, lc_e1, lc_e2, &mid_edge,
				&path, &path_len);
		if (path_type == SIMPLE_PATH){
			bridge_type = TRIVIAL_CASE;
			join_bridge_by_path(g->edges[e1], g->edges[e2], *lg,
					path, path_len, gpos1, lpos1, gpos2,
					lpos2, &bridge_seq);
			goto path_found;
		} else {
			bridge_type = NON_TRIVIAL_CASE;
			goto path_not_found;
			/*if (mid_edge != -1){
				decode_seq(&bridge_seq, lg->edges[mid_edge].seq,
						lg->edges[mid_edge].seq_len);
				goto path_found;
			}
			goto path_not_found;*/
		}
	}
path_not_found:
	*ret_seq = NULL;
	*seq_len = 0;
	bridge_type = PATH_NOT_FOUND;
	goto end_function;
path_found:
	*seq_len = strlen(bridge_seq);
	encode_seq(ret_seq, bridge_seq);
	free(bridge_seq);
end_function:
	return bridge_type;
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
		struct asm_graph_t lg, int local_edge,
		struct subseq_pos_t gpos1, struct subseq_pos_t lpos1,
		struct subseq_pos_t gpos2, struct subseq_pos_t lpos2,
		char **res_seq)
{
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
		struct subseq_pos_t gpos1, struct subseq_pos_t lpos1,
		struct subseq_pos_t gpos2, struct subseq_pos_t lpos2,
		char **res_seq)
{
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
	char *tmp;
	decode_seq(&tmp, g->edges[path[0]].seq, g->edges[path[0]].seq_len);
	join_seq(contig, tmp);
	free(tmp);
	for (int i = 1; i < path_len; ++i){
		int u = path[i - 1];
		int v = path[i];
		struct asm_graph_t lg = get_local_assembly(opt, g,
				g->edges[u].rc_id, v);
		__VERBOSE("\n+------------------------------------------------------------------------------+\n");
		__VERBOSE_LOG("INFO", "Processing %d on %d bridges\n", i, path_len - 1);
		__VERBOSE_LOG("PATH", "Building bridge from %d to %d\n", u, v);
		uint32_t *seq;
		uint32_t leng;
		int res = get_bridge(g, &lg, u, v, &seq, &leng);

		if (res == TRIVIAL_CASE){
			__VERBOSE_LOG("", "Trivial path found\n");
			char *tmp;
			decode_seq(&tmp, seq, leng);
			join_seq(contig, tmp + g->edges[path[i - 1]].seq_len);
			free(tmp);
		} else if (res == NON_TRIVIAL_CASE && seq != NULL){
			__VERBOSE_LOG("", "Middle edge found\n");
			char *dump_N;
			get_dump_N(&dump_N);
			join_seq(contig, dump_N);

			char *tmp;
			decode_seq(&tmp, seq, leng);
			join_seq(contig, tmp);
			free(tmp);

			join_seq(contig, dump_N);

			free(dump_N);
		} else {
			char *dump_N;
			get_dump_N(&dump_N);
			join_seq(contig, dump_N);

			char *tmp;
			decode_seq(&tmp, g->edges[path[i]].seq,
					g->edges[path[i]].seq_len);
			join_seq(contig, tmp);

			free(tmp);
			free(dump_N);
		}
	}
}

