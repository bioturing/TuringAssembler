#include "build_bridge.h"
#include "helper.h"

int find_path(struct asm_graph_t lg, gint_t u, int start_edge, int end_edge,
		khash_t(gint_t_int) *visited, int *is_disable, int depth,
		int **path, int *path_leng)
{
	if (u == end_edge){
		*path_leng = depth + 1;
		*path = (int *) calloc(*path_leng, sizeof(int));
		goto path_found;
	}
	int tg = lg.edges[u].target;
	for (int i = 0; i < lg.nodes[tg].deg; ++i){
		gint_t v = lg.nodes[tg].adj[i];
		if (v == start_edge)
			continue;
		if (is_disable[v])
			continue;
		if (lg.edges[v].rc_id == end_edge)
			continue;
		gint_t edge_code = get_edge_code(u, v);
		if (kh_get(gint_t_int, visited, edge_code) != kh_end(visited))
			continue;
		int ret;
		kh_put(gint_t_int, visited, edge_code, &ret);
		if (find_path(lg, v, start_edge, end_edge, visited,
				is_disable, depth + 1, path, path_leng))
			goto path_found;
	}
	return 0;
path_found:
	(*path)[depth] = u;
	return 1;
}

void find_all_paths(struct asm_graph_t g, gint_t u, int start_edge, int end_edge,
		khash_t(gint_t_int) *visited, int *is_disable, int depth,
		int *path, FILE *record, int *count_paths)
{
	path[depth] = u;
	if (u == end_edge){
		/*uint32_t *seq;
		uint32_t leng;
		combine_edges(g, &seq, &leng, g.edges[start_edge].seq_len,
				g.edges[start_edge].seq_len - 1,
				g.edges[end_edge].seq_len, 0, path, depth + 1);
		fprintf(record, ">%d_%d_path_%d\n", start_edge, end_edge,
				*count_paths);
		__VERBOSE("%d %d\n", start_edge, end_edge);
		__VERBOSE("%d %d %d\n", g.edges[start_edge].seq_len,
				g.edges[end_edge].seq_len, leng);
		for (int i = 0; i < leng; ++i)
			fprintf(record, "%c", int_to_base(__binseq_get(seq, i)));
		fprintf(record, "\n");*/
		++(*count_paths);
		return;
	}
	int tg = g.edges[u].target;
	for (int i = 0; i < g.nodes[tg].deg; ++i){
		gint_t v = g.nodes[tg].adj[i];
		if (v == start_edge)
			continue;
		if (is_disable[v])
			continue;
		if (g.edges[v].rc_id == end_edge)
			continue;
		gint_t edge_code = get_edge_code(u, v);
		khiter_t it = kh_get(gint_t_int, visited, edge_code);
		if (it == kh_end(visited)){
			int ret;
			it = kh_put(gint_t_int, visited, edge_code, &ret);
			kh_val(visited, it) = 0;
		}
		if (kh_val(visited, it) == 1)
			continue;
		++kh_val(visited, it);
		find_all_paths(g, v, start_edge, end_edge, visited, is_disable,
			depth + 1, path, record, count_paths);
		it = kh_get(gint_t_int, visited, edge_code);
		--kh_val(visited, it);
	}
}

void cov_filter(struct asm_graph_t g, int start_edge, int end_edge,
		int *is_disable)
{
	__VERBOSE_LOG("COV FILTER", "+----------------------------------+\n");
	__VERBOSE_LOG("", "Before filter: %ld edges\n", g.n_e);
	int thresh = (int) (MIN_DEPTH_RATIO * max(get_cov(g, start_edge),
						get_cov(g, end_edge)));
	int disabled = 0;
	for (int i = 0; i < g.n_e; ++i){
		is_disable[i] = get_cov(g, i) < thresh;
		disabled += is_disable[i];
	}
	__VERBOSE_LOG("", "After filter: %ld edges\n", g.n_e - disabled);
}

int dfs(struct asm_graph_t g, int u, int start_edge, int end_edge, int *visited,
		int *is_disable)
{
	if (u == end_edge)
		return 1;
	if (visited[u]++)
		return 0;
	int tg = g.edges[u].target;
	//__VERBOSE_LOG("", "%d ", u);
	for (int i = 0; i < (int) g.nodes[tg].deg; ++i){
		int v = g.nodes[tg].adj[i];
		if (is_disable[v])
			continue;
		if (v == g.edges[start_edge].rc_id ||
			v == g.edges[end_edge].rc_id)
			continue;
		if (dfs(g, v, start_edge, end_edge, visited, is_disable))
			return 1;
	}
	return 0;
}

int reachable(struct asm_graph_t g, int start_edge, int end_edge,
		int *is_disable)
{
	int *visited = (int *) calloc(g.n_e, sizeof(int));
	int res = dfs(g, start_edge, start_edge, end_edge, visited, is_disable);
	free(visited);
	return res;
}

void bfs(struct asm_graph_t g, int start_edge, int *is_disable, int **path_leng)
{
	int *queue = (int *) calloc(g.n_e, sizeof(int));
	*path_leng = (int *) calloc(g.n_e, sizeof(int));
	for (int i = 0; i < (int) g.n_e; ++i)
		(*path_leng)[i] = -1;
	(*path_leng)[start_edge] = 0;
	int front = 0, back = 1;
	queue[0] = start_edge;
	while (front < back){
		int u = queue[front++];
		int tg = g.edges[u].target;
		for (int i = 0; i < (int) g.nodes[tg].deg; ++i){
			int v = g.nodes[tg].adj[i];
			if (is_disable[v])
				continue;
			if ((*path_leng)[v] != -1)
				continue;
			float a = get_cov(g, u);
			float b = get_cov(g, v);
			if (a / b < MIN_RELATIVE_DEPTH_RATIO ||
				b / a < MIN_RELATIVE_DEPTH_RATIO)
				continue;
			(*path_leng)[v] = (*path_leng)[u] + 1;
			queue[back++] = v;
		}
	}
	free(queue);
}

void connection_filter(struct asm_graph_t g, int start_edge, int end_edge,
			int *is_disable)
{
	__VERBOSE_LOG("CONNECTION FILTER", "+----------------------------------+\n");
	int disabled = 0;
	for (int i = 0; i < g.n_e; ++i)
		disabled += is_disable[i];
	__VERBOSE_LOG("", "Before filter: %ld edges\n", g.n_e - disabled);
	int *forward_leng, *backward_leng;
	bfs(g, start_edge, is_disable, &forward_leng);
	bfs(g, g.edges[end_edge].rc_id, is_disable, &backward_leng);
	for (int i = 0; i < g.n_e; ++i){
		if (is_disable[i])
			continue;
		int l1 = forward_leng[i];
		int l2 = backward_leng[g.edges[i].rc_id];
		if (l1 == -1 || l2 == -1 || l1 + l2 > MIN_PATH_LENGTH){
			is_disable[i] = 1;
			++disabled;
		}
	}
	__VERBOSE_LOG("", "After filter: %ld edges\n", g.n_e - disabled);
}


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

void filter_edges(struct asm_graph_t g, int start_edge, int end_edge,
		int **is_disable)
{
	*is_disable = (int *) calloc(g.n_e, sizeof(int));
	cov_filter(g, start_edge, end_edge, *is_disable);
	connection_filter(g, start_edge, end_edge, *is_disable);
}

void get_path(struct asm_graph_t lg, int start_edge, int end_edge,
		int **path, int *path_leng)
{
	int *is_disable;
	filter_edges(lg, start_edge, end_edge, &is_disable);
	khash_t(gint_t_int) *visited = kh_init(gint_t_int);
	*path_leng = 0;
	*path = NULL;
	int found = find_path(lg, start_edge, start_edge, end_edge, visited,
			is_disable, 0, path, path_leng);
	kh_destroy(gint_t_int, visited);
	free(is_disable);
}

void print_path(int *path, int path_leng)
{
	for (int i = 0; i < path_leng; ++i)
		__VERBOSE_LOG("", "%d ", path[i]);
	__VERBOSE_LOG("", "\n");
}

gint_t get_edge_code(gint_t e1, gint_t e2)
{
	return (e1 << 32) | e2;
}

int check_simple_path(struct asm_graph_t lg, int start_edge, int end_edge)
{
	khash_t(gint_t_int) *visited = kh_init(gint_t_int);
	int *path = NULL;
	int path_leng;
	int *new_path = NULL;
	int new_path_leng;
	int *is_disable;
	get_path(lg, start_edge, end_edge, &path, &path_leng);
	__VERBOSE_LOG("PATH", "++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	__VERBOSE_LOG("", "Found first path:\n");
	print_path(path, path_leng);
	__VERBOSE_LOG("", "\n");

	int res = 1;
	int ban_edge1 = -1;
	int ban_edge2 = -1;
	filter_edges(lg, start_edge, end_edge, &is_disable);
	for (int i = 1; i < path_leng; ++i){
		kh_destroy(gint_t_int, visited);
		visited = kh_init(gint_t_int);
		ban_edge1 = path[i - 1];
		ban_edge2 = path[i];
		gint_t edge_code = get_edge_code(ban_edge1, ban_edge2);
		int ret;
		kh_put(gint_t_int, visited, edge_code, &ret);
		int tmp = find_path(lg, start_edge, start_edge, end_edge,
				visited, is_disable, 0, &new_path,
				&new_path_leng);
		if (tmp){
			res = 0;
			__VERBOSE_LOG("", "Found second path:\n");
			print_path(new_path, new_path_leng);
			__VERBOSE_LOG("", "\n");
			goto end_check_simple_path;
		}
	}
end_check_simple_path:
	kh_destroy(gint_t_int, visited);
	if (path != NULL)
		free(path);
	if (new_path != NULL)
		free(new_path);
	free(is_disable);
	return res;
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

void get_midair_bridge(struct asm_graph_t lg, int start_edge, int end_edge,
		int *midair_edge)
{
	int *is_disable;
	filter_edges(lg, start_edge, end_edge, &is_disable);
	int *bad = (int *) calloc(lg.n_e, sizeof(int));
	for (int i = 0; i < lg.n_e; ++i)
		bad[i] = is_disable[i];
	bad[start_edge] = bad[end_edge] = bad[lg.edges[start_edge].rc_id]
		= bad[lg.edges[end_edge].rc_id] = 1;
	for (int i = 0; i < lg.n_e; ++i){
		if (is_disable[i] || is_disable[lg.edges[i].rc_id])
			continue;
		is_disable[lg.edges[i].rc_id] = 1;
		if (reachable(lg, start_edge, i, is_disable) &&
			reachable(lg, i, end_edge, is_disable))
			bad[lg.edges[i].rc_id] = 1;
		is_disable[lg.edges[i].rc_id] = 0;
	}
	int res = -1;
	for (int i = 0; i < lg.n_e; ++i){
		if (bad[i])
			continue;
		if (res == -1 || lg.edges[res].seq_len < lg.edges[i].seq_len)
			res = i;
	}
	free(bad);
	free(is_disable);
	*midair_edge = res;
}

void get_local_edge_head(struct asm_graph_t g, struct asm_graph_t lg,
		struct asm_edge_t e, int *edge_id, struct subseq_pos_t *gpos,
		struct subseq_pos_t *lpos)
{
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
	} else if (lg->edges[lc_e1].rc_id == lc_e2){
		get_local_edge_tail(*g, *lg, g->edges[g->edges[e2].rc_id],
				&lc_e2, &gpos2, &lpos2);
		bridge_type = MIS_SCAFFOLD;
		join_trivial_bridge(g->edges[e1], g->edges[e2], *lg, lc_e1,
				gpos1, lpos1, gpos2, lpos2, &bridge_seq);
		goto path_found;
	} else {
		int *path;
		int path_len;
		get_path(*lg, lc_e1, lc_e2, &path, &path_len);
		/*if (path != NULL){
			get_path_max_cov(*lg, lc_e1, lc_e2, &path, &path_len);
		} else {
			goto path_not_found;
		}*/
		if (path != NULL && check_simple_path(*lg, lc_e1, lc_e2)){
			bridge_type = TRIVIAL_CASE;
			join_bridge_by_path(g->edges[e1], g->edges[e2], *lg,
					path, path_len, gpos1, lpos1, gpos2,
					lpos2, &bridge_seq);
			goto path_found;
		} else {
			bridge_type = COMPLEX_PATH;
			__VERBOSE_LOG("", "Complex path, searching for middle edge\n");
			goto path_not_found;
			int mid_edge;
			get_midair_bridge(*lg, lc_e1, lc_e2, &mid_edge);
			__VERBOSE_LOG("", "Middle edge found: %d\n", mid_edge);
			if (mid_edge != -1){
				decode_seq(&bridge_seq, lg->edges[mid_edge].seq,
						lg->edges[mid_edge].seq_len);
				goto path_found;
			}
			goto path_not_found;
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
	__VERBOSE("Joinging form %d to %d\n", lc_e1, lc_e2);
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
		struct asm_graph_t lg = test_local_assembly(opt, g,
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
		} else if (res == MIS_SCAFFOLD){
			__VERBOSE_LOG("", "Mis scaffold\n");
			char *tmp;
			decode_seq(&tmp, seq, leng);
			join_seq(contig, tmp + g->edges[path[i - 1]].seq_len);
			free(tmp);
			path[i] = g->edges[path[i]].rc_id;
		} else if (res == COMPLEX_PATH){
			char *dump_N;
			get_dump_N(&dump_N);
			join_seq(contig, dump_N);

			char *tmp;
			decode_seq(&tmp, seq, leng);
			join_seq(contig, tmp);
			free(tmp);

			join_seq(contig, dump_N);

			decode_seq(&tmp, g->edges[path[i]].seq,
					g->edges[path[i]].seq_len);
			join_seq(contig, tmp);

			free(tmp);
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

void find_path_max_cov(struct asm_graph_t lg, gint_t u, int start_edge, int end_edge,
		khash_t(gint_t_int) *visited, int *is_disable, int depth,
		int *tmp_path, int cov_sum, int *max_cov, int **path,
		int *path_leng)
{
	tmp_path[depth] = u;
	cov_sum += get_cov(lg, u);
	if (u == end_edge){
		if (*max_cov < cov_sum){
			*path_leng = depth + 1;
			*path = (int *) realloc(*path, sizeof(int) * *path_leng);
			for (int i = 0; i < *path_leng; ++i)
				(*path)[i] = tmp_path[i];
			*max_cov = cov_sum;
		}
		return;
	}
	int tg = lg.edges[u].target;
	for (int i = 0; i < lg.nodes[tg].deg; ++i){
		gint_t v = lg.nodes[tg].adj[i];
		if (v == start_edge)
			continue;
		if (is_disable[v])
			continue;
		if (lg.edges[v].rc_id == end_edge)
			continue;
		gint_t edge_code = get_edge_code(u, v);
		khiter_t it = kh_get(gint_t_int, visited, edge_code);
		if (it == kh_end(visited)){
			int ret;
			khiter_t it = kh_put(gint_t_int, visited, edge_code,
					&ret);
			kh_val(visited, it) = 0;
		}
		if (kh_val(visited, it) == 1)
			continue;
		++kh_val(visited, it);
		find_path_max_cov(lg, v, start_edge, end_edge, visited,
				is_disable, depth + 1, tmp_path, cov_sum,
				max_cov, path, path_leng);
		it = kh_get(gint_t_int, visited, edge_code);
		--kh_val(visited, it);
	}
}

void get_path_max_cov(struct asm_graph_t lg, int start_edge, int end_edge,
		int **path, int *path_len)
{
	int *tmp_path = (int *) calloc(lg.n_e, sizeof(int));
	khash_t(gint_t_int) *visited = kh_init(gint_t_int);
	int *is_disable;
	int max_cov = 0;
	filter_edges(lg, start_edge, end_edge, &is_disable);
	find_path_max_cov(lg, start_edge, start_edge, end_edge, visited,
			is_disable, 0, tmp_path, 0, &max_cov, path, path_len);
	free(is_disable);
	free(tmp_path);
	kh_destroy(gint_t_int, visited);
}
