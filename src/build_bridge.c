#include "build_bridge.h"

int min(int a, int b)
{
	return a < b ? a : b;
}

int max(int a, int b)
{
	return a > b ? a : b;
}

int match_head(struct asm_edge_t P, struct asm_edge_t T)
{
	int Tleng = T.seq_len;
	int Pleng = P.seq_len;
	for (int i = Tleng - 1; i >= MIN_MATCH_LENG; --i){
		int p = 0;
		int unmatch = 0;
		int overlap_leng = min(i + 1, Pleng);
		int flag = 1;
		while (p < overlap_leng){
			int found = 0;
			for (int j = -MIN_RADIUS; j <= MIN_RADIUS; ++j){
				int new_p = i - p + j;
				if (new_p < 0 || new_p >= T.seq_len)
					continue;
				if(__binseq_get(T.seq, new_p) ==
					__binseq_get(P.seq, Pleng - p - 1)){
					found = 1;
					break;
				}
			}
			if (!found)
				++unmatch;
			if ((float) unmatch / overlap_leng > MIN_UNMATCHED_RATIO){
				flag = 0;
				break;
			}
			++p;
		}
		if (flag)
			return i;
	}
	return -1;
}

int match_tail(struct asm_edge_t P, struct asm_edge_t T)
{
	int Tleng = T.seq_len;
	int Pleng = P.seq_len;
	for (int i = 0; i < Tleng - MIN_MATCH_LENG; ++i){
		int p = 0;
		int unmatch = 0;
		int overlap_leng = min(Tleng - i, Pleng);
		int flag = 1;
		while (p < overlap_leng){
			int found = 0;
			for (int j = -MIN_RADIUS; j <= MIN_RADIUS; ++j){
				int new_p = i + p + j;
				if (new_p < 0 || new_p >= T.seq_len)
					continue;
				if (__binseq_get(T.seq, new_p) ==
					__binseq_get(P.seq, p)){
					found = 1;
					break;
				}
			}
			if (!found)
				++unmatch;
			if ((float) unmatch / overlap_leng > MIN_UNMATCHED_RATIO){
				flag = 0;
				break;
			}
			++p;
		}
		if (flag)
			return i;
	}
	return -1;
}

void get_local_edge_id_head(struct asm_graph_t lg, struct asm_edge_t e,
				int *edge_id, int *pos)
{
	*pos = -1;
	*edge_id = -1;
	int max_match = 0;
	for (int i = 0; i < lg.n_e; ++i){
		/*if (lg.edges[i].seq_len < MIN_EDGE_LENGTH)
			continue;*/
		/*if ((float) lg.edges[i].seq_len / e.seq_len < MIN_EDGE_LENGTH_RATIO)
			continue;*/
		int tmp_pos = match_head(e, lg.edges[i]);
		int match_leng = min((int) e.seq_len, tmp_pos + 1);
		/*if ((float) match_leng / e.seq_len < MIN_EDGE_LENGTH_RATIO)
			continue;*/
		//if (max_match < match_leng){
		if (i == 24){
			__VERBOSE("bla %d\n", lg.edges[i].seq_len);
			__VERBOSE("%d\n", tmp_pos);
		}
		if (*pos < tmp_pos){
			max_match = match_leng;
			*pos = tmp_pos;
			*edge_id = i;
		}
	}
}

void get_local_edge_id_tail(struct asm_graph_t lg, struct asm_edge_t e,
				int *edge_id, int *pos)
{
	*pos = -1;
	*edge_id = -1;
	int max_match = 0;
	for (int i = 0; i < lg.n_e; ++i){
		/*if (lg.edges[i].seq_len < MIN_EDGE_LENGTH)
			continue;*/
		/*if ((float) lg.edges[i].seq_len / e.seq_len < MIN_EDGE_LENGTH_RATIO)
			continue;*/
		int tmp_pos = match_tail(e, lg.edges[i]);
		int match_leng = min((int) e.seq_len, (int) lg.edges[i].seq_len - tmp_pos);
		/*if ((float) match_leng / e.seq_len < MIN_EDGE_LENGTH_RATIO)
			continue;*/
		//if (max_match < match_leng){
		if (tmp_pos != -1 && (*pos == -1 || *pos > tmp_pos)){
			max_match = match_leng;
			*pos = tmp_pos;
			*edge_id = i;
		}
	}
}

/*int find_path(struct asm_graph_t lg, int u, int end_edge, int *visited,
		int *trace, int ban_edge1, int ban_edge2)
{
	if (u == end_edge)
		return 1;
	visited[u] = 1;
	visited[lg.edges[u].rc_id] = 1;
	int tg = lg.edges[u].target;
	for (int i = 0; i < lg.nodes[tg].deg; ++i){
		int v = lg.nodes[tg].adj[i];
		if (u == ban_edge1 && v == ban_edge2)
			continue;
		if (lg.edges[v].rc_id == end_edge)
			continue;
		if (visited[v])
			continue;
		trace[v] = u;
		if (find_path(lg, v, end_edge, visited, trace, ban_edge1,
					ban_edge2))
			return 1;
	}
	return 0;
}*/

float get_cov(struct asm_graph_t g, int edge)
{
	return (float) g.edges[edge].count / (g.edges[edge].seq_len - g.ksize);
}

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

/*int find_path(struct asm_graph_t lg, gint_t u, int start_edge, int end_edge,
		khash_t(gint_t_int) *visited, int depth, int **path,
		int *path_leng)
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
		if (lg.edges[v].rc_id == end_edge)
			continue;
		gint_t edge_code = get_edge_code(u, v);
		if (kh_get(gint_t_int, visited, edge_code) != kh_end(visited))
			continue;
		int ret;
		kh_put(gint_t_int, visited, edge_code, &ret);
		if (find_path(lg, v, start_edge, end_edge, visited, depth + 1,
					path, path_leng))
			goto path_found;
	}
	return 0;
path_found:
	(*path)[depth] = u;
	return 1;
}*/

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
	//__VERBOSE_LOG("", "PATH\n");
	int res = dfs(g, start_edge, start_edge, end_edge, visited, is_disable);
	//__VERBOSE_LOG("", "END PATH\n");
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
	/*__VERBOSE_LOG("CONNECTION FILTER", "+----------------------------------+\n");
	int disabled = 0;
	for (int i = 0; i < g.n_e; ++i)
		disabled += is_disable[i];
	__VERBOSE_LOG("", "Before filter: %ld edges\n", g.n_e - disabled);
	for (int i = 0; i < g.n_e; ++i){
		if (is_disable[i])
			continue;
		if (reachable(g, start_edge, i, is_disable) == 0 ||
			reachable(g, i, end_edge, is_disable) == 0)
			is_disable[i] = 1;
		disabled += is_disable[i];
	}
	__VERBOSE_LOG("", "After filter: %ld edges\n", g.n_e - disabled);*/
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

char int_to_base(int x)
{
	if (x == 0)
		return 'A';
	if (x == 1)
		return 'C';
	if (x == 2)
		return 'G';
	return 'T';
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

/*void get_path(struct asm_graph_t lg, int start_edge, int end_edge,
		int **path, int *path_leng)
{
	khash_t(gint_t_int) *visited = kh_init(gint_t_int);
	*path_leng = 0;
	*path = NULL;
	int found = find_path(lg, start_edge, start_edge, end_edge, visited, 0,
			path, path_leng);
	kh_destroy(gint_t_int, visited);
}*/

void filter_edges(struct asm_graph_t g, int start_edge, int end_edge,
		int **is_disable)
{
	*is_disable = (int *) calloc(g.n_e, sizeof(int));
	//cov_filter(g, start_edge, end_edge, *is_disable);
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

/*void get_path(struct asm_graph_t lg, int start_edge, int end_edge,
		int **path, int *path_leng)
{
	int *visited = (int *) calloc(lg.n_e, sizeof(int));
	int *trace = (int *) calloc(lg.n_e, sizeof(int));
	trace[start_edge] = -1;
	int found = find_path(lg, start_edge, end_edge, visited, trace, -1, -1);
	if (!found)
		goto path_not_found;

	*path_leng = 0;
	for (int u = end_edge; u != -1; u = trace[u])
		++(*path_leng);

	*path = (int *) calloc(*path_leng, sizeof(int));
	int p = *path_leng - 1;
	for (int u = end_edge; u != -1; u = trace[u]){
		(*path)[p] = u;
		--p;
	}
	goto path_found;
path_not_found:
	*path = NULL;
	*path_leng = 0;
path_found:
	free(visited);
	free(trace);
}*/

/*void print_path(struct asm_graph_t lg, int start_edge, int u, int *trace)
{
	if (u == trace[start_edge])
		return;
	print_path(lg, start_edge, trace[u], trace);
	__VERBOSE_LOG("", "%d ", u);
}*/

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

/*int check_simple_path(struct asm_graph_t lg, int start_edge, int end_edge)
{
	int *visited = (int *) calloc(lg.n_e, sizeof(int));
	int *trace = (int *) calloc(lg.n_e, sizeof(int));
	int *new_trace = (int *) calloc(lg.n_e, sizeof(int));

	trace[start_edge] = -1;
	find_path(lg, start_edge, end_edge, visited, trace, -1, -1);
	__VERBOSE("++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	__VERBOSE("Found first path:\n");
	print_path(lg, start_edge, end_edge, trace);
	__VERBOSE("\n");

	int res = 1;
	int ban_edge1 = -1;
	int ban_edge2 = -1;
	for (int u = end_edge; u != start_edge; u = trace[u]){
		for (int i = 0; i < lg.n_e; ++i)
			visited[i] = 0;
		ban_edge1 = trace[u];
		ban_edge2 = u;
		new_trace[start_edge] = -1;
		int tmp = find_path(lg, start_edge, end_edge, visited,
				new_trace, ban_edge1, ban_edge2);
		if (tmp){
			res = 0;
			__VERBOSE("Found second path:\n");
			print_path(lg, start_edge, end_edge, new_trace);
			__VERBOSE("\n");
			goto end_check_simple_path;
		}
	}
end_check_simple_path:
	free(visited);
	free(trace);
	free(new_trace);
	return res;
}*/

void combine_edges(struct asm_graph_t lg, uint32_t **seq, uint32_t *leng,
			int head_leng, int head_pos, int tail_leng,
			int tail_pos, int *path, int path_leng)
{
	int start_leng = lg.edges[path[0]].seq_len;
	int end_leng = lg.edges[path[path_leng - 1]].seq_len;
	*leng = head_leng + start_leng - head_pos - 1;
	for (int i = 1; i < path_leng - 1; ++i)
		*leng += lg.edges[path[i]].seq_len - lg.ksize;
	*leng += tail_pos + tail_leng - lg.ksize;
	*seq = (uint32_t *) calloc((*leng + 3) / 4, sizeof(uint32_t));

	int tmp = 0;
	for (int i = head_pos - head_leng + 1; i < start_leng; ++i){
		__binseq_set(*seq, tmp, __binseq_get(lg.edges[path[0]].seq, i));
		++tmp;
	}
	for (int i = 1; i < path_leng - 1; ++i){
		struct asm_edge_t e = lg.edges[path[i]];
		for (int j = lg.ksize; j < (int) e.seq_len; ++j){
			__binseq_set(*seq, tmp, __binseq_get(e.seq, j));
			++tmp;
		}
	}
	for (int i = lg.ksize; i < tail_pos + tail_leng; ++i){
		__binseq_set(*seq, tmp,
			__binseq_get(lg.edges[path[path_leng - 1]].seq, i));
		++tmp;
	}
}

/*void combine_edges(struct asm_graph_t lg, uint32_t **seq, uint32_t *leng,
			int head_leng, int head_pos, int tail_leng,
			int tail_pos, int *path, int path_leng)
{
	int start_leng = lg.edges[path[0]].seq_len;
	int end_leng = lg.edges[path[path_leng - 1]].seq_len;
	*leng = head_leng + start_leng - head_pos - 1;
	for (int i = 1; i < path_leng - 1; ++i)
		*leng += lg.edges[path[i]].seq_len;
	*leng += tail_pos + tail_leng;
	*seq = (uint32_t *) calloc((*leng + 3) / 4, sizeof(uint32_t));

	int tmp = 0;
	for (int i = head_pos - head_leng + 1; i < start_leng; ++i){
		__binseq_set(*seq, tmp, __binseq_get(lg.edges[path[0]].seq, i));
		++tmp;
	}
	for (int i = 1; i < path_leng - 1; ++i){
		struct asm_edge_t e = lg.edges[path[i]];
		for (int j = 0; j < (int) e.seq_len; ++j){
			__binseq_set(*seq, tmp, __binseq_get(e.seq, j));
			++tmp;
		}
	}
	for (int i = 0; i < tail_pos + tail_leng; ++i){
		__binseq_set(*seq, tmp,
			__binseq_get(lg.edges[path[path_leng - 1]].seq, i));
		++tmp;
	}
}*/

void get_midair_bridge(struct asm_graph_t lg, int start_edge, int end_edge,
		int *is_disable, int *midair_edge)
{
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
	*midair_edge = res;
}

int get_bridge(struct asm_graph_t *g, struct asm_graph_t *lg, int e1, int e2,
		uint32_t **ret_seq, uint32_t *seq_len)
{
	int lc_e1, pos1;
	get_local_edge_id_head(*lg, g->edges[e1], &lc_e1, &pos1);

	int lc_e2, pos2;
	get_local_edge_id_tail(*lg, g->edges[e2], &lc_e2, &pos2);
	__VERBOSE_LOG("", "Local edge 1: %d, pos: %d\n", lc_e1, pos1);
	__VERBOSE_LOG("", "Local edge 2: %d, pos: %d\n", lc_e2, pos2);

	int bridge_type;
	if (lc_e1 == -1 || lc_e2 == -1){
		goto path_not_found;
	} else if (lc_e1 == lc_e2){
		pos2 = max(pos2, pos1 + 1);
		bridge_type = 1;
		*seq_len = g->edges[e1].seq_len + g->edges[e2].seq_len
				+ pos2 - pos1 - 1;
		*ret_seq = (uint32_t *) calloc((*seq_len + 3) / 4, 
						sizeof(uint32_t));

		int tmp = 0;
		for (int i = 0; i < (int) g->edges[e1].seq_len; ++i){
			__binseq_set(*ret_seq, tmp,
					__binseq_get(g->edges[e1].seq, i));
			tmp++;
		}
		for (int i = pos1 + 1; i < pos2; ++i){
			__binseq_set(*ret_seq, tmp,
					__binseq_get(lg->edges[lc_e1].seq, i));
			tmp++;
		}
		for (int i = 0; i < (int) g->edges[e2].seq_len; ++i){
			__binseq_set(*ret_seq, tmp,
					__binseq_get(g->edges[e2].seq, i));
			tmp++;
		}
		goto path_found;
	} else if (lg->edges[lc_e1].rc_id == lc_e2){
		__VERBOSE_LOG("", "Same edge but in reverse, continuing\n");
		goto path_not_found;
	} else {
		__VERBOSE("Printing graph %d %d\n", e1, e2);
		int *is_disable;
		filter_edges(*lg, lc_e1, lc_e2, &is_disable);
		char file_path[10000];
		sprintf(file_path, "local_assembly_filtered_%d_%d.gfa", e1, e2);
		print_graph(*lg, is_disable, file_path);
		int *path;
		int path_leng;
		get_path(*lg, lc_e1, lc_e2, &path, &path_leng);
		if (path == NULL){
			goto path_not_found;
		} else {
			__VERBOSE("Path exists, checking for simple path\n");
			if (check_simple_path(*lg, lc_e1, lc_e2)){
				bridge_type = 1;
				__VERBOSE("Simple path\n");
				combine_edges(*lg, ret_seq, seq_len,
						g->edges[e1].seq_len, pos1,
						g->edges[e2].seq_len, pos2,
						path, path_leng);
			} else {
				bridge_type = 0;
				__VERBOSE("Complex path\n");
				int midair_edge;
				get_midair_bridge(*lg, lc_e1, lc_e2, is_disable,
						&midair_edge);
				if (midair_edge != -1){
					*ret_seq = lg->edges[midair_edge].seq;
					*seq_len = lg->edges[midair_edge].seq_len;
				} else {
					*ret_seq = NULL;
					*seq_len = 0;
				}
				/*combine_edges(*lg, ret_seq, seq_len,
						g->edges[e1].seq_len, pos1,
						g->edges[e2].seq_len, pos2,
						path, path_leng);*/
			}
			goto path_found;
		}
	}
path_not_found:
	*ret_seq = NULL;
	*seq_len = 0;
	bridge_type = -1;
path_found:
	return bridge_type;
}
