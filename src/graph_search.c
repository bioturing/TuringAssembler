#include "graph_search.h"
#include "verbose.h"
#include "helper.h"

void graph_info_init(struct graph_info_t *ginfo, int n_edges, int start_edge,
		int end_edge)
{
	ginfo->n_edges = n_edges;
	ginfo->start_edge = start_edge;
	ginfo->end_edge = end_edge;
	ginfo->is_edge_trash = (int *) calloc(n_edges, sizeof(int));
	ginfo->is_edge_vst = (int *) calloc(n_edges, sizeof(int));
	ginfo->is_link_trash = kh_init(gint_int);
	ginfo->is_link_vst = kh_init(gint_int);
}

void copy_static_info(struct graph_info_t *dest, struct graph_info_t *source)
{
	dest->n_edges = source->n_edges;
	dest->start_edge = source->start_edge;
	dest->end_edge = source->end_edge;

	dest->is_edge_trash = (int *) calloc(dest->n_edges, sizeof(int));
	memcpy(dest->is_edge_trash, source->is_edge_trash, sizeof(int) *
			dest->n_edges);
	dest->is_edge_vst = (int *) calloc(dest->n_edges, sizeof(int));

	dest->is_link_trash = kh_init(gint_int);
	for (khiter_t it = kh_begin(source->is_link_trash);
		it != kh_end(source->is_link_trash); ++it){
		gint_t key = kh_key(source->is_link_trash, it);
		insert_key(dest->is_link_trash, key);
	}
	dest->is_link_vst = kh_init(gint_int);
}

int check_edge_trash(struct graph_info_t *ginfo, int e)
{
	return ginfo->is_edge_trash[e];
}

int check_link_trash(struct graph_info_t *ginfo, int e1, int e2)
{
	gint_t edge_code = get_edge_code(e1, e2);
	return check_key_exist(ginfo->is_link_trash, edge_code);
}

void graph_info_destroy(struct graph_info_t *ginfo)
{
	free(ginfo->is_edge_trash);
	free(ginfo->is_edge_vst);
	kh_destroy(gint_int, ginfo->is_link_trash);
	kh_destroy(gint_int, ginfo->is_link_vst);
}

gint_t get_edge_code(gint_t u, gint_t v)
{
	return (u << 32) | v;
}

int check_edge_visited(struct graph_info_t *ginfo, int e)
{
	return ginfo->is_edge_vst[e];
}

int check_link_visited(struct graph_info_t *ginfo, int e1, int e2)
{
	gint_t edge_code = get_edge_code(e1, e2);
	return check_key_exist(ginfo->is_link_vst, edge_code);
}

void mark_edge_visited(struct graph_info_t *ginfo, int e)
{
	ginfo->is_edge_vst[e] = 1;
}

void mark_edge_trash(struct graph_info_t *ginfo, int e)
{
	ginfo->is_edge_trash[e] = 1;
}

void mark_link_trash(struct graph_info_t *ginfo, int e1, int e2)
{
	gint_t edge_code = get_edge_code(e1, e2);
	insert_key(ginfo->is_link_trash, edge_code);
}

void mark_link_visited(struct graph_info_t *ginfo, int e1, int e2)
{
	gint_t edge_code = get_edge_code(e1, e2);
	insert_key(ginfo->is_link_vst, edge_code);
}

void unmark_link_visited(struct graph_info_t *ginfo, int e1, int e2)
{
	gint_t edge_code = get_edge_code(e1, e2);
	remove_key(ginfo->is_link_vst, edge_code);
}

void unmark_edge_visited(struct graph_info_t *ginfo, int e)
{
	ginfo->is_edge_vst[e] = 0;
}

void remove_key(khash_t(gint_int) *h, gint_t key)
{
	khiter_t it = kh_get(gint_int, h, key);
	if (it != kh_end(h))
		kh_del(gint_int, h, it);
}

void insert_key(khash_t(gint_int) *h, gint_t key)
{
	khiter_t it = kh_get(gint_int, h, key);
	if (it == kh_end(h)){
		int ret;
		it = kh_put(gint_int, h, key, &ret);
		kh_val(h, it) = 0;
	}
	++kh_val(h, it);
}

int check_key_exist(khash_t(gint_int) *h, gint_t key)
{
	khiter_t it = kh_get(gint_int, h, key);
	return it != kh_end(h);
}

void print_graph(struct asm_graph_t *lg, struct graph_info_t *ginfo)
{
	char path[1024];
	sprintf(path, "%d_%d_local.gfa", ginfo->start_edge, ginfo->end_edge);
	FILE *f = fopen(path, "w");
	int *mark = (int *) calloc(ginfo->n_edges, sizeof(int));
	for (int u = 0; u < ginfo->n_edges; ++u){
		if (check_edge_trash(ginfo, u))
			continue;
		if (mark[u] || mark[lg->edges[u].rc_id])
			continue;
		mark[u] = mark[lg->edges[u].rc_id] = 1;
		fprintf(f, "S\t%d_%ld_%.3f\t", u, lg->edges[u].rc_id,
				get_cov(*lg, u));
		for (int i = 0; i < (int) lg->edges[u].seq_len; ++i)
			fprintf(f, "%c", int_to_base(__binseq_get(lg->edges[u].seq, i)));
		fprintf(f, "\tKC:i:%ld\n", lg->edges[u].count);
	}
	for (int u = 0; u < ginfo->n_edges; ++u)
		mark[u] = 0;
	for (int u = 0; u < ginfo->n_edges; ++u){
		if (check_edge_trash(ginfo, u))
			continue;
		if (mark[u] || mark[lg->edges[u].rc_id])
			continue;
		mark[u] = mark[lg->edges[u].rc_id] = 1;
		int tg = lg->edges[u].target;
		for (int i = 0; i < (int) lg->nodes[tg].deg; ++i){
			int v = lg->nodes[tg].adj[i];
			if (check_edge_trash(ginfo, v))
				continue;
			if (check_link_trash(ginfo, u, v))
				continue;
			fprintf(f, "L\t%d_%ld_%.3f\t+\t%d_%ld_%.3f\t+\t%dM\n", u,
					lg->edges[u].rc_id, get_cov(*lg, u),
					v, lg->edges[v].rc_id, get_cov(*lg, v),
					lg->ksize);
		}
	}
	fclose(f);
}

int get_path(struct asm_graph_t *lg, int start_edge, int end_edge,
		int *middle_edge, int **path, int *path_len)
{
	int res;
	__VERBOSE_LOG("PATH", "Finding path from %d to %d\n", start_edge,
			end_edge);
	struct graph_info_t ginfo;
	graph_info_init(&ginfo, lg->n_e, start_edge, end_edge);
	filter_edges(lg, &ginfo);
	print_graph(lg, &ginfo);
	*path_len = 0;
	*path = NULL;
	int found = find_path_hao(lg, &ginfo, start_edge, 0, path, path_len);

	if (found == 0){
		__VERBOSE_LOG("PATH", "Path not found, exit\n");
		res = PATH_NOT_FOUND;
		goto path_not_found;
	}

	__VERBOSE_LOG("PATH", "Path found, checking if the path is simple\n");
	struct graph_info_t new_ginfo;
	copy_static_info(&new_ginfo, &ginfo);
	int *new_path;
	int new_path_len;
	int simple_path = check_simple_path(lg, &new_ginfo, *path, *path_len,
			&new_path, &new_path_len);
	if (simple_path){
		__VERBOSE_LOG("PATH", "The path is simple:\n");
		print_path(*path, *path_len);
		res = SIMPLE_PATH;
	} else {
		__VERBOSE_LOG("PATH", "Found first path:\n");
		print_path(*path, *path_len);
		__VERBOSE_LOG("PATH", "found second path:\n");
		print_path(new_path, new_path_len);

		free(new_path);
		res = COMPLEX_PATH;
		__VERBOSE_LOG("PATH", "Finding middle edge\n");
		*middle_edge = get_best_middle_edge(lg, &ginfo);
	}
	graph_info_destroy(&new_ginfo);
	goto end_function;
path_not_found:
	*middle_edge = -1;
	goto end_function;
end_function:
	graph_info_destroy(&ginfo);
	return res;
}

void filter_edges(struct asm_graph_t *lg, struct graph_info_t *ginfo)
{
	cov_filter(lg, ginfo);
	//link_filter(lg, ginfo);
	connection_filter(lg, ginfo);
}

void cov_filter(struct asm_graph_t *lg, struct graph_info_t *ginfo)
{
	__VERBOSE_LOG("COV FILTER", "+----------------------------------+\n");
	__VERBOSE_LOG("", "Before filter: %ld edges\n", lg->n_e);
	int thresh = (int) (MIN_DEPTH_RATIO * 
			max(get_cov(*lg, ginfo->start_edge),
				get_cov(*lg, ginfo->end_edge)));
	int disabled = 0;
	for (int i = 0; i < ginfo->n_edges; ++i){
		if (get_cov(*lg, i) < thresh){
			mark_edge_trash(ginfo, i);
			++disabled;
		}
	}
	__VERBOSE_LOG("", "After filter: %d edges\n",
			ginfo->n_edges - disabled);
}

void link_filter(struct asm_graph_t *lg, struct graph_info_t *ginfo)
{
	for (int u = 0; u < ginfo->n_edges; ++u){
		int tg = lg->edges[u].target;
		for (int i = 0; i < lg->nodes[tg].deg; ++i){
			int v = lg->nodes[tg].adj[i];
			float a = get_cov(*lg, u);
			float b = get_cov(*lg, v);
			if (a / b < MIN_RELATIVE_COV_RATIO ||
				b / a < MIN_RELATIVE_COV_RATIO)
				mark_link_trash(ginfo, u, v);
		}
	}
}

void connection_filter(struct asm_graph_t *lg, struct graph_info_t *ginfo)
{
	__VERBOSE_LOG("CONNECTION FILTER", "+----------------------------------+\n");
	int disabled = 0;
	for (int i = 0; i < ginfo->n_edges; ++i)
		disabled += check_edge_trash(ginfo, i);
	__VERBOSE_LOG("", "Before filter: %d edges\n",
			ginfo->n_edges - disabled);
	int *forward_len, *backward_len;
	bfs(lg, ginfo, ginfo->start_edge, &forward_len);
	bfs(lg, ginfo, lg->edges[ginfo->end_edge].rc_id, &backward_len);
	for (int i = 0; i < ginfo->n_edges; ++i){
		if (check_edge_trash(ginfo, i))
			continue;
		int l1 = forward_len[i];
		int l2 = backward_len[lg->edges[i].rc_id];
		if (l1 == -1 || l2 == -1 || l1 + l2 > MIN_PATH_LENGTH){
			mark_edge_trash(ginfo, i);
			++disabled;
		}
	}
	__VERBOSE_LOG("", "After filter: %d edges\n",
			ginfo->n_edges - disabled);
}

void bfs(struct asm_graph_t *lg, struct graph_info_t *ginfo, int start_edge,
		int **bfs_len)
{
	int *queue = (int *) calloc(ginfo->n_edges, sizeof(int));
	*bfs_len = (int *) calloc(ginfo->n_edges, sizeof(int));
	for (int i = 0; i < (int) ginfo->n_edges; ++i)
		(*bfs_len)[i] = -1;
	(*bfs_len)[start_edge] = 0;
	int front = 0, back = 1;
	queue[0] = start_edge;
	while (front < back){
		int u = queue[front++];
		int tg = lg->edges[u].target;
		for (int i = 0; i < (int) lg->nodes[tg].deg; ++i){
			int v = lg->nodes[tg].adj[i];
			if (check_edge_trash(ginfo, v))
				continue;
			if (check_link_trash(ginfo, u, v))
				continue;
			if ((*bfs_len)[v] != -1)
				continue;
			(*bfs_len)[v] = (*bfs_len)[u] + 1;
			queue[back++] = v;
		}
	}
	free(queue);
}

int find_path_hao(struct asm_graph_t *lg, struct graph_info_t *ginfo, int u,
		int depth, int **path, int *path_len)
{
	*path = NULL;
	*path_len = 0;
	if (u == ginfo->end_edge){
		*path_len = depth + 1;
		*path = (int *) calloc(*path_len, sizeof(int));
		goto path_found;
	}
	int tg = lg->edges[u].target;
	for (int i = 0; i < lg->nodes[tg].deg; ++i){
		gint_t v = lg->nodes[tg].adj[i];
		if (v == ginfo->start_edge)
			continue;
		if (lg->edges[v].rc_id == ginfo->end_edge)
			continue;
		if (check_edge_trash(ginfo, v))
			continue;
		if (check_link_visited(ginfo, u, v))
			continue;
		mark_link_visited(ginfo, u, v);
		if (find_path_hao(lg, ginfo, v, depth + 1, path, path_len))
			goto path_found;
	}
	return 0;
path_found:
	(*path)[depth] = u;
	return 1;
}

int check_simple_path(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int *path, int path_len, int **new_path, int *new_path_len)
{
	int res = 1;
	for (int i = 1; i < path_len; ++i){
		kh_destroy(gint_int, ginfo->is_link_vst);
		ginfo->is_link_vst = kh_init(gint_int);
		mark_link_visited(ginfo, path[i - 1], path[i]);
		int tmp = find_path_hao(lg, ginfo, ginfo->start_edge, 0, new_path,
				new_path_len);
		if (tmp){
			res = 0;
			break;
		}
	}
	return res;
}

void print_path(int *path, int path_len)
{
	for (int i = 0; i < path_len; ++i)
		__VERBOSE_LOG("", "%d ", path[i]);
	__VERBOSE_LOG("", "\n");
}

int get_best_middle_edge(struct asm_graph_t *lg, struct graph_info_t *ginfo)
{
	struct graph_info_t new_ginfo;
	copy_static_info(&new_ginfo, ginfo);

	int res = -1;
	int deg_sum = 0;
	for (int i = 0; i < lg->n_v; ++i)
		deg_sum += lg->nodes[i].deg;
	int *path = (int *) calloc(deg_sum, sizeof(int));
	int *mark = (int *) calloc(ginfo->n_edges, sizeof(int));
	find_middle_edge_candidates(lg, &new_ginfo, new_ginfo.start_edge, path,
			0, mark);
	for (int i = 0; i < ginfo->n_edges; ++i){
		if (!mark[i] || mark[lg->edges[i].rc_id])
			continue;
		if (res == -1 || lg->edges[res].seq_len < lg->edges[i].seq_len)
			res = i;
	}
	free(path);
	free(mark);
	graph_info_destroy(&new_ginfo);
	return res;
}

void find_middle_edge_candidates(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int u, int *path, int depth, int *mark)
{
	path[depth] = u;
	if (u == ginfo->end_edge){
		for (int i = 1; i < depth; ++i)
			mark[path[i]] = 1;
		return;
	}
	int tg = lg->edges[u].target;
	for (int i = 0; i < lg->nodes[tg].deg; ++i){
		int v = lg->nodes[tg].adj[i];
		if (v == lg->edges[ginfo->start_edge].rc_id ||
			v == lg->edges[ginfo->end_edge].rc_id ||
			v == ginfo->start_edge)
			continue;
		if (check_edge_trash(ginfo, v))
			continue;
		if (check_link_trash(ginfo, u, v))
			continue;
		if (check_link_visited(ginfo, u, v))
			continue;
		mark_link_visited(ginfo, u, v);
		find_middle_edge_candidates(lg, ginfo, v, path, depth + 1, mark);
		unmark_link_visited(ginfo, u, v);
	}
}
