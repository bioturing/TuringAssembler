#include "graph_search.h"
#include "verbose.h"
#include "helper.h"
#include <math.h>

void graph_info_init(struct asm_graph_t *lg, struct graph_info_t *ginfo,
			int start_edge, int end_edge)
{
	ginfo->g = lg;
	ginfo->start_edge = start_edge;
	ginfo->end_edge = end_edge;
	ginfo->is_edge_trash = (int *) calloc(lg->n_e, sizeof(int));
	ginfo->is_edge_vst = (int *) calloc(lg->n_e, sizeof(int));
	ginfo->is_link_trash = kh_init(gint_int);
	ginfo->is_link_vst = kh_init(gint_int);
	graph_info_init_max_vst(ginfo);
}

void graph_info_init_max_vst(struct graph_info_t *ginfo)
{
	ginfo->link_max_vst = kh_init(gint_int);
	float avg_cov = 0;
	int n_e = 0;
	/*int *mark = (int *) calloc(ginfo->g->n_e, sizeof(int));
	for (int i = 0; i < ginfo->g->n_e; ++i){
		if (check_edge_trash(ginfo, i))
			continue;
		if (mark[i])
			continue;
		mark[i] = mark[ginfo->g->edges[i].rc_id] = 1;
		++n_e;
		avg_cov += get_cov(*(ginfo->g), i);
	}
	avg_cov /= n_e;*/
	float init_cov = (float) (get_cov(*ginfo->g, ginfo->start_edge) +
			get_cov(*ginfo->g, ginfo->end_edge)) / 2;
	for (int u = 0; u < ginfo->g->n_e; ++u){
		int tg = ginfo->g->edges[u].target;
		for (int i = 0; i < (int) ginfo->g->nodes[tg].deg; ++i){
			int v = ginfo->g->nodes[tg].adj[i];
			gint_t edge_code = get_edge_code(u, v);
			int ret;
			khiter_t it = kh_put(gint_int, ginfo->link_max_vst,
					edge_code, &ret);
			float tmp = ceil(1.0 * get_cov(*ginfo->g, v) / init_cov);
			/*if (tmp > avg_cov)
				tmp = avg_cov;*/
			kh_val(ginfo->link_max_vst, it) = (int) tmp;
		}
	}
}

void copy_static_info(struct graph_info_t *dest, struct graph_info_t *source)
{
	dest->g = source->g;
	dest->start_edge = source->start_edge;
	dest->end_edge = source->end_edge;

	dest->is_edge_trash = (int *) calloc(dest->g->n_e, sizeof(int));
	memcpy(dest->is_edge_trash, source->is_edge_trash, sizeof(int) *
			dest->g->n_e);
	dest->is_edge_vst = (int *) calloc(dest->g->n_e, sizeof(int));

	dest->is_link_trash = kh_init(gint_int);
	for (khiter_t it = kh_begin(source->is_link_trash);
		it != kh_end(source->is_link_trash); ++it){
		gint_t key = kh_key(source->is_link_trash, it);
		insert_key(dest->is_link_trash, key);
	}
	dest->is_link_vst = kh_init(gint_int);
	dest->link_max_vst = kh_init(gint_int);
	graph_info_init_max_vst(dest);
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
	kh_destroy(gint_int, ginfo->link_max_vst);
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
	khiter_t it = kh_get(gint_int, ginfo->is_link_vst, edge_code);
	int n_vst;
	if (it == kh_end(ginfo->is_link_vst))
		n_vst = 0;
	else
		n_vst = kh_val(ginfo->is_link_vst, it);
	int max_vst;
	it = kh_get(gint_int, ginfo->link_max_vst, edge_code);
	if (it == kh_end(ginfo->link_max_vst))
		max_vst = 0;
	else
		max_vst = kh_val(ginfo->link_max_vst, it);
	return n_vst == max_vst;
	//return check_key_exist(ginfo->is_link_vst, edge_code);
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
	--kh_val(h, it);
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
	int *mark = (int *) calloc(ginfo->g->n_e, sizeof(int));
	for (int u = 0; u < ginfo->g->n_e; ++u){
		if (check_edge_trash(ginfo, u))
			continue;
		if (mark[u])
			continue;
		mark[u] = 1;
		mark[lg->edges[u].rc_id] = 2;
		fprintf(f, "S\t%d_%ld_%.3f\t", u, lg->edges[u].rc_id,
				get_cov(*lg, u));
		for (int i = 0; i < (int) lg->edges[u].seq_len; ++i)
			fprintf(f, "%c", int_to_base(__binseq_get(lg->edges[u].seq, i)));
		fprintf(f, "\tKC:i:%ld\n", lg->edges[u].count);
	}
	for (int it = 0; it < ginfo->g->n_e; ++it){
		int u = it;
		if (!mark[u])
			continue;
		int tg = lg->edges[u].target;
		for (int i = 0; i < (int) lg->nodes[tg].deg; ++i){
			int v = lg->nodes[tg].adj[i];
			if (check_edge_trash(ginfo, v))
				continue;
			if (check_link_trash(ginfo, u, v))
				continue;
			char dir1;
			if (mark[u] == 1){
				dir1 = '+';
			} else {
				dir1 = '-';
				u = lg->edges[u].rc_id;
			}
			char dir2;
			if (mark[v] == 1){
				dir2 = '+';
			} else {
				dir2 = '-';
				v = lg->edges[v].rc_id;
			}
			fprintf(f, "L\t%d_%ld_%.3f\t%c\t%d_%ld_%.3f\t%c\t%dM\n",
				u, lg->edges[u].rc_id, get_cov(*lg, u), dir1,
				v, lg->edges[v].rc_id, get_cov(*lg, v), dir2,
				lg->ksize);
		}
	}
	fclose(f);
}

void get_all_paths(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		struct path_info_t *pinfo)
{
	filter_edges(lg, ginfo);
	print_graph(lg, ginfo);
	int deg_sum = 0;
	for (int i = 0; i < ginfo->g->n_v; ++i)
		deg_sum += ginfo->g->nodes[i].deg;
	int *path = (int *) calloc(deg_sum, sizeof(int));
	find_all_paths(lg, ginfo, ginfo->start_edge, 0, path, pinfo);
	free(path);
}

void find_all_paths(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int u, int depth, int *cur_path, struct path_info_t *pinfo)
{
	if (pinfo->n_paths == MAX_PATH_COUNT)
		return;
	cur_path[depth] = u;
	if (u == ginfo->end_edge){
		path_info_push(pinfo, cur_path, depth + 1);
		return;
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
		find_all_paths(lg, ginfo, v, depth + 1, cur_path, pinfo);
		unmark_link_visited(ginfo, u, v);
	}
}

int get_path(struct asm_graph_t *lg, int start_edge, int end_edge,
		int *middle_edge, int **path, int *path_len)
{
	int res;
	__VERBOSE_LOG("PATH", "Finding path from %d to %d\n", start_edge,
			end_edge);
	struct graph_info_t ginfo;
	graph_info_init(lg, &ginfo, start_edge, end_edge);
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
	graph_info_init_max_vst(ginfo);
}

void cov_filter(struct asm_graph_t *lg, struct graph_info_t *ginfo)
{
	__VERBOSE_LOG("COV FILTER", "+----------------------------------+\n");
	int is_disabled = 0;
	for (int i = 0; i < lg->n_e; ++i)
		is_disabled += check_edge_trash(ginfo, i);
	__VERBOSE_LOG("", "Before filter: %ld edges\n", lg->n_e - is_disabled);
	int thresh = (int) (MIN_DEPTH_RATIO *
			max(get_cov(*lg, ginfo->start_edge),
				get_cov(*lg, ginfo->end_edge)));
	int disabled = 0;
	for (int i = 0; i < ginfo->g->n_e; ++i){
		if (get_cov(*lg, i) < thresh){
			mark_edge_trash(ginfo, i);
			++disabled;
		}
	}
	__VERBOSE_LOG("", "After filter: %d edges\n",
			ginfo->g->n_e - disabled);
}

void link_filter(struct asm_graph_t *lg, struct graph_info_t *ginfo)
{
	for (int u = 0; u < ginfo->g->n_e; ++u){
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
	for (int i = 0; i < ginfo->g->n_e; ++i)
		disabled += check_edge_trash(ginfo, i);
	__VERBOSE_LOG("", "Before filter: %d edges\n",
			ginfo->g->n_e - disabled);
	int *forward_len, *backward_len;
	bfs(lg, ginfo, ginfo->start_edge, &forward_len);
	bfs(lg, ginfo, lg->edges[ginfo->end_edge].rc_id, &backward_len);
	for (int i = 0; i < ginfo->g->n_e; ++i){
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
			ginfo->g->n_e - disabled);
}

void bfs(struct asm_graph_t *lg, struct graph_info_t *ginfo, int start_edge,
		int **bfs_len)
{
	int *queue = (int *) calloc(ginfo->g->n_e, sizeof(int));
	*bfs_len = (int *) calloc(ginfo->g->n_e, sizeof(int));
	for (int i = 0; i < (int) ginfo->g->n_e; ++i)
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
	int *mark = (int *) calloc(ginfo->g->n_e, sizeof(int));
	find_middle_edge_candidates(lg, &new_ginfo, new_ginfo.start_edge, path,
			0, mark);
	for (int i = 0; i < ginfo->g->n_e; ++i){
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

void path_info_init(struct path_info_t *path)
{
	path->n_paths = 0;
	path->m_paths = 1;
	path->path_lens = (int *) calloc(1, sizeof(int));
	path->paths = (int **) calloc(1, sizeof(int *));
}

void path_info_push(struct path_info_t *pinfo, int *path, int len)
{
	if (pinfo->n_paths == pinfo->m_paths){
		pinfo->m_paths <<= 1;
		pinfo->paths = (int **) realloc(pinfo->paths,
				sizeof(int *) * pinfo->m_paths);
		pinfo->path_lens = (int *) realloc(pinfo->path_lens,
				sizeof(int) * pinfo->m_paths);
	}
	pinfo->path_lens[pinfo->n_paths] = len;
	pinfo->paths[pinfo->n_paths] = (int *) calloc(len,
			sizeof(int));
	memcpy(pinfo->paths[pinfo->n_paths], path, sizeof(int) * len);
	++pinfo->n_paths;
}

void path_info_destroy(struct path_info_t *pinfo)
{
	free(pinfo->path_lens);
	for (int i = 0; i < pinfo->n_paths; ++i)
		free(pinfo->paths[i]);
	free(pinfo->paths);
}
