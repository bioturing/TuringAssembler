#include "graph_search.h"
#include "verbose.h"
#include "helper.h"
#include <math.h>

void graph_info_init(struct asm_graph_t *lg, struct graph_info_t *ginfo,
			int lc_e1, int lc_e2)
{
	ginfo->g = lg;
	ginfo->lc_e1 = lc_e1;
	ginfo->lc_e2 = lc_e2;
	ginfo->is_edge_trash = (int *) calloc(lg->n_e, sizeof(int));
	ginfo->edge_vst_count = (int *) calloc(lg->n_e, sizeof(int));
	ginfo->is_link_trash = kh_init(gint_int);
	ginfo->edge_max_vst = NULL;
	graph_info_init_max_vst(ginfo);
}

void graph_info_init_max_vst(struct graph_info_t *ginfo)
{
	if (ginfo->edge_max_vst != NULL)
		free(ginfo->edge_max_vst);
	ginfo->edge_max_vst = (int *) calloc(ginfo->g->n_e, sizeof(int));
	float init_cov = (float) (get_cov(*ginfo->g, ginfo->lc_e1) +
			get_cov(*ginfo->g, ginfo->lc_e2)) / 2;
	for (int i = 0; i < ginfo->g->n_e; ++i){
		int cov = get_cov(*(ginfo->g), i);
		ginfo->edge_max_vst[i] = (int) max(1, round(1.0 * cov / init_cov));
	}
}

void copy_static_info(struct graph_info_t *dest, struct graph_info_t *source)
{
	dest->g = source->g;
	dest->lc_e1 = source->lc_e1;
	dest->lc_e2 = source->lc_e2;

	dest->is_edge_trash = (int *) calloc(dest->g->n_e, sizeof(int));
	memcpy(dest->is_edge_trash, source->is_edge_trash, sizeof(int) *
			dest->g->n_e);
	dest->edge_vst_count = (int *) calloc(dest->g->n_e, sizeof(int));

	dest->is_link_trash = kh_init(gint_int);
	for (khiter_t it = kh_begin(source->is_link_trash);
		it != kh_end(source->is_link_trash); ++it){
		gint_t key = kh_key(source->is_link_trash, it);
		insert_key(dest->is_link_trash, key);
	}
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
	free(ginfo->edge_vst_count);
	free(ginfo->edge_max_vst);
	kh_destroy(gint_int, ginfo->is_link_trash);
}

gint_t get_edge_code(gint_t u, gint_t v)
{
	return (u << 32) | v;
}

int check_edge_visited(struct graph_info_t *ginfo, int e)
{
	return ginfo->edge_vst_count[e] == ginfo->edge_max_vst[e];
}

void mark_edge_visited(struct graph_info_t *ginfo, int e)
{
	++ginfo->edge_vst_count[e];
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
	--ginfo->edge_vst_count[e];
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

void print_graph(struct asm_graph_t *lg, int lc_e1, int lc_e2)
{
	char path[1024];
	sprintf(path, "%d_%d_local.gfa", lc_e1, lc_e2);
	FILE *f = fopen(path, "w");
	int *mark = (int *) calloc(lg->n_e, sizeof(int));
	for (int u = 0; u < lg->n_e; ++u){
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
	for (int it = 0; it < lg->n_e; ++it){
		int u = it;
		if (!mark[u])
			continue;
		int tg = lg->edges[u].target;
		for (int i = 0; i < (int) lg->nodes[tg].deg; ++i){
			int v = lg->nodes[tg].adj[i];
			int new_u;
			char dir1;
			if (mark[u] == 1){
				dir1 = '+';
				new_u = u;
			} else {
				dir1 = '-';
				new_u = lg->edges[u].rc_id;
			}
			int new_v;
			char dir2;
			if (mark[v] == 1){
				dir2 = '+';
				new_v = v;
			} else {
				dir2 = '-';
				new_v = lg->edges[v].rc_id;
			}
			fprintf(f, "L\t%d_%ld_%.3f\t%c\t%d_%ld_%.3f\t%c\t%dM\n",
					new_u, lg->edges[new_u].rc_id,
					get_cov(*lg, new_u), dir1,
					new_v, lg->edges[new_v].rc_id,
					get_cov(*lg, new_v), dir2,
					lg->ksize);
		}
	}
	free(mark);
	fclose(f);
}

void get_all_paths(struct asm_graph_t *lg, struct edge_map_info_t *emap1,
		struct edge_map_info_t *emap2, struct path_info_t *pinfo)
{
	struct graph_info_t ginfo;
	graph_info_init(lg, &ginfo, emap1->lc_e, emap2->lc_e);
	mark_edge_trash(&ginfo, emap1->lc_e);
	mark_edge_trash(&ginfo, lg->edges[emap1->lc_e].rc_id);
	mark_edge_trash(&ginfo, lg->edges[emap2->lc_e].rc_id);
	int path[1024];
	find_all_paths(lg, &ginfo, emap1->lc_e, 0, path, pinfo);
	graph_info_destroy(&ginfo);
}

void get_all_paths_kmer_check(struct asm_graph_t *lg, struct edge_map_info_t *emap1,
		struct edge_map_info_t *emap2, struct path_info_t *pinfo,
		int ksize, khash_t(kmer_int) *h)
{
	struct graph_info_t ginfo;
	graph_info_init(lg, &ginfo, emap1->lc_e, emap2->lc_e);
	mark_edge_trash(&ginfo, emap1->lc_e);
	mark_edge_trash(&ginfo, lg->edges[emap1->lc_e].rc_id);
	mark_edge_trash(&ginfo, lg->edges[emap2->lc_e].rc_id);
	int path[1024];
	find_all_paths_kmer_check(lg, &ginfo, emap1->lc_e, 0, path, pinfo,
			ksize, h);
	graph_info_destroy(&ginfo);
}

void find_all_paths(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int u, int depth, int *cur_path, struct path_info_t *pinfo)
{
	if (pinfo->n_paths == MAX_PATH_COUNT)
		return;
	cur_path[depth] = u;
	if (u == ginfo->lc_e2){
		path_info_push(pinfo, cur_path, depth + 1);
		return;
	}
	int tg = lg->edges[u].target;
	for (int i = 0; i < lg->nodes[tg].deg; ++i){
		gint_t v = lg->nodes[tg].adj[i];
		if (check_edge_trash(ginfo, v))
			continue;
		if (check_edge_visited(ginfo, v))
			continue;
		//__VERBOSE("%d %d %d\n", u, v, depth);
		mark_edge_visited(ginfo, v);
		find_all_paths(lg, ginfo, v, depth + 1, cur_path, pinfo);
		unmark_edge_visited(ginfo, v);
	}
}

void find_all_paths_kmer_check(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int u, int depth, int *cur_path, struct path_info_t *pinfo,
		int ksize, khash_t(kmer_int) *h)
{
	if (pinfo->n_paths == MAX_PATH_COUNT)
		return;
	cur_path[depth] = u;
	if (u == ginfo->lc_e2){
		path_info_push(pinfo, cur_path, depth + 1);
		__VERBOSE("%d paths found\n", pinfo->n_paths);
		return;
	}
	int tg = lg->edges[u].target;
	char *first;
	decode_seq(&first, lg->edges[u].seq, lg->edges[u].seq_len);
	for (int i = 0; i < lg->nodes[tg].deg; ++i){
		gint_t v = lg->nodes[tg].adj[i];
		if (check_edge_trash(ginfo, v))
			continue;
		if (check_edge_visited(ginfo, v))
			continue;
		char *second;
		decode_seq(&second, lg->edges[v].seq, lg->edges[v].seq_len);
		int kmer_res = kmer_check(first, second, lg->ksize,
				ksize, h);
		/*int zero = count_zero_kmer_map(first, second, lg->ksize,
				ksize, h);*/
		int max_con = count_max_consecutive_zero_kmer(first, second,
				lg->ksize, ksize, h);
		//__VERBOSE_LOG("", "%d %d %d\n", u, v, zero);
		//printf("%d %d %d %d %d\n", depth, lg->nodes[tg].deg, u, v, max_con);
		free(second);
		//if (lg->nodes[tg].deg > 1 && max_con >= 1)
		if (lg->nodes[tg].deg > 1 && max_con >= 1)
		//if (kmer_res == 0)
			continue;
		mark_edge_visited(ginfo, v);
		find_all_paths_kmer_check(lg, ginfo, v, depth + 1, cur_path,
				pinfo, ksize, h);
		unmark_edge_visited(ginfo, v);
	}
	free(first);
}

int get_path(struct asm_graph_t *lg, int lc_e1, int lc_e2,
		int *middle_edge, int **path, int *path_len)
{
	int res;
	__VERBOSE_LOG("PATH", "Finding path from %d to %d\n", lc_e1,
			lc_e2);
	struct graph_info_t ginfo;
	graph_info_init(lg, &ginfo, lc_e1, lc_e2);
	//print_graph(lg, &ginfo); //DEBUG only
	*path_len = 0;
	*path = NULL;
	int found = find_path_hao(lg, &ginfo, lc_e1, 0, path, path_len);
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

void bfs(struct asm_graph_t *lg, struct graph_info_t *ginfo, int lc_e1,
		int **bfs_len)
{
	int *queue = (int *) calloc(ginfo->g->n_e, sizeof(int));
	*bfs_len = (int *) calloc(ginfo->g->n_e, sizeof(int));
	for (int i = 0; i < (int) ginfo->g->n_e; ++i)
		(*bfs_len)[i] = -1;
	(*bfs_len)[lc_e1] = 0;
	int front = 0, back = 1;
	queue[0] = lc_e1;
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
	if (u == ginfo->lc_e2){
		*path_len = depth + 1;
		*path = (int *) calloc(*path_len, sizeof(int));
		goto path_found;
	}
	int tg = lg->edges[u].target;
	for (int i = 0; i < lg->nodes[tg].deg; ++i){
		gint_t v = lg->nodes[tg].adj[i];
		if (v == ginfo->lc_e1)
			continue;
		if (lg->edges[v].rc_id == ginfo->lc_e2)
			continue;
		if (check_edge_trash(ginfo, v))
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
		int tmp = find_path_hao(lg, ginfo, ginfo->lc_e1, 0, new_path,
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
	find_middle_edge_candidates(lg, &new_ginfo, new_ginfo.lc_e1, path,
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
	if (u == ginfo->lc_e2){
		for (int i = 1; i < depth; ++i)
			mark[path[i]] = 1;
		return;
	}
	int tg = lg->edges[u].target;
	for (int i = 0; i < lg->nodes[tg].deg; ++i){
		int v = lg->nodes[tg].adj[i];
		if (v == lg->edges[ginfo->lc_e1].rc_id ||
			v == lg->edges[ginfo->lc_e2].rc_id ||
			v == ginfo->lc_e1)
			continue;
		if (check_edge_trash(ginfo, v))
			continue;
		if (check_link_trash(ginfo, u, v))
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

void get_nearby_edges(struct asm_graph_t *g, int e, struct graph_info_t *ginfo,
		int radius, int **res, int *n_nb)
{
	*res = (int *) calloc(g->n_e, sizeof(int));
	*n_nb = 0;
	int *dis = (int *) calloc(g->n_e, sizeof(int));
	for (int i = 0; i < g->n_e; ++i)
		dis[i] = -1;
	dis[e] = 0;
	int *queue = (int *) calloc(g->n_e, sizeof(int));
	int fr = 0, bk = 1;
	queue[0] = e;
	while (fr < bk){
		int u = queue[fr++];
		(*res)[*n_nb] = u;
		++(*n_nb);
		if (dis[u] == radius)
			continue;
		int tg = g->edges[u].target;
		for (int i = 0; i < g->nodes[tg].deg; ++i){
			int v = g->nodes[tg].adj[i];
			if (check_edge_trash(ginfo, v))
				continue;
			if (dis[v] == -1){
				dis[v] = dis[u] + 1;
				queue[bk++] = v;
			}
		}
	}
	free(dis);
	free(queue);
	*res = (int *) realloc(*res, *n_nb * sizeof(int));
}

