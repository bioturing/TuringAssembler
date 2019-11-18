#include "split_molecules.h"
#include "verbose.h"
#define MAX_EDGE_VISITED 100000
#define MAX_EDGE_DEPTH 10
#define MAX_SEARCH_LEN 4000

void init_line_graph(struct line_graph_t *lig, int n_e, int *edges)
{
	lig->n_v = n_e;
	lig->line_v = calloc(n_e, sizeof(int));
	lig->vertices = kh_init(edge_line);
	for (int i = 0; i < n_e; ++i){
		lig->line_v[i] = edges[i];
		int ret;
		struct line_vertex_t *lv = calloc(1, sizeof(struct line_vertex_t));
		khiter_t it = kh_put(edge_line, lig->vertices, edges[i], &ret);
		kh_val(lig->vertices, it) = lv;
	}
}

void construct_line_graph(struct asm_graph_t *g, struct line_graph_t *lig)
{
	for (int i = 0; i < lig->n_v; ++i){
		int e = lig->line_v[i];
		khash_t(set_int) *nearby = kh_init(set_int);
		get_nearby_edges(g, e, nearby);
		__VERBOSE("Nearby of %d: ", e);
		for (int j = 0; j < lig->n_v; ++j){
			int next_e = lig->line_v[j];
			if (next_e != e){
				__VERBOSE("%d ", next_e);
				add_line_edge(lig, e, next_e);
			}
		}
		__VERBOSE("\n");
	}
}

void add_line_edge(struct line_graph_t *lig, int v, int u)
{
	khiter_t it = kh_get(edge_line, lig->vertices, v);
	struct line_vertex_t *lv = kh_val(lig->vertices, it);
	lv->children = realloc(lv->children, sizeof(int) * (lv->deg_out + 1));
	lv->children[lv->deg_out++] = u;

	it = kh_get(edge_line, lig->vertices, u);
	lv = kh_val(lig->vertices, it);
	lv->parents = realloc(lv->parents, sizeof(int) * (lv->deg_in + 1));
	lv->parents[lv->deg_in++] = v;
}

void get_nearby_edges_dfs(struct asm_graph_t *g, int e, int len,
		khash_t(set_int) *visited)
{
	if (len > MAX_SEARCH_LEN)
		return;
	if (check_in_set(visited, e))
		return;
	put_in_set(visited, e);
	len += g->edges[e].seq_len;
	int target = g->edges[e].target;
	for (int i = 0; i < g->nodes[target].deg; ++i){
		int next_e = g->nodes[target].adj[i];
		get_nearby_edges_dfs(g, next_e, len, visited);
	}
}

void get_nearby_edges(struct asm_graph_t *g, int e, khash_t(set_int) *nearby)
{
	get_nearby_edges_dfs(g, e, 0, nearby);
}

//int get_edges_order_dfs(struct asm_graph_t *g, int e, int p, int total, int *path,
//		int depth, khash_t(int_int) *multi, struct edge_ordering_t *order)
//{
//	//__VERBOSE("%d %d\n", e, *n_visit);
//	int res = 0;
//	if (depth > MAX_EDGE_DEPTH)
//		goto ignore_stage_1;
//	int id = g->edges[e].rc_id;
//	if (id > e)
//		id = e;
//
//	int early_stop = 1;
//	int is_edge_in_list = check_in_map(multi, id);
//	if (!is_edge_in_list || get_in_map(multi, id) > 0)
//		early_stop = 0;
//	if (early_stop)
//		goto ignore_stage_1;
//
//
//	if (is_edge_in_list){
//		__VERBOSE("HAHA %d\n", e);
//		increase_in_map(multi, id, -1);
//		path[p++] = e;
//		if (p == total){
//			order->n = total;
//			order->edges = calloc(total, sizeof(int));
//			for (int i = 0; i < total; ++i)
//				order->edges[i] = path[i];
//			res = 1;
//			goto ignore_stage_2;
//		}
//	}
//
//	int v = g->edges[e].target;
//	for (int i = 0; i < g->nodes[v].deg; ++i){
//		int next_e = g->nodes[v].adj[i];
//		if (get_edges_order_dfs(g, next_e, p, total, path, depth + 1,
//					multi, order)){
//			res = 1;
//			break;
//		}
//	}
//ignore_stage_2:
//	if (is_edge_in_list){
//		increase_in_map(multi, id, 1);
//		--p;
//	}
//ignore_stage_1:
//	return res;
//}
//
//void get_edges_order(struct asm_graph_t *g, int *edges, int n_e,
//		struct edge_ordering_t *order)
//{
//	float unit_cov = get_genome_coverage(g);
//	khash_t(int_int) *multi = kh_init(int_int);
//	int total = 0;
//	for (int i = 0; i < n_e; ++i){
//		int e = edges[i];
//		int rc = g->edges[e].rc_id;
//		if (e > rc)
//			e = rc;
//		float cov = __get_edge_cov(g->edges + e, g->ksize);
//		int v = (int) (cov / unit_cov + 0.5);
//		total += v;
//		put_in_map(multi, e, v);
//	}
//
//	int *path = calloc(g->n_e, sizeof(int));
//	for (int i = 0; i < n_e; ++i){
//		int e = edges[i];
//		__VERBOSE("START FROM %d\n", e);
//		if (get_edges_order_dfs(g, e, 0, total, path, 0, multi, order) ||
//			get_edges_order_dfs(g, g->edges[e].rc_id, 0, total, path,
//				0, multi, order)){
//			__VERBOSE("FOUND: ");
//			for (int i = 0; i < total; ++i)
//				__VERBOSE("%d ", path[i]);
//			__VERBOSE("\n");
//		}
//	}
//	free(path);
//	kh_destroy(int_int, multi);
//}

