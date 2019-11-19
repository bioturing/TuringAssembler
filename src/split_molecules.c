#include "split_molecules.h"
#include "verbose.h"
#define MAX_EDGE_VISITED 100000
#define MAX_PATH_LEN 2

void count_hits_on_edges(struct mm_hits_t *hits)
{
	/*for (khiter_t it = kh_begin(hits->aln); it != kh_end(hits->aln); ++it){
		if (!kh_exist(hits->aln, it))
			continue;
		int e = kh_key(hits->aln, it);
		int 
	}*/
}

void init_line_graph(struct line_graph_t *lig, struct asm_graph_t *g,
		struct mm_hits_t *hits)
{
	khash_t(set_int) *mark = kh_init(set_int);
	for (khiter_t it = kh_begin(hits->edges); it != kh_end(hits->edges);
			++it){
		if (!kh_exist(hits->edges, it))
			continue;
		int e = kh_key(hits->edges, it);
		//__VERBOSE("e=%d\n", e);
		int rc = g->edges[e].rc_id;
		if (g->edges[e].seq_len < 500)
			continue;
		/*if (g->edges[e].seq_len > 50000)
			continue;*/
		put_in_set(mark, e);
		put_in_set(mark, rc);
	}
	lig->line_v = calloc(kh_size(mark), sizeof(int));
	lig->vertices = kh_init(edge_line);
	for (khiter_t it = kh_begin(mark); it != kh_end(mark); ++it){
		if (!kh_exist(mark, it))
			continue;
		int e = kh_key(mark, it);
		lig->line_v[lig->n_v++] = e;
		int ret;
		struct line_vertex_t *lv = calloc(1, sizeof(struct line_vertex_t));
		khiter_t it2 = kh_put(edge_line, lig->vertices, e, &ret);
		kh_val(lig->vertices, it2) = lv;
	}
	kh_destroy(set_int, mark);
}

void construct_line_graph(struct asm_graph_t *g, struct line_graph_t *lig)
{
	for (int i = 0; i < lig->n_v; ++i){
		int e = lig->line_v[i];
		khash_t(set_int) *nearby = kh_init(set_int);
		get_edges_in_radius(g, e, nearby);
		//__VERBOSE("Nearby of %d: ", e);
		for (int j = 0; j < lig->n_v; ++j){
			int next_e = lig->line_v[j];
			if (next_e != e && check_in_set(nearby, next_e)){
				//__VERBOSE("%d ", next_e);
				add_line_edge(lig, e, next_e);
			}
		}
		//__VERBOSE("\n");
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

void get_edges_in_radius_dfs(struct asm_graph_t *g, int e, int len,
		khash_t(set_int) *visited, khash_t(set_int) *nearby)
{
	if (len > MAX_PATH_LEN)
		return;
	if (check_in_set(visited, e))
		return;
	put_in_set(visited, e);
	put_in_set(nearby, e);
	//__VERBOSE("%d %d\n", e, len);
	int target = g->edges[e].target;
	for (int i = 0; i < g->nodes[target].deg; ++i){
		int next_e = g->nodes[target].adj[i];
		get_edges_in_radius_dfs(g, next_e, len + 1, visited, nearby);
	}
	erase_from_set(visited, e);
}

void get_edges_in_radius(struct asm_graph_t *g, int e, khash_t(set_int) *nearby)
{
	struct queue_t q;
	init_queue(&q, 1024);
	push_queue(&q, pointerize(&e, sizeof(int)));
	khash_t(int_int) *L = kh_init(int_int);
	put_in_map(L, e, 0);
	//__VERBOSE("BFS FROM %d\n", e);
	khash_t(int_int) *P = kh_init(int_int);
	put_in_map(P, e, -1);
	while (!is_queue_empty(&q)){
		int *cur_e = get_queue(&q);
		pop_queue(&q);
		int tg = g->edges[*cur_e].target;
		int cur_len = get_in_map(L, *cur_e);
		put_in_set(nearby, *cur_e);
		//__VERBOSE("%d %d\n", *cur_e, cur_len);
		if (cur_len < MAX_PATH_LEN){
			for (int i = 0; i < g->nodes[tg].deg; ++i){
				int next_e = g->nodes[tg].adj[i];
				if (check_in_map(L, next_e) == 0){
					put_in_map(L, next_e, cur_len + 1);
					put_in_map(P, next_e, *cur_e);
					push_queue(&q, pointerize(&next_e, sizeof(int)));
				}
			}
		}
		free(cur_e);
	}
	/*if (e == 633){
		__VERBOSE("path: ");
		for (int i = 1559; i != -1; i = get_in_map(P, i)){
			__VERBOSE("%d, ", i);
		}
		exit(0);
	}*/
	kh_destroy(int_int, L);
	destroy_queue(&q);
}

//void get_edges_in_radius(struct asm_graph_t *g, int e, khash_t(set_int) *nearby)
//{
//	khash_t(set_int) *visited = kh_init(set_int);
//	get_edges_in_radius_dfs(g, e, 0, visited, nearby);
//	kh_destroy(set_int, visited);
//}

void order_edges(struct opt_proc_t *opt, struct asm_graph_t *g, struct mm_hits_t *hits)
{
	FILE *f = fopen(opt->lc, "a");
	struct line_graph_t *lig = calloc(1, sizeof(struct line_graph_t));
	init_line_graph(lig, g, hits);
	construct_line_graph(g, lig);
	for (khiter_t it = kh_begin(lig->vertices); it != kh_end(lig->vertices);
			++it){
		if (!kh_exist(lig->vertices, it))
			continue;
		int e = kh_key(lig->vertices, it);
		struct line_vertex_t *lv = kh_val(lig->vertices, it);
		if (lv->deg_in != 0)
			continue;
		int *tmp = calloc(lig->n_v, sizeof(int));
		int n = 0;
		//__VERBOSE("start at %d:\n", e);
		while (lv->deg_out == 1){
			tmp[n++] = e;
			e = lv->children[0];
			khiter_t it2 = kh_get(edge_line, lig->vertices, e);
			lv = kh_val(lig->vertices, it2);
			if (lv->deg_in != 1)
				break;
			//__VERBOSE("%d %d\n", e, lv->deg_out);
		}
		if (lv->deg_out == 0 && lv->deg_in == 1){
			tmp[n++] = e;
			if (n > 1){
				fprintf(f, "%s\n", opt->bx_str);
				for (int i = 0; i < n; ++i){
					fprintf(f, "%d ", tmp[i]);
					__VERBOSE("%d ", tmp[i]);
				}
				__VERBOSE("\n");
				fprintf(f, "\n");
			}
		}
	}
}

