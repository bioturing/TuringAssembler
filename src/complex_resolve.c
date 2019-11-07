#include "complex_resolve.h"
#include "khash.h"
#include "log.h"
#include "resolve.h"
#include "verbose.h"

KHASH_SET_INIT_INT(set_int);

void asm_graph_to_virtual(struct asm_graph_t *g, struct virtual_graph_t *vg)
{
	vg->n_v = g->n_v;
	vg->vertices = calloc(vg->n_v, sizeof(struct vertex_t));
	for (int i = 0; i < g->n_v; ++i){
		vg->vertices[i].deg_out = g->nodes[i].deg;
		vg->vertices[i].children = calloc(vg->vertices[i].deg_out,
				sizeof(int));
		for (int j = 0; j < vg->vertices[i].deg_out; ++j){
			struct asm_edge_t e = g->edges[g->nodes[i].adj[j]];
			vg->vertices[i].children[j] = e.target;
		}

		int rc = g->nodes[i].rc_id;
		vg->vertices[i].deg_in = g->nodes[rc].deg;
		vg->vertices[i].parent = calloc(vg->vertices[i].deg_in,
				sizeof(int));
		for (int j = 0; j < vg->vertices[i].deg_in; ++j){
			struct asm_edge_t e = g->edges[g->nodes[rc].adj[j]];
			int v = g->nodes[e.target].rc_id;
			vg->vertices[i].parent[j] = v;
		}
	}

	vg->exist = calloc(vg->n_v, sizeof(int));
	for (int i = 0; i < vg->n_v; ++i)
		vg->exist[i] = 1;
}

void virtual_graph_destroy(struct virtual_graph_t *vg)
{
	for (int i = 0; i < vg->n_v; ++i){
		free(vg->vertices[i].children);
		free(vg->vertices[i].parent);
	}
	free(vg->vertices);
	free(vg->exist);
}

//void bfs_by_vertice(struct asm_graph_t *g, int v, int **path_len)
//{
//	int *q = calloc(g->n_v, sizeof(int));
//	int front = 0;
//	int back = 0;
//	q[0] = v;
//	*path_len = calloc(g->n_v, sizeof(int));
//	for (int i = 0; i < g->n_v; ++i)
//		(*path_len)[i] = -1;
//	(*path_len)[v] = 0;
//	while (front <= back){
//		int v = q[front++];
//		for (int i = 0; i < g->nodes[v].deg; ++i){
//			int e = g->nodes[v].adj[i];
//			int u = g->edges[e].target;
//			if ((*path_len)[u] != -1)
//				continue;
//			(*path_len)[u] = (*path_len[v]) + 1;
//			q[back++] = u;
//		}
//	}
//	free(q);
//}
//
//void get_dominated_vertices(struct virtual_graph_t *vg, int v, int **dom, int *n_dom)
//{
//	int *deg_in = calloc(vg->n_v, sizeof(int));
//	int *is_v_parents = calloc(vg->n_v, sizeof(int));
//	for (int i = 0; i < vg->n_v; ++i){
//		for (int j = 0; j < vg->vertices[i].deg; ++j){
//			int u = vg->vertices[i].adj[j];
//			++deg_in[u];
//			if (u == v)
//				is_v_parents[i] = 1;
//		}
//	}
//
//	khash_t(set_int) *h = kh_init(set_int);
//	int ret;
//	kh_put(set_int, h, v, &ret);
//	*dom = calloc(vg->n_v, sizeof(int));
//	*n_dom = 0;
//	while (kh_size(h) > 0){
//		int v = -1;
//		for (khiter_t it = kh_begin(h); it != kh_end(h); ++it){
//			if (kh_exist(h, it)){
//				v = kh_key(h, it);
//				kh_del(set_int, h, it);
//				break;
//			}
//		}
//		if (v == -1)
//			log_error("Something went wrong");
//		(*dom)[(*n_dom)++] = v;
//		for (int i = 0; i < vg->vertices[v].deg; ++i){
//			int u = vg->vertices[v].adj[i];
//			--deg_in[u];
//			if (deg_in[u] == 0 && is_v_parents[u] == 0)
//				kh_put(set_int, h, u, &ret);
//		}
//	}
//	*dom = realloc(*dom, *n_dom * sizeof(int));
//	kh_destroy(set_int, h);
//	free(is_v_parents);
//	free(deg_in);
//}
//
//void get_closure(struct asm_graph_t *g, int *B_cand, int n_cand, int **B, int *n_B)
//{
//}
//
//void asm_resolve_complex_bulges_ite(struct opt_proc_t *opt, struct asm_graph_t *g)
//{
//	struct virtual_graph_t vg;
//	asm_graph_to_virtual(g, &vg);
//	int *dom;
//	int n_dom;
//	int v = opt->lk;
//	get_dominated_vertices(&vg, v, &dom, &n_dom);
//	/*for (int i = 0; i < n_dom; ++i)
//		__VERBOSE("%d ", dom[i]);
//	__VERBOSE("\n");*/
//	int *mark = calloc(g->n_v, sizeof(int));
//	for (int i = 0; i < n_dom; ++i){
//		/*__VERBOSE("node %d deg %d\n", dom[i], g->nodes[dom[i]].deg);
//		for (int j = 0; j < g->nodes[dom[i]].deg; ++j)
//			__VERBOSE("%d->%d\n", g->nodes[dom[i]].adj[j],
//					g->edges[g->nodes[dom[i]].adj[j]].target);*/
//		mark[dom[i]] = 1;
//	}
//	/*__VERBOSE("\n%d\n", g->edges[266711].target);
//	for (int i = 0; i < g->nodes[v].deg; ++i){
//		__VERBOSE("%d ", g->nodes[v].adj[i]);
//	}
//	__VERBOSE("\n");*/
//	int *keep = calloc(g->n_e, sizeof(int));
//	for (int i = 0; i < g->n_e; ++i){
//		int u = g->edges[i].source;
//		int v = g->edges[i].target;
//		if (u == -1)
//			continue;
//		if (mark[u] && mark[v])
//			keep[i] = keep[g->edges[i].rc_id] = 1;
//	}
//	for (int i = 0; i < g->n_e; ++i){
//		int rc = g->edges[i].rc_id;
//		if (i > rc)
//			continue;
//		if (!keep[i]){
//			asm_remove_edge(g, i);
//			asm_remove_edge(g, rc);
//		}
//	}
//	free(dom);
//	free(mark);
//	free(keep);
//}
