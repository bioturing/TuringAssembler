#include "complex_resolve.h"
#include "khash.h"
#include "log.h"
#include "resolve.h"
#include "verbose.h"
#include "helper.h"
#include "dqueue.h"
#include "process.h"
#include "utils.h"

KHASH_MAP_INIT_INT(int_int, int);
KHASH_SET_INIT_INT(set_int);
void *pointerize(void *data, int size)
{
	void *res = malloc(size);
	memcpy(res, data, size);
	return res;
}

void init_queue(struct queue_t *q, int cap)
{
	q->data = calloc(cap, sizeof(void *));
	q->cap = cap;
	q->front = 0;
	q->back = 0;
}

void push_queue(struct queue_t *q, void *ptr)
{
	if (q->back == q->cap){
		q->cap = q->cap == 0 ? 1 : (q->cap << 1);
		q->data = realloc(q->data, sizeof(void *) * q->cap);
	}
	q->data[q->back++] = ptr;
}

void *get_queue(struct queue_t *q)
{
	return q->data[q->front];
}

void pop_queue(struct queue_t *fq)
{
	fq->front++;
}

int is_queue_empty(struct queue_t *q)
{
	return q->front >= q->back;
}

void destroy_queue(struct queue_t *q)
{
	free(q->data);
}

void free_queue_content(struct queue_t *q)
{
	for (int i = q->front; i < q->back; ++i)
		free(q->data[i]);
}

void init_resolve_bulges(struct asm_graph_t *g, struct resolve_bulges_bundle_t *bundle)
{
	bundle->graph = g;

	bundle->B_vertices = calloc(1, sizeof(struct queue_t));
	init_queue(bundle->B_vertices, bundle->graph->n_v);

	bundle->dom_vertices = calloc(1, sizeof(struct queue_t));
	init_queue(bundle->dom_vertices, bundle->graph->n_v);

	bundle->closest = calloc(1, sizeof(struct queue_t));
	init_queue(bundle->closest, bundle->graph->n_v);

	/*bundle->dom_vertices = calloc(1, sizeof(struct queue_t));
	init_queue(bundle->dom_vertices, bundle->graph->n_v);*/

	bundle->source = -1;
	bundle->dom = calloc(g->n_v, sizeof(int));
	bundle->B = calloc(g->n_v, sizeof(int));
	bundle->S = calloc(g->n_v, sizeof(int));
	bundle->T = calloc(g->n_v, sizeof(int));
	bundle->g = calloc(g->n_v, sizeof(int));
	bundle->j = calloc(g->n_v, sizeof(int));
	bundle->height = calloc(g->n_v, sizeof(int));
	bundle->PE = calloc(g->n_v, sizeof(int));
	bundle->L = calloc(g->n_v, sizeof(int));
}

void reset_source(struct resolve_bulges_bundle_t *bundle, int s)
{
	bundle->B_vertices->front = bundle->B_vertices->back = 0;
	bundle->dom_vertices->front = bundle->dom_vertices->back = 0;
	bundle->closest->front = bundle->closest->back = 0;
	bundle->source = s;
	memset(bundle->dom, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->B, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->S, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->T, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->g, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->j, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->height, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->PE, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->L, 0, sizeof(int) * bundle->graph->n_v);
}

void bulges_bundle_destroy(struct resolve_bulges_bundle_t *bundle)
{
	destroy_queue(bundle->B_vertices);
	free(bundle->B_vertices);
	destroy_queue(bundle->dom_vertices);
	free(bundle->dom_vertices);
	destroy_queue(bundle->closest);
	free(bundle->closest);
	free(bundle->dom);
	free(bundle->B);
	free(bundle->S);
	free(bundle->T);
	free(bundle->g);
	free(bundle->j);
	free(bundle->height);
	free(bundle->PE);
	free(bundle->L);
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

void get_dominated_vertices(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	int s = bundle->source;
	int *dom = bundle->dom;

	int *is_s_parents = calloc(graph->n_v, sizeof(int));
	int s_rc = graph->nodes[s].rc_id;
	for (int i = 0; i < graph->nodes[s_rc].deg; ++i){
		int e = graph->edges[graph->nodes[s_rc].adj[i]].rc_id;
		int p = graph->edges[e].source;
		is_s_parents[p] = 1;
	}

	struct queue_t *queue = calloc(1, sizeof(struct queue_t));
	init_queue(queue, graph->n_v);
	push_queue(queue, pointerize(&s, sizeof(int)));

	khash_t(int_int) *deg_in = kh_init(int_int);
	while (is_queue_empty(queue) == 0){
		int *v = get_queue(queue);
		pop_queue(queue);
		dom[*v] = 1;
		push_queue(bundle->dom_vertices, pointerize(v, sizeof(int)));
		for (int i = 0; i < graph->nodes[*v].deg; ++i){
			int u = get_adj_node(graph, *v, i);
			//__VERBOSE("%d %d\n", *v, u);
			khiter_t it = kh_get(int_int, deg_in, u);
			if (it == kh_end(deg_in)){
				int ret;
				it = kh_put(int_int, deg_in, u, &ret);
				kh_val(deg_in, it) = 0;
			}
			int used_deg = ++kh_val(deg_in, it);
			int u_rc = graph->nodes[u].rc_id;
			//__VERBOSE("%d %d %d\n", u, kh_val(deg_in, it), graph->nodes[u_rc].deg);
			if (used_deg == graph->nodes[u_rc].deg &&
					is_s_parents[u] == 0)
				push_queue(queue, pointerize(&u, sizeof(int)));
		}
		free(v);
	}
	kh_destroy(int_int, deg_in);
	destroy_queue(queue);
	free(queue);
	free(is_s_parents);

}

void add_vertex_to_B(struct resolve_bulges_bundle_t *bundle, int v)
{
	//__VERBOSE("Add %d\n", v);
	bundle->B[v] = 1;
	push_queue(bundle->B_vertices, pointerize(&v, sizeof(int)));
}

void add_vertex_to_B_dfs(struct resolve_bulges_bundle_t *bundle, int v,
		int *in_queue, struct queue_t *q, int depth)
{
	struct asm_graph_t *graph = bundle->graph;
	int int_vertex = 0;
	//__VERBOSE("vertex %d\n", v);
	if (depth == 0){
		for (int i = 0; i < graph->nodes[v].deg; ++i){
			int u = get_adj_node(graph, v, i);
			if (bundle->B[u] != 0)
				int_vertex = 1;
		}
	} else {
		int_vertex = 1;
	}
	if (int_vertex && in_queue[v] == 0){
		in_queue[v] = 1;
		push_queue(q, pointerize(&v, sizeof(int)));
		//__VERBOSE("queue in %d\n", v);
	}
	if (bundle->B[v] != 0)
		return;
	add_vertex_to_B(bundle, v);
	//__VERBOSE("add %d\n", v);
	int v_rc = graph->nodes[v].rc_id;
	for (int i = 0; i < graph->nodes[v_rc].deg; ++i){
		int pe = graph->edges[graph->nodes[v_rc].adj[i]].rc_id;
		int p = graph->edges[pe].source;
		add_vertex_to_B_dfs(bundle, p, in_queue, q, depth + 1);
	}
}

int get_closure(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	int s = bundle->source;
	int res = 1;
	struct queue_t q;
	init_queue(&q, graph->n_v);
	
	int *in_queue = calloc(graph->n_v, sizeof(int));
	struct queue_t *B_vertices = bundle->B_vertices;
	/*__VERBOSE("before: ");
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		__VERBOSE("%d ", *v);
	}
	__VERBOSE("\n");*/
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		//__VERBOSE("%d %d\n", i, vg->vertices[i].deg_out);
		for (int i = 0; i < graph->nodes[*v].deg; ++i){
			int u = get_adj_node(graph, *v, i);
			//__VERBOSE("%d %d\n", i, v);
			if (bundle->B[u]){
				in_queue[*v] = 1;
				push_queue(&q, pointerize(v, sizeof(int)));
				//__VERBOSE("queue in %d\n", *v);
				break;
			}
		}
	}

	while (res && is_queue_empty(&q) == 0){
		int *v = get_queue(&q);
		pop_queue(&q);
		//__VERBOSE("queue out %d\n", *v);
		for (int i = 0; i < graph->nodes[*v].deg; ++i){
			int u = get_adj_node(graph, *v, i);
			//__VERBOSE("hihi %d %d\n", *v, u);
			if (bundle->dom[u] == 0){
				res = 0;
				break;
			}
			if (bundle->B[u] != 0)
				continue;
			add_vertex_to_B_dfs(bundle, u, in_queue, &q, 0);
		}
		free(v);
	}
	free_queue_content(&q);
	free(in_queue);
	destroy_queue(&q);
	return res;
}

void print_dom_debug(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	int v = opt->lk;
	struct resolve_bulges_bundle_t bundle;
	init_resolve_bulges(g, &bundle);
	reset_source(&bundle, v);
	get_dominated_vertices(&bundle);
	int *keep = calloc(g->n_e, sizeof(int));
	for (int i = 0; i < g->n_e; ++i){
		int u = g->edges[i].source;
		int v = g->edges[i].target;
		if (u == -1)
			continue;
		if (bundle.dom[u] != 0 && bundle.dom[v] != 0)
			keep[i] = keep[g->edges[i].rc_id] = 1;
	}
	for (int i = 0; i < g->n_e; ++i){
		int rc = g->edges[i].rc_id;
		if (i > rc)
			continue;
		if (!keep[i]){
			asm_remove_edge(g, i);
			asm_remove_edge(g, rc);
		}
	}
	free(keep);
	bulges_bundle_destroy(&bundle);
}

void print_closure_debug(struct opt_proc_t *opt, struct asm_graph_t *g)
{
//	struct virtual_graph_t vg;
//	asm_graph_to_virtual(g, &vg);
//	struct virtual_graph_t B;
//	clone_virtual_graph(&vg, &B);
//	for (int i = 0; i < B.n_v; ++i)
//		B.f[i] = -1;
//
//	FILE *f = fopen(opt->in_fasta, "r");
//	int n;
//	fscanf(f, "%d\n", &n);
//	for (int i = 0; i < n; ++i){
//		int v;
//		fscanf(f, "%d", &v);
//		B.f[v] = v;
//	}
//	fclose(f);
//
//	int s = opt->lk;
//	int ret = get_closure(&vg, s, &B);
//	for (int i = 0; i < B.n_v; ++i)
//		if (B.f[i] != -1)
//			__VERBOSE("%d ", i);
//	int *keep = calloc(g->n_e, sizeof(int));
//	for (int i = 0; i < g->n_e; ++i){
//		int u = g->edges[i].source;
//		int v = g->edges[i].target;
//		if (u == -1)
//			continue;
//		if (B.f[u] != -1 && B.f[v] != -1)
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
//	free(keep);
//	virtual_graph_destroy(&B);
//	virtual_graph_destroy(&vg);
}

void get_height_dfs(struct resolve_bulges_bundle_t *bundle, int v)
{
	struct asm_graph_t *graph = bundle->graph;
	int *height = bundle->height;
	height[v] = 0;
	for (int i = 0; i < graph->nodes[v].deg; ++i){
		int u = get_adj_node(graph, v, i);
		if (bundle->B[u] == 0)
			continue;
		if(height[u] == -1)
			get_height_dfs(bundle, u);
		height[v] = max(height[v], height[u] + 1);
	}
}

void get_height(struct resolve_bulges_bundle_t *bundle)
{
	memset(bundle->height, -1, sizeof(int) * bundle->graph->n_v);
	get_height_dfs(bundle, bundle->source);
}

int is_able_to_map(struct resolve_bulges_bundle_t *bundle, int u, int v)
{
	if (bundle->height[u] == 0)
		return 0;
	struct asm_graph_t *graph = bundle->graph;
	for (int i = 0; i < graph->nodes[u].deg; ++i){
		int w = get_adj_node(graph, u, i);
		if (bundle->B[w] == 0)
			continue;
		int ok = 0;
		int w_rc = graph->nodes[w].rc_id;
		for (int j = 0; j < graph->nodes[w_rc].deg; ++j){
			int pe = graph->edges[graph->nodes[w_rc].adj[j]].rc_id;
			int t = graph->edges[pe].source;
			if (t == v){
				ok = 1;
				break;
			}
		}
		if (!ok)
			return 0;
	}
	return 1;
}

int get_skeleton(struct resolve_bulges_bundle_t *bundle)
{
	get_height(bundle);
	struct asm_graph_t *graph = bundle->graph;
	struct queue_t *q = calloc(1, sizeof(struct queue_t));
	init_queue(q, graph->n_v);
	for (int i = 0; i < graph->n_v; ++i){
		if (bundle->height[i] == 0){
			struct vertex_height_t pair = {
				.vertex = i,
				.height = 0
			};
			push_queue(q, pointerize(&pair, sizeof(struct vertex_height_t)));
		}
	}

	int *height = bundle->height;
	while (!is_queue_empty(q)){
		int p = q->front;
		while (p < q->back &&
			((struct vertex_height_t *)q->data[p])->height ==
			((struct vertex_height_t *)q->data[q->front])->height)
			++p;
		int n = p - q->front;
		for (int i = 0; i < n - 1; ++i){
			struct vertex_height_t *u = q->data[q->front + i];
			if (u->height != height[u->vertex]){
				struct vertex_height_t new_u = {
					.vertex = u->vertex,
					.height = u->height - 1
				};
				//int pr = 
			}
			for (int j = i + 1; j < n; ++j){
			}
		}
	}
	free(q);
}

void bfs_to_sinks(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	int *PE = bundle->PE;

	struct queue_t *q = calloc(1, sizeof(struct queue_t));
	push_queue(q, pointerize(&bundle->source, sizeof(int)));
	int *L = calloc(graph->n_v, sizeof(int));
	memset(L, -1, sizeof(int) * graph->n_v);
	L[bundle->source] = 0;


	while (!is_queue_empty(q)){
		int *v = get_queue(q);
		pop_queue(q);
		for (int i = 0; i < graph->nodes[*v].deg; ++i){
			int e = graph->nodes[*v].adj[i];
			int u = get_adj_node(graph, *v, i);
			if (bundle->B[u] == 0)
				continue;
			if (L[u] == -1){
				L[u] = L[*v] + 1;
				PE[u] = e;
				push_queue(q, pointerize(&u, sizeof(int)));
			}
		}
		free(v);
	}
	destroy_queue(q);
	free(q);
	free(L);
}

void get_tree(struct resolve_bulges_bundle_t *bundle)
{
}

void get_distance(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	int *L = bundle->L;

	struct queue_t *q = calloc(1, sizeof(struct queue_t));
	init_queue(q, graph->n_v);
	push_queue(q, pointerize(&bundle->source, sizeof(int)));
	memset(L, -1, sizeof(int) * graph->n_v);
	L[bundle->source] = 0;

	while (!is_queue_empty(q)){
		int *v = get_queue(q);
		pop_queue(q);
		push_queue(bundle->closest, pointerize(v, sizeof(int)));
		for (int i = 0; i < graph->nodes[*v].deg; ++i){
			int u = get_adj_node(graph, *v, i);
			if (bundle->dom[u] == 0)
				continue;
			if (L[u] == -1){
				L[u] = L[*v] + 1;
				push_queue(q, pointerize(&u, sizeof(int)));
			}
		}
		free(v);
	}
	destroy_queue(q);
	free(q);
}

int is_complex_closure(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	struct queue_t *B_vertices = bundle->B_vertices;
	int res = 0;
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		for (int j = 0; j < graph->nodes[*v].deg; ++j){
			int e = graph->nodes[*v].adj[j];
			int u = graph->edges[e].target;
			if (bundle->B[u] != 0)
				res = MAX(res, graph->edges[e].seq_len);
		}
	}
	return res >= 1000;
}

int is_closure_tree(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	struct queue_t *B_vertices = bundle->B_vertices;
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		//__VERBOSE("HAHA %d\n", *v);
		int C = 0;
		int v_rc = graph->nodes[*v].rc_id;
		for (int i = 0; i < graph->nodes[v_rc].deg; ++i){
			int pe = graph->edges[graph->nodes[v_rc].adj[i]].rc_id;
			int w = graph->edges[pe].source;
			//__VERBOSE("HAHA %d %d\n", *v, w);
			if (bundle->B[w] != 0)
				++C;
		}
		if (C > 1)
			return 0;
	}
	return 1;
}

int get_next_B_candidate(struct resolve_bulges_bundle_t *bundle)
{
	int res = -1;
	struct queue_t *closest = bundle->closest;
	while (res == -1 && !is_queue_empty(closest)){
		int *v = get_queue(closest);
		pop_queue(closest);
		if (bundle->B[*v] == 0)
			res = *v;
		free(v);
	}
	return res;
}

void supress_bulge(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	struct queue_t *B_vertices = bundle->B_vertices;
	int *mark = calloc(graph->n_v, sizeof(int));
	mark[bundle->source] = 1;
	//__VERBOSE("sink: ");
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		int is_sink = 1;
		for (int j = 0; j < graph->nodes[*v].deg; ++j){
			int u = get_adj_node(graph, *v, j);
			if (bundle->B[u] != 0){
				is_sink = 0;
				break;
			}
		}
		if (is_sink){
			//__VERBOSE("%d ", *v);
			int w = *v;
			while (mark[w] == 0){
				mark[w] = 1;
				int e = bundle->PE[w];
				w = graph->edges[e].source;
			}
		}
	}
	//__VERBOSE("\n");
	khash_t(set_int) *rm_edges = kh_init(set_int);
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		for (int j = 0; j < graph->nodes[*v].deg; ++j){
			int u = get_adj_node(graph, *v, j);
			int e = graph->nodes[*v].adj[j];
			int rc = graph->edges[e].rc_id;
			//__VERBOSE("HEHE %d %d %d %d\n", *v, u, mark[*v], mark[u]);
			if (bundle->B[u] == 0)
				continue;
			//__VERBOSE("HAHA %d %d %d %d\n", *v, u, mark[*v], mark[u]);
			if (!mark[*v] || !mark[u] ||
				(bundle->PE[u] != e && bundle->PE[u] != rc)){
				int ret;
				kh_put(set_int, rm_edges, e, &ret);
				kh_put(set_int, rm_edges, rc, &ret);
				//__VERBOSE("remove %d %d\n", e, rc);
			}
		}
	}
	for (khiter_t it = kh_begin(rm_edges); it != kh_end(rm_edges); ++it){
		if (!kh_exist(rm_edges, it))
			continue;
		int e = kh_key(rm_edges, it);
		//__VERBOSE("remove %d\n", e);
		asm_remove_edge(graph, e);
	}
	kh_destroy(set_int, rm_edges);
	free(mark);
}

int resolve_bulges(struct asm_graph_t *g)
{
	int res = 0;
	struct resolve_bulges_bundle_t bundle;
	init_resolve_bulges(g, &bundle);
	struct asm_graph_t *graph = bundle.graph;
	for (int i = 0; i < graph->n_v; ++i){
		/*if (i != 36027)
			continue;*/
		reset_source(&bundle, i);
		get_dominated_vertices(&bundle);
		/*__VERBOSE("dom: ");
		for (int i = 0; i < bundle.dom_vertices->back; ++i)
			__VERBOSE("%d ", *(int *)bundle.dom_vertices->data[i]);
		__VERBOSE("\n");*/
		get_distance(&bundle);

		add_vertex_to_B(&bundle, i);
		free(get_queue(bundle.closest));
		pop_queue(bundle.closest);
		while(1){
			int next_cand = get_next_B_candidate(&bundle);
			if (next_cand == -1)
				break;

			add_vertex_to_B(&bundle, next_cand);
			int ret = get_closure(&bundle);
			if (!ret)
				break;
			if (is_complex_closure(&bundle))
				break;
			if (is_closure_tree(&bundle))
				continue;
			bfs_to_sinks(&bundle);
			/*__VERBOSE("B\n");
			for (int j = bundle.B_vertices->front; j < bundle.B_vertices->back; ++j)
				__VERBOSE("%d ", *(int *)bundle.B_vertices->data[j]);
			__VERBOSE("\n");*/
			/*for (int i = 0; i < g->n_v; ++i){
				for (int j = 0; j < g->nodes[i].deg; ++j){
					int e = g->nodes[i].adj[j];
					int sr = g->edges[e].source;
					int tg = g->edges[e].target;
					if (bundle.B[sr] == 0 || bundle.B[tg] == 0){
						int rc = g->edges[e].rc_id;
						asm_remove_edge(g, rc);
						asm_remove_edge(g, e);
					}
				}
			}*/
			supress_bulge(&bundle);
			log_debug("Bulge detected at %d", i);
			++res;
			break;
		}
		free_queue_content(bundle.B_vertices);
		free_queue_content(bundle.dom_vertices);
		free_queue_content(bundle.closest);
	}
	bulges_bundle_destroy(&bundle);
	test_asm_graph(g);
	return res;
}

int asm_resolve_complex_bulges_ite(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	int ite = 0;
	int res = 0;
	do{
		int resolved = resolve_bulges(g);
		if (!resolved)
			break;
		res += resolved;
		++ite;
		log_debug("%d-th iteration: %d bulge(s) resolved", ite, resolved);
		char path[1024];
		sprintf(path, "%d_ite", ite);
		save_graph_info(opt->out_dir, g, path);
		struct asm_graph_t g1;
		asm_condense(g, &g1);
		asm_graph_destroy(g);
		*g = g1;

	} while(1);
	log_info("%d bulge(s) resolved after %d iterations", res, ite);
	return res;
}

int get_adj_node(struct asm_graph_t *g, int v, int id)
{
	int e = g->nodes[v].adj[id];
	return g->edges[e].target;
}

