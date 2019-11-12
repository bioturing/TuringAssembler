#include "complex_resolve.h"
#include "khash.h"
#include "log.h"
#include "resolve.h"
#include "verbose.h"
#include "helper.h"
#include "dqueue.h"

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

void asm_graph_to_virtual(struct asm_graph_t *g, struct virtual_graph_t *vg)
{
	vg->n_v = g->n_v;
	vg->vertices = calloc(vg->n_v, sizeof(struct vertex_t));
	for (int i = 0; i < g->n_v; ++i){
		vg->vertices[i].deg_out = g->nodes[i].deg;
		vg->vertices[i].child = calloc(vg->vertices[i].deg_out,
				sizeof(int));
		for (int j = 0; j < vg->vertices[i].deg_out; ++j){
			struct asm_edge_t e = g->edges[g->nodes[i].adj[j]];
			vg->vertices[i].child[j] = e.target;
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

void clone_virtual_graph(struct virtual_graph_t *org, struct virtual_graph_t *clone)
{
	clone->n_v = org->n_v;
	clone->vertices = calloc(clone->n_v, sizeof(struct vertex_t));
	for (int i = 0; i < clone->n_v; ++i){
		clone->vertices[i].deg_out = org->vertices[i].deg_out;
		clone->vertices[i].child = calloc(clone->vertices[i].deg_out,
				sizeof(int));
		memcpy(clone->vertices[i].child, org->vertices[i].child,
				sizeof(int) * clone->vertices[i].deg_out);

		clone->vertices[i].deg_in = org->vertices[i].deg_in;
		clone->vertices[i].parent = calloc(clone->vertices[i].deg_in,
				sizeof(int));
		memcpy(clone->vertices[i].parent, org->vertices[i].parent,
				sizeof(int) * clone->vertices[i].deg_in);
	}
	clone->exist = calloc(clone->n_v, sizeof(int));
	memcpy(clone->exist, org->exist, sizeof(int) * clone->n_v);
}

void virtual_graph_destroy(struct virtual_graph_t *vg)
{
	for (int i = 0; i < vg->n_v; ++i){
		free(vg->vertices[i].child);
		free(vg->vertices[i].parent);
	}
	free(vg->vertices);
	free(vg->exist);
}

void init_resolve_bulges(struct asm_graph_t *g, struct resolve_bulges_bundle_t *bundle)
{
	bundle->graph = calloc(1, sizeof(struct virtual_graph_t));
	asm_graph_to_virtual(g, bundle->graph);

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
	bundle->P = calloc(g->n_v, sizeof(int));
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
	memset(bundle->P, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->L, 0, sizeof(int) * bundle->graph->n_v);
}

void bulges_bundle_destroy(struct resolve_bulges_bundle_t *bundle)
{
	virtual_graph_destroy(bundle->graph);
	destroy_queue(bundle->B_vertices);
	destroy_queue(bundle->dom_vertices);
	destroy_queue(bundle->closest);
	free(bundle->dom);
	free(bundle->B);
	free(bundle->S);
	free(bundle->T);
	free(bundle->g);
	free(bundle->j);
	free(bundle->height);
	free(bundle->P);
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
	struct virtual_graph_t *graph = bundle->graph;
	int s = bundle->source;
	int *dom = bundle->dom;

	int *descresed_deg = calloc(graph->n_v, sizeof(int));

	int *is_s_parents = calloc(graph->n_v, sizeof(int));
	for (int i = 0; i < graph->vertices[s].deg_in; ++i){
		int p = graph->vertices[s].parent[i];
		is_s_parents[p] = 1;
	}

	struct queue_t *queue = calloc(1, sizeof(struct queue_t));
	init_queue(queue, graph->n_v);
	push_queue(queue, pointerize(&s, sizeof(int)));

	while (is_queue_empty(queue) == 0){
		int *v = get_queue(queue);
		pop_queue(queue);
		dom[*v] = 1;
		push_queue(bundle->dom_vertices, pointerize(v, sizeof(int)));
		for (int i = 0; i < graph->vertices[*v].deg_out; ++i){
			int u = graph->vertices[*v].child[i];
			int remain_deg = --graph->vertices[u].deg_in;
			++descresed_deg[u];
			if (remain_deg == 0 && is_s_parents[u] == 0)
				push_queue(queue, pointerize(&u, sizeof(int)));
		}
		free(v);
	}
	for (int i = bundle->dom_vertices->front; i < bundle->dom_vertices->back;
			++i){
		int *v = bundle->dom_vertices->data[i];
		bundle->graph->vertices[*v].deg_in += descresed_deg[*v];
	}
	destroy_queue(queue);
	free(queue);
	free(is_s_parents);
	free(descresed_deg);
}

void add_vertex_to_B(struct resolve_bulges_bundle_t *bundle, int v)
{
	bundle->B[v] = 1;
	push_queue(bundle->B_vertices, pointerize(&v, sizeof(int)));
}

void add_vertex_to_B_dfs(struct resolve_bulges_bundle_t *bundle, int v,
		int *in_queue, struct queue_t *q, int depth)
{
	int int_vertex = 0;
	if (depth == 0){
		for (int i = 0; i < bundle->graph->vertices[v].deg_out; ++i){
			int u = bundle->graph->vertices[v].child[i];
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
	for (int i = 0; i < bundle->graph->vertices[v].deg_in; ++i){
		int p = bundle->graph->vertices[v].parent[i];
		add_vertex_to_B_dfs(bundle, p, in_queue, q, depth + 1);
	}
}

int get_closure(struct resolve_bulges_bundle_t *bundle)
{
	struct virtual_graph_t *graph = bundle->graph;
	int s = bundle->source;
	int res = 0;
	struct queue_t q;
	init_queue(&q, graph->n_v);
	
	int *in_queue = calloc(graph->n_v, sizeof(int));
	struct queue_t *B_vertices = bundle->B_vertices;
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		//__VERBOSE("%d %d\n", i, vg->vertices[i].deg_out);
		for (int i = 0; i < graph->vertices[*v].deg_out; ++i){
			int u = graph->vertices[*v].child[i];
			//__VERBOSE("%d %d\n", i, v);
			if (bundle->B[u]){
				in_queue[*v] = 1;
				push_queue(&q, pointerize(v, sizeof(int)));
				//__VERBOSE("queue in %d\n", i);
				break;
			}
		}
	}

	while (is_queue_empty(&q) == 0){
		int *v = get_queue(&q);
		pop_queue(&q);
		//__VERBOSE("queue out %d\n", v);
		for (int i = 0; i < graph->vertices[*v].deg_out; ++i){
			int u = graph->vertices[*v].child[i];
			if (bundle->dom[u] == 0){
				res = 0;
				goto end_function;
			}
			if (bundle->B[*v] != 0)
				continue;
			add_vertex_to_B_dfs(bundle, u, in_queue, &q, 0);
		}
		free(v);
	}
	res = 1;
end_function:
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
	struct virtual_graph_t *graph = bundle->graph;
	int *height = bundle->height;
	height[v] = 0;
	for (int i = 0; i < graph->vertices[v].deg_out; ++i){
		int u = graph->vertices[v].child[i];
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
	struct virtual_graph_t *graph = bundle->graph;
	for (int i = 0; i < graph->vertices[u].deg_out; ++i){
		int w = graph->vertices[u].child[i];
		if (bundle->B[w] == 0)
			continue;
		int ok = 0;
		for (int j = 0; j < graph->vertices[w].deg_in; ++j){
			int t = graph->vertices[w].parent[j];
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
	struct virtual_graph_t *graph = bundle->graph;
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
	struct virtual_graph_t *graph = bundle->graph;
	int *P = bundle->P;

	struct queue_t *q = calloc(1, sizeof(struct queue_t));
	push_queue(q, pointerize(&bundle->source, sizeof(int)));
	int *L = calloc(graph->n_v, sizeof(int));
	memset(L, -1, sizeof(int) * graph->n_v);
	L[bundle->source] = 0;


	while (!is_queue_empty(q)){
		int *v = get_queue(q);
		pop_queue(q);
		for (int i = 0; i < graph->vertices[*v].deg_out; ++i){
			int u = graph->vertices[*v].child[i];
			if (bundle->B[u] == 0)
				continue;
			if (L[u] == -1){
				L[u] = L[*v] + 1;
				P[u] = *v;
				push_queue(q, pointerize(&u, sizeof(int)));
			}
		}
		free(v);
	}
	free(L);
}

void get_tree(struct resolve_bulges_bundle_t *bundle)
{
}

void get_distance(struct resolve_bulges_bundle_t *bundle)
{
	struct virtual_graph_t *graph = bundle->graph;
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
		for (int i = 0; i < graph->vertices[*v].deg_out; ++i){
			int u = graph->vertices[*v].child[i];
			if (bundle->dom[u] == 0)
				continue;
			if (L[u] == -1){
				L[u] = L[*v] + 1;
				push_queue(q, pointerize(&u, sizeof(int)));
			}
		}
		free(v);
	}
}

int is_complex_closure(struct resolve_bulges_bundle_t *bundle)
{
	struct virtual_graph_t *graph = bundle->graph;
	struct queue_t *B_vertices = bundle->B_vertices;
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		for (int j = 0; j < graph->vertices[*v].deg_out; ++j){
			int u = graph->vertices[*v].child[j];
			if (bundle->B[u] != 0){
			}
		}
	}
}

int is_closure_tree(struct resolve_bulges_bundle_t *bundle)
{
	struct virtual_graph_t *graph = bundle->graph;
	struct queue_t *B_vertices = bundle->B_vertices;
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		int C = 0;
		for (int i = 0; i < graph->vertices[*v].deg_in; ++i){
			int w = graph->vertices[*v].parent[i];
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

void remove_edge_virtual(struct asm_graph_t *g, int e,
		struct resolve_bulges_bundle_t *bundle)
{
	int v = g->edges[e].source;
	int u = g->edges[e].target;
	asm_remove_edge(g, e);
	struct virtual_graph_t *graph = bundle->graph;
	int deg_out = graph->vertices[v].deg_out;
	for (int i = 0; i < deg_out; ++i){
		if (graph->vertices[v].child[i] == u){
			swap(graph->vertices[v].child + i,
				graph->vertices[v].child + deg_out - 1,
				sizeof(int));
			break;
		}
	}
	--(graph->vertices[v].deg_out);

	int deg_in = graph->vertices[u].deg_in;
	for (int i = 0; i < deg_in; ++i){
		if (graph->vertices[u].parent[i] == v){
			swap(graph->vertices[u].parent + i,
				graph->vertices[u].parent + deg_in - 1,
				sizeof(int));
			break;
		}
	}
	--(graph->vertices[u].deg_in);
}

int resolve_bulges(struct asm_graph_t *g)
{
	int res = 0;
	struct resolve_bulges_bundle_t bundle;
	init_resolve_bulges(g, &bundle);
	struct virtual_graph_t *graph = bundle.graph;
	for (int i = 0; i < graph->n_v; ++i){
		reset_source(&bundle, i);
		get_dominated_vertices(&bundle);
		/*if (i == 1229){
			for (int i = 0; i < bundle.dom_vertices->back; ++i)
				__VERBOSE("%d ", *(int *)bundle.dom_vertices->data[i]);
			__VERBOSE("\n");
		}*/
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
			if (is_closure_tree(&bundle))
				continue;
			bfs_to_sinks(&bundle);
			/*for (int j = bundle.B_vertices->front; j < bundle.B_vertices->back; ++j)
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
			}
			if (i == 1229){
				write_gfa(g, "./tmp.gfa");
				exit(0);
			}*/
			while (!is_queue_empty(bundle.B_vertices)){
				int *v = get_queue(bundle.B_vertices);
				pop_queue(bundle.B_vertices);
				int tmp[4];
				int C = 0;
				for (int j = 0; j < g->nodes[*v].deg; ++j){
					int e = g->nodes[*v].adj[j];
					int u = g->edges[e].target;
					int del = 0;
					if (bundle.P[u] != *v){
						del = 1;
					} else {
						for (int k = 0; k < C; ++k)
							if (tmp[k] == u)
								del = 1;
					}
					if (del){
						int rc = g->edges[e].rc_id;
						remove_edge_virtual(g, e, &bundle);
						remove_edge_virtual(g, rc, &bundle);
					} else {
						tmp[C++] = u;
					}
				}
				free(v);
			}
			log_debug("Bulge detected at %d", i);
			++res;
			break;
		}
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
	} while(1);
	log_info("%d bulge(s) resolved after %d iterations", res, ite);
	return res;
}

