#include "complex_resolve.h"
#include "khash.h"
#include "log.h"
#include "resolve.h"
#include "verbose.h"
#include "helper.h"
#include "dqueue.h"

void init_queue(struct fixed_size_q_t *fq, int size)
{
	fq->q = calloc(size, sizeof(int));
	fq->front = 0;
	fq->back = 0;
}

void push_queue(struct fixed_size_q_t *fq, int v)
{
	fq->q[fq->back++] = v;
}

int get_queue(struct fixed_size_q_t *fq)
{
	return fq->q[fq->front];
}

void pop_queue(struct fixed_size_q_t *fq)
{
	fq->front++;
}

int is_queue_empty(struct fixed_size_q_t *fq)
{
	return fq->front >= fq->back;
}

void destroy_queue(struct fixed_size_q_t *fq)
{
	free(fq->q);
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
	bundle->source = -1;
	bundle->dom = calloc(g->n_v, sizeof(int));
	bundle->B = calloc(g->n_v, sizeof(int));
	bundle->S = calloc(g->n_v, sizeof(int));
	bundle->T = calloc(g->n_v, sizeof(int));
	bundle->g = calloc(g->n_v, sizeof(int));
	bundle->j = calloc(g->n_v, sizeof(int));
	bundle->height = calloc(g->n_v, sizeof(int));
}

void reset_source(struct resolve_bulges_bundle_t *bundle, int s)
{
	bundle->source = s;
	memset(bundle->dom, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->B, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->S, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->T, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->g, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->j, 0, sizeof(int) * bundle->graph->n_v);
	memset(bundle->height, 0, sizeof(int) * bundle->graph->n_v);
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

	int *deg_in = calloc(graph->n_v, sizeof(int));
	for (int i = 0; i < graph->n_v; ++i)
		deg_in[i] = graph->vertices[i].deg_in;

	int *is_s_parents = calloc(graph->n_v, sizeof(int));
	for (int i = 0; i < graph->vertices[s].deg_in; ++i){
		int p = graph->vertices[s].parent[i];
		is_s_parents[p] = 1;
	}

	struct fixed_size_q_t *queue = calloc(1, sizeof(struct fixed_size_q_t));
	init_queue(queue, graph->n_v);
	push_queue(queue, s);

	while (is_queue_empty(queue) == 0){
		int v = get_queue(queue);
		pop_queue(queue);
		dom[v] = 1;
		for (int i = 0; i < graph->vertices[v].deg_out; ++i){
			int u = graph->vertices[v].child[i];
			--deg_in[u];
			if (deg_in[u] == 0 && is_s_parents[u] == 0)
				push_queue(queue, u);
		}
	}
	destroy_queue(queue);
	free(queue);
	free(is_s_parents);
	free(deg_in);
}

void add_vertex_to_B_dfs(struct resolve_bulges_bundle_t *bundle, int v,
		int *in_queue, struct fixed_size_q_t *interested, int depth)
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
	if (int_vertex && !in_queue[v]){
		in_queue[v] = 1;
		push_queue(interested, v);
		//__VERBOSE("queue in %d\n", v);
	}
	if (bundle->B[v] != 0)
		return;
	bundle->B[v] = 1;
	//__VERBOSE("add %d\n", v);
	for (int i = 0; i < bundle->graph->vertices[v].deg_in; ++i){
		int p = bundle->graph->vertices[v].parent[i];
		add_vertex_to_B_dfs(bundle, p, in_queue, interested, depth + 1);
	}
}

int get_closure(struct resolve_bulges_bundle_t *bundle)
{
	struct virtual_graph_t *graph = bundle->graph;
	int s = bundle->source;
	int res = 0;
	struct fixed_size_q_t interested;
	init_queue(&interested, graph->n_v);
	
	int *in_queue = calloc(graph->n_v, sizeof(int));
	for (int i = 0; i < graph->n_v; ++i){
		if (bundle->B[i] == 0)
			continue;
		//__VERBOSE("%d %d\n", i, vg->vertices[i].deg_out);
		for (int j = 0; j < graph->vertices[i].deg_out; ++j){
			int v = graph->vertices[i].child[j];
			//__VERBOSE("%d %d\n", i, v);
			if (bundle->B[i]){
				in_queue[i] = 1;
				push_queue(&interested, i);
				//__VERBOSE("queue in %d\n", i);
				break;
			}
		}
	}

	while (is_queue_empty(&interested) == 0){
		int v = get_queue(&interested);
		pop_queue(&interested);
		//__VERBOSE("queue out %d\n", v);
		for (int i = 0; i < graph->vertices[v].deg_out; ++i){
			int u = graph->vertices[v].child[i];
			if (bundle->dom[u] == 0){
				res = 0;
				goto end_function;
			}
			add_vertex_to_B_dfs(bundle, u, in_queue, &interested, 0);
		}
	}
	res = 1;
end_function:
	free(in_queue);
	destroy_queue(&interested);
	return res;
}

//void print_dom_debug(struct opt_proc_t *opt, struct asm_graph_t *g)
//{
//	struct virtual_graph_t vg;
//	asm_graph_to_virtual(g, &vg);
//
//	struct virtual_graph_t dom;
//	int v = opt->lk;
//	get_dominated_vertices(&vg, v, &dom);
//	int *keep = calloc(g->n_e, sizeof(int));
//	for (int i = 0; i < g->n_e; ++i){
//		int u = g->edges[i].source;
//		int v = g->edges[i].target;
//		if (u == -1)
//			continue;
//		if (dom.f[u] != -1 && dom.f[v] != -1)
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
//	virtual_graph_destroy(&dom);
//	virtual_graph_destroy(&vg);
//}
//
//void print_closure_debug(struct opt_proc_t *opt, struct asm_graph_t *g)
//{
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
//}
//
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

int get_skeleton(struct resolve_bulges_bundle_t *bundle)
{
	get_height(bundle);
	struct virtual_graph_t *graph = bundle->graph;
	struct dqueue_t *q = init_dqueue(bundle->graph->n_v);
	for (int i = 0; i < graph->n_v; ++i)
		if (bundle->height[i] == 0){
			struct vertex_height_t *pair = calloc(1,
					sizeof(struct vertex_height_t));
			d_enqueue_in(q, pair);
		}
}

void resolve_bulges(struct asm_graph_t *g)
{
	struct resolve_bulges_bundle_t bundle;
	init_resolve_bulges(g, &bundle);
	struct virtual_graph_t *graph = bundle.graph;
	for (int i = 0; i < graph->n_v; ++i){
		reset_source(&bundle, i);
		get_dominated_vertices(&bundle);
		int closure_exist = get_closure(&bundle);
	}
}

