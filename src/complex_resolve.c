#include "complex_resolve.h"
#include "khash.h"
#include "log.h"
#include "resolve.h"
#include "verbose.h"

KHASH_SET_INIT_INT(set_int);

void init_queue(struct fixed_size_queue_t *fq, int size)
{
	fq->q = calloc(size, sizeof(int));
	fq->front = 0;
	fq->back = 0;
}

void push_queue(struct fixed_size_queue_t *fq, int v)
{
	fq->q[fq->back++] = v;
}

int get_queue(struct fixed_size_queue_t *fq)
{
	return fq->q[fq->front];
}

void pop_queue(struct fixed_size_queue_t *fq)
{
	fq->front++;
}

int is_queue_empty(struct fixed_size_queue_t *fq)
{
	return fq->front >= fq->back;
}

void destroy_queue(struct fixed_size_queue_t *fq)
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

	vg->f = calloc(vg->n_v, sizeof(int));
	for (int i = 0; i < vg->n_v; ++i)
		vg->f[i] = i;
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
	clone->f = calloc(clone->n_v, sizeof(int));
	memcpy(clone->f, org->f, sizeof(int) * clone->n_v);
}

void virtual_graph_destroy(struct virtual_graph_t *vg)
{
	for (int i = 0; i < vg->n_v; ++i){
		free(vg->vertices[i].child);
		free(vg->vertices[i].parent);
	}
	free(vg->vertices);
	free(vg->f);
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

void get_dominated_vertices(struct virtual_graph_t *vg, int v,
		struct virtual_graph_t *dom)
{
	int *deg_in = calloc(vg->n_v, sizeof(int));
	for (int i = 0; i < vg->n_v; ++i)
		deg_in[i] = vg->vertices[i].deg_in;

	int *is_v_parents = calloc(vg->n_v, sizeof(int));
	for (int i = 0; i < vg->vertices[v].deg_in; ++i){
		int p = vg->vertices[v].parent[i];
		is_v_parents[p] = 1;
	}

	struct fixed_size_queue_t *queue = calloc(1, sizeof(struct fixed_size_queue_t));
	init_queue(queue, vg->n_v);
	push_queue(queue, v);

	clone_virtual_graph(vg, dom);
	for (int i = 0; i < dom->n_v; ++i)
		dom->f[i] = -1;
	while (is_queue_empty(queue) == 0){
		int v = get_queue(queue);
		pop_queue(queue);
		dom->f[v] = v;
		for (int i = 0; i < vg->vertices[v].deg_out; ++i){
			int u = vg->vertices[v].child[i];
			--deg_in[u];
			if (deg_in[u] == 0 && is_v_parents[u] == 0)
				push_queue(queue, u);
		}
	}
	destroy_queue(queue);
	free(queue);
	free(is_v_parents);
	free(deg_in);
}

void add_vertex_to_B_dfs(struct virtual_graph_t *vg, int v, struct virtual_graph_t *B,
		int *in_queue, struct fixed_size_queue_t *interested)
{
	if (!in_queue[v]){
		in_queue[v] = 1;
		push_queue(interested, v);
	}
	if (B->f[v] != -1)
		return;
	B->f[v] = v;
	for (int i = 0; i < vg->vertices[i].deg_in; ++i){
		int p = vg->vertices[v].parent[i];
		add_vertex_to_B_dfs(vg, p, B, in_queue, interested);
	}
}

int get_closure(struct virtual_graph_t *vg, int s, struct virtual_graph_t *B)
{
	int res = 0;
	struct fixed_size_queue_t interested;
	init_queue(&interested, vg->n_v);
	
	int *in_queue = calloc(vg->n_v, sizeof(int));
	for (int i = 0; i < vg->n_v; ++i){
		if (B->f[i] == -1)
			continue;
		for (int j = 0; j < vg->vertices[i].deg_out; ++j){
			int v = vg->vertices[i].child[j];
			if (B->f[v] != -1){
				in_queue[i] = 1;
				push_queue(&interested, i);
				break;
			}
		}
	}

	struct virtual_graph_t dom;
	get_dominated_vertices(vg, s, &dom);

	while (is_queue_empty(&interested) == 0){
		int v = get_queue(&interested);
		pop_queue(&interested);
		for (int i = 0; i < vg->vertices[v].deg_out; ++i){
			int u = vg->vertices[v].child[i];
			if (dom.f[u] == -1){
				res = 0;
				goto end_function;
			}
			add_vertex_to_B_dfs(vg, u, B, in_queue, &interested);
		}
	}
	res = 1;
end_function:
	virtual_graph_destroy(&dom);
	free(in_queue);
	destroy_queue(&interested);
	return res;
}

void asm_resolve_complex_bulges_ite(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	struct virtual_graph_t vg;
	asm_graph_to_virtual(g, &vg);

	struct virtual_graph_t dom;
	int v = opt->lk;
	get_dominated_vertices(&vg, v, &dom);
	int *mark = calloc(g->n_v, sizeof(int));
	for (int i = 0; i < dom.n_v; ++i){
		if (dom.f[i] != -1){
			mark[i] = 1;
		}
	}
	int *keep = calloc(g->n_e, sizeof(int));
	for (int i = 0; i < g->n_e; ++i){
		int u = g->edges[i].source;
		int v = g->edges[i].target;
		if (u == -1)
			continue;
		if (mark[u] && mark[v])
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
	free(mark);
	virtual_graph_destroy(&dom);
	virtual_graph_destroy(&vg);
}
