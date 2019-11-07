#ifndef __COMPLEX_RESOLVE__
#define __COMPLEX_RESOLVE__
#include "assembly_graph.h"

struct fixed_size_queue_t{
	int *q;
	int front;
	int back;
};

void init_queue(struct fixed_size_queue_t *fq, int size);
void push_queue(struct fixed_size_queue_t *fq, int v);
int get_queue(struct fixed_size_queue_t *fq);
void pop_queue(struct fixed_size_queue_t *fq);
int is_queue_empty(struct fixed_size_queue_t *fq);
void destroy_queue(struct fixed_size_queue_t *fq);

struct vertex_t{
	int deg_in;
	int deg_out;
	int *child;
	int *parent;
};

struct virtual_graph_t {
	int n_v;
	struct vertex_t *vertices;
	int *f;
};

void asm_graph_to_virtual(struct asm_graph_t *g, struct virtual_graph_t *vg);
void clone_virtual_graph(struct virtual_graph_t *org, struct virtual_graph_t *clone);
void virtual_graph_destroy(struct virtual_graph_t *vg);
//void bfs_by_vertice(struct asm_graph_t *g, int v, int **path_len);
void get_dominated_vertices(struct virtual_graph_t *vg, int v,
		struct virtual_graph_t *dom);
int get_closure(struct virtual_graph_t *vg, int s, struct virtual_graph_t *B);
void asm_resolve_complex_bulges_ite(struct opt_proc_t *opt, struct asm_graph_t *g);
void add_vertex_to_B_dfs(struct virtual_graph_t *vg, int v, struct virtual_graph_t *B,
		int *in_queue, struct fixed_size_queue_t *interested);
#endif
