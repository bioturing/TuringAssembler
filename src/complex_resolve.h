#ifndef __COMPLEX_RESOLVE__
#define __COMPLEX_RESOLVE__
#include "assembly_graph.h"

struct fixed_size_q_t{
	int *q;
	int front;
	int back;
};

void init_queue(struct fixed_size_q_t *fq, int size);
void push_queue(struct fixed_size_q_t *fq, int v);
int get_queue(struct fixed_size_q_t *fq);
void pop_queue(struct fixed_size_q_t *fq);
int is_queue_empty(struct fixed_size_q_t *fq);
void destroy_queue(struct fixed_size_q_t *fq);

struct vertex_t{
	int deg_in;
	int deg_out;
	int *child;
	int *parent;
};

struct virtual_graph_t {
	int n_v;
	struct vertex_t *vertices;
	int *exist;
};

struct vertex_height_t{
	int v;
	int height;
};

struct resolve_bulges_bundle_t{
	struct virtual_graph_t *graph;
	int source;
	int *dom;
	int *B;
	int *S;
	int *T;
	int *g;
	int *j;
	int *height;
};

void asm_graph_to_virtual(struct asm_graph_t *g, struct virtual_graph_t *vg);
void clone_virtual_graph(struct virtual_graph_t *org, struct virtual_graph_t *clone);
void virtual_graph_destroy(struct virtual_graph_t *vg);
void init_resolve_bulges(struct asm_graph_t *g, struct resolve_bulges_bundle_t *bundle);
void reset_source(struct resolve_bulges_bundle_t *bundle, int s);
//void bfs_by_vertice(struct asm_graph_t *g, int v, int **path_len);
void get_dominated_vertices(struct resolve_bulges_bundle_t *bundle);
int get_closure(struct resolve_bulges_bundle_t *bundle);

void get_height_dfs(struct resolve_bulges_bundle_t *bundle, int v);
void get_height(struct resolve_bulges_bundle_t *bundle);
int get_skeleton(struct resolve_bulges_bundle_t *bundle);
void print_dom_debug(struct opt_proc_t *opt, struct asm_graph_t *g);
void print_closure_debug(struct opt_proc_t *opt, struct asm_graph_t *g);
void add_vertex_to_B_dfs(struct resolve_bulges_bundle_t *bundle, int v,
		int *in_queue, struct fixed_size_q_t *interested, int depth);
void resolve_bulges(struct asm_graph_t *g);
#endif
