#ifndef __COMPLEX_RESOLVE__
#define __COMPLEX_RESOLVE__
#include "assembly_graph.h"

void *pointerize(void *data, int size);

struct queue_t{
	void **data;
	int cap;
	int front;
	int back;
};

void init_queue(struct queue_t *q, int cap);
void push_queue(struct queue_t *q, void *ptr);
void *get_queue(struct queue_t *q);
void pop_queue(struct queue_t *q);
int is_queue_empty(struct queue_t *q);
void destroy_queue(struct queue_t *q);

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
	int vertex;
	int height;
};

struct resolve_bulges_bundle_t{
	struct asm_graph_t *graph;
	struct queue_t *B_vertices;
	struct queue_t *dom_vertices;
	struct queue_t *closest;
	int source;
	int *dom;
	int *B;
	int *S;
	int *T;
	int *g;
	int *j;
	int *height;
	int *P;
	int *L;
};

void asm_graph_to_virtual(struct asm_graph_t *g, struct virtual_graph_t *vg);
void clone_virtual_graph(struct virtual_graph_t *org, struct virtual_graph_t *clone);
void virtual_graph_destroy(struct virtual_graph_t *vg);

void init_resolve_bulges(struct asm_graph_t *g, struct resolve_bulges_bundle_t *bundle);
void reset_source(struct resolve_bulges_bundle_t *bundle, int s);
void bulges_bundle_destroy(struct resolve_bulges_bundle_t *bundle);

//void bfs_by_vertice(struct asm_graph_t *g, int v, int **path_len);
void get_dominated_vertices(struct resolve_bulges_bundle_t *bundle);
int get_closure(struct resolve_bulges_bundle_t *bundle);


int is_able_to_map(struct resolve_bulges_bundle_t *bundle, int u, int v);
void get_height_dfs(struct resolve_bulges_bundle_t *bundle, int v);
void get_height(struct resolve_bulges_bundle_t *bundle);
int get_skeleton(struct resolve_bulges_bundle_t *bundle);
void bfs_to_sinks(struct resolve_bulges_bundle_t *bundle);
void get_tree(struct resolve_bulges_bundle_t *bundle);
void get_distance(struct resolve_bulges_bundle_t *bundle);
int get_next_B_candidate(struct resolve_bulges_bundle_t *bundle);
int is_complex_closure(struct resolve_bulges_bundle_t *bundle);
int is_closure_tree(struct resolve_bulges_bundle_t *bundle);
void free_queue_content(struct queue_t *q);

void print_dom_debug(struct opt_proc_t *opt, struct asm_graph_t *g);
void print_closure_debug(struct opt_proc_t *opt, struct asm_graph_t *g);
void add_vertex_to_B(struct resolve_bulges_bundle_t *bundle, int v);
void add_vertex_to_B_dfs(struct resolve_bulges_bundle_t *bundle, int v,
		int *in_queue, struct queue_t *q, int depth);
int resolve_bulges(struct asm_graph_t *g);
void remove_edge_virtual(struct asm_graph_t *g, int e,
		struct resolve_bulges_bundle_t *bundle);
int asm_resolve_complex_bulges_ite(struct opt_proc_t *opt, struct asm_graph_t *g);
#endif
