#ifndef __COMPLEX_RESOLVE__
#define __COMPLEX_RESOLVE__
#include "assembly_graph.h"
#include "simple_queue.h"
#include "khash_operations.h"

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
	khash_t(set_int) *dom;
	khash_t(set_int) *B;
	khash_t(int_int) *PE;
	khash_t(int_int) *L;
};

void init_resolve_bulges(struct asm_graph_t *g, struct resolve_bulges_bundle_t *bundle);
void reset_source(struct resolve_bulges_bundle_t *bundle, int s);
void bulges_bundle_destroy(struct resolve_bulges_bundle_t *bundle);

void get_dominated_vertices(struct resolve_bulges_bundle_t *bundle);
int get_closure(struct resolve_bulges_bundle_t *bundle);

void bfs_to_sinks(struct resolve_bulges_bundle_t *bundle);
void get_distance(struct resolve_bulges_bundle_t *bundle);
int get_next_B_candidate(struct resolve_bulges_bundle_t *bundle);
int is_complex_closure(struct resolve_bulges_bundle_t *bundle);
int is_closure_tree(struct resolve_bulges_bundle_t *bundle);

void add_vertex_to_B(struct resolve_bulges_bundle_t *bundle, int v);
void add_vertex_to_B_dfs(struct resolve_bulges_bundle_t *bundle, int v,
		khash_t(set_int) *in_queue, struct queue_t *q, int depth);
void supress_bulge(struct resolve_bulges_bundle_t *bundle);
int resolve_bulges(struct asm_graph_t *g);
int asm_resolve_complex_bulges_ite(struct asm_graph_t *g);
int get_adj_node(struct asm_graph_t *g, int v, int id);
#endif
