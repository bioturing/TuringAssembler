#ifndef __COMPLEX_RESOLVE__
#define __COMPLEX_RESOLVE__
#include "assembly_graph.h"

struct vertex_t{
	int deg_in;
	int deg_out;
	int *children;
	int *parent;
};

struct virtual_graph_t {
	int n_v;
	struct vertex_t *vertices;
	int *exist;
};

void asm_graph_to_virtual(struct asm_graph_t *g, struct virtual_graph_t *vg);
void virtual_graph_destroy(struct virtual_graph_t *vg);
//void bfs_by_vertice(struct asm_graph_t *g, int v, int **path_len);
//void get_dominated_vertices(struct virtual_graph_t *vg, int v, int **dom, int *n_dom);
//void get_closure(struct asm_graph_t *g, int *B_cand, int n_cand, int **B, int *n_B);
//void asm_resolve_complex_bulges_ite(struct opt_proc_t *opt, struct asm_graph_t *g);
#endif
