#ifndef BASIC_RESOLVE_H
#define BASIC_RESOLVE_H
#define MAX_JUNCTION_LEN 2000
void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void asm_condense_h(struct asm_graph_t *g);

int check_junction_cov(struct asm_graph_t *g, int e0, int e1, int e2);
void get_junction_edges(struct asm_graph_t *g, int v, int *e0, int *e1, int *e2);
int check_junction(struct asm_graph_t *g, int v);
int asm_resolve_1_2_junctions_ite(struct asm_graph_t *g);
int asm_resolve_1_2_junctions(struct asm_graph_t *g);
#endif
