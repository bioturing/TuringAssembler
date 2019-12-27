#ifndef BASIC_RESOLVE_H
#define BASIC_RESOLVE_H
void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void asm_condense_h(struct asm_graph_t *g);

int check_junction(struct asm_graph_t *g, int v);
int asm_resolve_1_2_junctions_ite(struct asm_graph_t *g);
int asm_resolve_1_2_junctions(struct asm_graph_t *g);
#endif
