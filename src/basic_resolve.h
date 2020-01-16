#ifndef BASIC_RESOLVE_H
#define BASIC_RESOLVE_H
void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void asm_condense_h(struct asm_graph_t *g);
void remove_duplicate_edge(struct asm_graph_t *g);
void resolve_graph_small_operation(struct asm_graph_t *g0, struct asm_graph_t *g);
#endif
