#ifndef __RESOLVE_H__
#define __RESOLVE_H__

#include "assembly_graph.h"

/* Condense graph by merging non-branching paths */
void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void remove_tips(struct asm_graph_t *g0, struct asm_graph_t *g);
void remove_tips_topology(struct asm_graph_t *g0, struct asm_graph_t *g);
void remove_bubble_and_loop(struct asm_graph_t *g0, struct asm_graph_t *g);
void resolve_bridge(struct asm_graph_t *g);
void detect_simple_tandem(struct asm_graph_t *g0);
void graph_expanding(struct asm_graph_t *g);

#endif /* __RESOLVE_H__ */
