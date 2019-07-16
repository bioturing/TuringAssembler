#ifndef __RESOLVE_H__
#define __RESOLVE_H__

#include "assembly_graph.h"

int check_simple_loop(struct asm_graph_t *g, gint_t e, double uni_cov);

void resolve_graph_operation(struct asm_graph_t *g0, struct asm_graph_t *g);

void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void asm_lazy_condense(struct asm_graph_t *g);
void resolve_chain(struct asm_graph_t *g0, struct asm_graph_t *g1);
void resolve_n_m_simple(struct asm_graph_t *g0, struct asm_graph_t *g1);
void resolve_complex(struct asm_graph_t *g, struct asm_graph_t *gd);

#endif /* __RESOLVE_H__ */
