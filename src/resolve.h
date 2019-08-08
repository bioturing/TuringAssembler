#ifndef __RESOLVE_H__
#define __RESOLVE_H__

#include "assembly_graph.h"

void resolve_graph_operation(struct asm_graph_t *g0, struct asm_graph_t *g);

void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void asm_lazy_condense(struct asm_graph_t *g);
void resolve_chain(struct asm_graph_t *g0, struct asm_graph_t *g1);
void resolve_n_m_simple(struct asm_graph_t *g0, struct asm_graph_t *g1);
void resolve_complex(struct asm_graph_t *g, struct asm_graph_t *gd);
void resolve_n_m_local(struct opt_proc_t *opt, struct read_path_t *rpath,
				struct asm_graph_t *g0, struct asm_graph_t *g1);
void resolve_loop(struct asm_graph_t *g0, struct asm_graph_t *g);

#endif /* __RESOLVE_H__ */
