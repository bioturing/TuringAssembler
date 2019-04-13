#ifndef __RESOLVE_H__
#define __RESOLVE_H__

#include "assembly_graph.h"

int check_simple_loop(struct asm_graph_t *g, gint_t e, double uni_cov);

/* Condense graph by merging non-branching paths */
void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void remove_tips(struct asm_graph_t *g0, struct asm_graph_t *g);
void remove_tips_topology(struct asm_graph_t *g0, struct asm_graph_t *g);
void remove_bubble_and_loop(struct asm_graph_t *g0, struct asm_graph_t *g);
void resolve_bridge(struct asm_graph_t *g);
void detect_simple_tandem(struct asm_graph_t *g0);
// void graph_expanding(struct asm_graph_t *g);
void detect_complex(struct asm_graph_t *g, uint32_t min_contig_size, uint32_t max_edge_count);
void resolve_chain(struct asm_graph_t *g0, struct asm_graph_t *g1);
void resolve_1_1_jungle(struct asm_graph_t *g0, struct asm_graph_t *g1);
void resolve_n_m_jungle(struct asm_graph_t *g0, struct asm_graph_t *g1);
void resolve_n_m_simple(struct asm_graph_t *g0, struct asm_graph_t *g1);

#endif /* __RESOLVE_H__ */
