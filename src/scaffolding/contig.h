#ifndef SCAFFOLD_CONTIG_H
#define SCAFFOLD_CONTIG_H
#include "scaffolding/score.h"
#include "assembly_graph.h"
void normalize_min_index(struct asm_graph_t *g, struct edges_score_type *edges_score);
void normalize_one_dir(struct asm_graph_t *g, struct edges_score_type *edges_score);
int is_long_contig(struct asm_edge_t *e);
int is_short_contig(struct asm_edge_t *e);
int is_very_short_contig(struct asm_edge_t *e);
#endif

