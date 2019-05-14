#ifndef SCAFFOLDING_OUTPUT_H
#define SCAFFOLDING_OUTPUT_H
#include "assembly_graph.h"
void print_contig(struct asm_graph_t *g, FILE *out_file, int index, int n_contig, int *list_contig);
gint_t dump_edge_seq_reduce_N(char **seq, uint32_t *m_seq, struct asm_edge_t *e);
void print_seq(FILE *fp, int index, char *seq, int len, int cov);
void print_gfa_from_E(struct asm_graph_t *g, int n_e, struct contig_edge *listE, int n_v, int *listV, FILE *out_graph);
#endif
