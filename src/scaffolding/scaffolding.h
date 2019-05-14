#ifndef __SCAFFOLDING_H__
#define __SCAFFOLDING_H__
void build_list_contig(struct asm_graph_t *g, FILE *out_file, struct opt_build_t *opt);
void connect_contig(FILE *fp, FILE *out_file, FILE *out_graph, struct asm_graph_t *g);
void check_contig(struct asm_graph_t *g);
#endif  /* __ASSEMBLY_GRAPH_H__ */
