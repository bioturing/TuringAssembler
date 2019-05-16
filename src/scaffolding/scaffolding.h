#ifndef __SCAFFOLDING_H__
#define __SCAFFOLDING_H__
void build_list_edges(struct asm_graph_t *g, FILE *out_file, struct opt_proc_t *opt);
void connect_contig(FILE *fp, FILE *out_file, FILE *out_graph, struct asm_graph_t *g);
void scaffolding_test(struct asm_graph_t *g, struct opt_proc_t *opt);
#endif  /* __ASSEMBLY_GRAPH_H__ */
