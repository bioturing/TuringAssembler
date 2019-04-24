#ifndef __HUU_H__
#define __HUU_H__
void build_list_contig(struct asm_graph_t *g, FILE *out_file);
void connect_contig(FILE *fp, FILE *out_file, FILE *out_graph, struct asm_graph_t *g);
void check_contig(struct asm_graph_t *g);
#endif  /* __ASSEMBLY_GRAPH_H__ */
