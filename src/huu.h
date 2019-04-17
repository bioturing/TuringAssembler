#ifndef __HUU_H__
#define __HUU_H__
void list_contig(struct asm_graph_t *g, FILE *out_file);
void build_graph_2(FILE *fp, FILE *out_file, struct asm_graph_t *g);
void check_contig(struct asm_graph_t *g);
#endif  /* __ASSEMBLY_GRAPH_H__ */
