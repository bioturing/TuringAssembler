#ifndef __HUU_H__
#define __HUU_H__
void build_list_contig(struct asm_graph_t *g, FILE *out_file, float huu_1_score, struct opt_proc_t *opt);
void connect_contig(FILE *fp, FILE *out_file, FILE *out_graph, struct asm_graph_t *g, float huu_1_score);
void check_contig(struct asm_graph_t *g, float huu_1_score);
#endif  /* __ASSEMBLY_GRAPH_H__ */
