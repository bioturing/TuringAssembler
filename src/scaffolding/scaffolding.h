#ifndef __SCAFFOLDING_H__
#define __SCAFFOLDING_H__

#include <stdio.h>
#include "../utils.h"
#include "../attribute.h"

#define MIN_EDGE_COV_SCAFFOLD 0.15
struct asm_graph_t* create_and_load_graph(struct opt_proc_t *opt);
void dirty(struct asm_graph_t *g, struct opt_proc_t *opt);
void dirty_code(struct opt_proc_t *opt);
void scaffolding(FILE *out_file, struct asm_graph_t *g,
		struct opt_proc_t *opt);
void test_sort_read(struct read_path_t *read_sorted_path, struct asm_graph_t *g);
int concate_edge(struct asm_graph_t *g, int n_path, int *path, int ksize,
		uint32_t **res_seq);
void get_long_contig(struct asm_graph_t *g, struct opt_proc_t *opt);
void resolve_212_by_cov(struct asm_graph_t *g, struct opt_proc_t *opt);
#endif  /* __ASSEMBLY_GRAPH_H__ */
