#ifndef __SCAFFOLDING_H__
#define __SCAFFOLDING_H__

#include <stdio.h>
#include "utils.h"
struct asm_graph_t* create_and_load_graph(struct opt_proc_t *opt);
void dirty_code(struct opt_proc_t *opt);
void scaffolding(FILE *out_file, struct asm_graph_t *g,
		struct opt_proc_t *opt);
void test_sort_read(struct read_path_t *read_sorted_path, struct asm_graph_t *g);
void dirty(struct asm_graph_t *g, struct opt_proc_t *opt);
char *concate_edge(struct asm_graph_t *g, int n_path, int *path, int ksize, int *ret_len,
		   int *total_count, double global_avg_cov, uint32_t **seq);
#endif  /* __ASSEMBLY_GRAPH_H__ */
