#ifndef __SCAFFOLDING_H__
#define __SCAFFOLDING_H__
void scaffolding_test(struct asm_graph_t *g, struct opt_proc_t *opt);
void scaffolding(FILE *out_file, struct asm_graph_t *g,
		struct opt_proc_t *opt);
void test_sort_read(struct read_path_t *read_sorted_path, struct asm_graph_t *g);
#endif  /* __ASSEMBLY_GRAPH_H__ */
