#ifndef __RESOLVE_H__
#define __RESOLVE_H__

#include "assembly_graph.h"
#include "barcode_resolve2.h"
void resolve_graph_operation(struct asm_graph_t *g0, struct asm_graph_t *g);
void resolve_local_graph_operation(struct asm_graph_t *g0,
		struct asm_graph_t *g);
void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void asm_condense_barcode(struct asm_graph_t *g0, struct asm_graph_t *g);
void asm_lazy_condense(struct asm_graph_t *g);
void resolve_chain(struct asm_graph_t *g0, struct asm_graph_t *g1);
void resolve_n_m_simple(struct asm_graph_t *g0, struct asm_graph_t *g1);
void resolve_complex(struct asm_graph_t *g, struct asm_graph_t *gd);
void resolve_n_m_local(struct opt_proc_t *opt, struct read_path_t *rpath,
				struct asm_graph_t *g0, struct asm_graph_t *g1);
void do_something_local(struct opt_proc_t *opt, struct asm_graph_t *g);
void do_some_resolve_bridge(struct asm_graph_t *g);
int resolve_loop(struct asm_graph_t *g0);
int asm_resolve_dump_loop(struct asm_graph_t *g);
int asm_resolve_dump_branch(struct asm_graph_t *g);
int detect_dump_jungle(struct asm_graph_t *g, int e, int **dump_edges, int *n_dump);
int asm_resolve_dump_jungle(struct opt_proc_t *opt, struct asm_graph_t *g);
void get_shared_barcode_reads(struct opt_proc_t *opt, struct asm_graph_t *g,
		int e1, int e2, struct read_path_t *local_read_path);
int asm_resolve_dump_loop_ite(struct asm_graph_t *g);
int asm_resolve_dump_jungle_ite(struct opt_proc_t *opt, struct asm_graph_t *g);
void get_union_barcode_reads(struct opt_proc_t *opt, struct asm_graph_t *g,
		int e1, int e2, khash_t(bcpos) *dict, struct read_path_t *read_sorted_path,
		struct read_path_t *local_read_path);
void destroy_read_path(struct read_path_t *reads);
#endif /* __RESOLVE_H__ */
