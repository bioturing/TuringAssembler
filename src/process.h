#ifndef __PROCESS_H__
#define __PROCESS_H__
#include "minimizers/minimizers.h"
#include "khash.h"
#include "assembly_graph.h"
#include "sort_read.h"
#include "attribute.h"
#include "barcode_builder.h"

void graph_convert_process(struct opt_proc_t *opt);
void clean_process(struct opt_proc_t *opt);
void build_0_process(struct opt_proc_t *opt);
void build_0_1_process(struct opt_proc_t *opt);
void build_0_1(struct opt_proc_t *opt, struct asm_graph_t *g0, struct asm_graph_t *g);
void build_local_0_1(struct asm_graph_t *g0, struct asm_graph_t *g);
void build_1_2_process(struct opt_proc_t *opt);
void build_2_3_process(struct opt_proc_t *opt);
void build_3_4_process(struct opt_proc_t *opt);
void build_3_4_no_bc_rebuild_process(struct opt_proc_t *opt);
void build_4_5_process(struct opt_proc_t *opt);
void build_scaffolding_1_2_process(struct opt_proc_t *opt);
void dirty_process(struct opt_proc_t *opt);
void build_barcode_process(struct opt_proc_t *opt);
void build_barcode_info(struct opt_proc_t *opt);
void build_barcode_scaffold(struct opt_proc_t *opt);
void build_barcode_process_fasta(struct opt_proc_t *opt);
void build_barcode_process_fastg(struct opt_proc_t *opt);
void graph_query_process(struct opt_proc_t *opt);
void assembly_process(struct opt_proc_t *opt);
void assembly2_process(struct opt_proc_t *opt);
void assembly3_process(struct opt_proc_t *opt);
void assembly_precount_process(struct opt_proc_t *opt);
void build_bridge_process(struct opt_proc_t *opt);
void resolve_local_process(struct opt_proc_t *opt);
void save_graph_info(const char *out_dir, struct asm_graph_t *g, const char *suffix);
void reduce_read_process(struct opt_proc_t *opt);
void resolve_bulges_process(struct opt_proc_t *opt);
void resolve_complex_bulges_process(struct opt_proc_t *opt);
void index_mm_process(struct opt_proc_t *opt);
void hits_barcode_process(struct opt_proc_t *opt);
void count_bx_process(struct opt_proc_t *opt);
void sort_read_process(struct opt_proc_t *opt);
void split_molecules_wrapper(struct opt_proc_t *opt);
void split_molecules_process(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct mm_db_edge_t *mm_edges, khash_t(bcpos) *bx_pos_dict);
void print_barcode_graph_process(struct opt_proc_t *opt);
void cluster_molecules_process(struct opt_proc_t *opt);
void debug_process(struct opt_proc_t *opt);
void build_coverage_process(struct opt_proc_t *opt);
void read_pairs_count_process(struct opt_proc_t *opt);
void build_barcode_coverage_info(struct opt_proc_t *opt);
void resolve_212_cov_process(struct opt_proc_t *opt);
void resolve_molecule_process(struct opt_proc_t *opt);
#endif /* __PROCESS_H__ */
