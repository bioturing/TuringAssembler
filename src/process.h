#ifndef __PROCESS_H__
#define __PROCESS_H__

void graph_convert_process(struct opt_proc_t *opt);
void clean_process(struct opt_proc_t *opt);
void build_0_process(struct opt_proc_t *opt);
void build_0_1_process(struct opt_proc_t *opt);
void build_1_2_process(struct opt_proc_t *opt);
void build_2_3_process(struct opt_proc_t *opt);
void build_3_4_process(struct opt_proc_t *opt);
void build_3_4_no_bc_rebuild_process(struct opt_proc_t *opt);
void build_4_5_process(struct opt_proc_t *opt);
void build_scaffolding_1_2_process(struct opt_proc_t *opt);
void build_scaffolding_test_process(struct opt_proc_t *opt);
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
#endif /* __PROCESS_H__ */
