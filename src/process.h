#ifndef __PROCESS_H__
#define __PROCESS_H__

void graph_convert_process(struct opt_proc_t *opt);
void clean_process(struct opt_proc_t *opt);
void build_0_process(struct opt_proc_t *opt);
void build_0_1_process(struct opt_proc_t *opt);
void build_1_2_process(struct opt_proc_t *opt);
void build_2_3_process(struct opt_proc_t *opt);
void build_3_4_process(struct opt_proc_t *opt);
void build_barcode_process(struct opt_proc_t *opt);
void graph_query_process(struct opt_proc_t *opt);
void assembly_process(struct opt_proc_t *opt);
void assembly2_process(struct opt_proc_t *opt);


#endif /* __PROCESS_H__ */
