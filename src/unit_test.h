#ifndef __UNIT_TEST__
#define __UNIT_TEST__
#include "attribute.h"
#include "assembly_graph.h"
#include "map_contig.h"
#define pre_test(test_function, ...) do{\
	log_debug("Testing pre-conditions for " #test_function);\
	precon_##test_function(__VA_ARGS__);\
	log_debug("Test passed");\
	} while(0)\

#define post_test(test_function, ...) do{\
	log_debug("Testing post-conditions for " #test_function);\
	postcon_##test_function(__VA_ARGS__);\
	log_debug("Test passed");\
	} while(0)\

void precon_get_all_local_graphs(struct opt_proc_t *opt);
void precon_get_local_assembly(struct opt_proc_t *opt);
void postcon_get_local_assembly(int lksize, char *work_dir, int is_reads_exist);
void test_original_reads_sorted(struct opt_proc_t *opt);

void test_bridge_result(char *bridge_seq, int seq_len, int bridge_type);
void postcon_get_bridge(char *bridge_seq, int seq_len, int bridge_type);

void precon_get_local_edge_head(struct asm_graph_t *g, int e);
void postcon_get_local_edge_head(struct asm_graph_t *g, struct asm_graph_t *lg,
		int e, struct edge_map_info_t *emap);

void test_edge_in_graph(int e, struct asm_graph_t *g);
void test_edge_map(struct edge_map_info_t *emap, struct asm_graph_t *g,
		struct asm_graph_t *lg);
void test_mapping_range(struct subseq_pos_t *ss_pos, struct asm_edge_t *e);

void precon_get_local_edge_tail(struct asm_graph_t *g, int e);
void postcon_get_local_edge_tail(struct asm_graph_t *g, struct asm_graph_t *lg,
		int e, struct edge_map_info_t *emap);

void precon_get_local_edge(struct asm_graph_t *g, int e);
void postcon_get_local_edge(struct asm_graph_t *g, struct asm_graph_t *lg,
		int e, struct edge_map_info_t *emap);

void precon_try_bridging(struct asm_graph_t *g, struct asm_graph_t *lg,
		int *scaffolds, int n_scaff, struct edge_map_info_t *emap1,
		struct edge_map_info_t *emap2);
void postcon_try_bridging(char *bridge_seq, int seq_len, int bridge_type);
#endif
