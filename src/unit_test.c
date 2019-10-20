#include "unit_test.h"
#include "log.h"
#include "helper.h"
#include "io_utils.h"
#include "build_bridge.h"

void test_original_reads_sorted(struct opt_proc_t *opt)
{
	if (opt->lib_type != LIB_TYPE_SORTED)
		log_error("Test failed, reads must be sorted");
	// TODO: Check if reads are actually sorted
}

void precon_get_all_local_graphs(struct opt_proc_t *opt)
{
	test_original_reads_sorted(opt);
}

void precon_get_local_assembly(struct opt_proc_t *opt)
{
	test_original_reads_sorted(opt);
}

void postcon_get_local_assembly(int lksize, char *work_dir, int is_reads_exist)
{
	char graph_path[MAX_PATH];
	sprintf(graph_path, "%s/graph_k_%d_local_lvl_1.bin", work_dir, lksize);
	int graph_exist = check_file_exist(graph_path);
	if (is_reads_exist && !graph_exist)
		log_error("Read files exist but no graph found");
	if (!is_reads_exist && graph_exist)
		log_error("Read files are empty but still found graph, please clean the folder and try again");
}

void test_bridge_result(char *bridge_seq, int seq_len, int bridge_type)
{
	if (seq_len == 0)
		log_error("The result bridge is empty");
	if ((int) strlen(bridge_seq) != seq_len)
		log_error("The length of the sequence is not correct: %d vs %d",
				strlen(bridge_seq), seq_len);
	if (bridge_type != BRIDGE_LOCAL_NOT_FOUND && bridge_type != BRIDGE_TRIVIAL_BRIDGE
		&& bridge_type != BRIDGE_PATH_NOT_FOUND && bridge_type != BRIDGE_MULTIPLE_PATH)
		log_error("Bridge type is not in the expected list");
}

void postcon_get_bridge(char *bridge_seq, int seq_len, int bridge_type)
{
	test_bridge_result(bridge_seq, seq_len, bridge_type);
}

void precon_get_local_edge_head(struct asm_graph_t *g, int e)
{
	precon_get_local_edge(g, e);
}

void postcon_get_local_edge_head(struct asm_graph_t *g, struct asm_graph_t *lg,
		int e, struct edge_map_info_t *emap)
{
	postcon_get_local_edge(g, lg, e, emap);
}

void precon_get_local_edge_tail(struct asm_graph_t *g, int e)
{
	precon_get_local_edge(g, e);
}

void postcon_get_local_edge_tail(struct asm_graph_t *g, struct asm_graph_t *lg,
		int e, struct edge_map_info_t *emap)
{
	postcon_get_local_edge(g, lg, e, emap);
}

void precon_get_local_edge(struct asm_graph_t *g, int e)
{
	test_edge_in_graph(e, g);
}

void postcon_get_local_edge(struct asm_graph_t *g, struct asm_graph_t *lg,
		int e, struct edge_map_info_t *emap)
{
	if (e != emap->gl_e)
		log_error("e (%d) != emap->gl_e (%d)", e, emap->gl_e);
	test_edge_map(emap, g, lg);
	// TODO: test if the mapping is actually correct 
}

void test_mapping_range(struct subseq_pos_t *ss_pos, struct asm_edge_t *e)
{
	int ok = (0 <= ss_pos->start) && (ss_pos->start <= ss_pos->end)
		&& (ss_pos->end < (int) e->seq_len);
	if (!ok)
		log_error("Edge mapping error");
}

void test_edge_in_graph(int e, struct asm_graph_t *g)
{
	if (e < 0 || e >= (int) g->n_e)
		log_error("Graph has %d edges but the queried edge id is %d",
				g->n_e, e);
}

void test_edge_map(struct edge_map_info_t *emap, struct asm_graph_t *g,
		struct asm_graph_t *lg)
{
	test_edge_in_graph(emap->gl_e, g);
	if (emap->lc_e != -1){
		test_edge_in_graph(emap->lc_e, lg);
		test_mapping_range(&(emap->gpos), g->edges + emap->gl_e);
		test_mapping_range(&(emap->lpos), lg->edges + emap->lc_e);
	}
}

void precon_try_bridging(struct asm_graph_t *g, struct asm_graph_t *lg,
		int *scaffolds, int n_scaff, struct edge_map_info_t *emap1,
		struct edge_map_info_t *emap2)
{
	for (int i = 0; i < n_scaff; ++i){
		int e = scaffolds[i];
		test_edge_in_graph(e, g);
	}
	test_edge_map(emap1, g, lg);
	test_edge_map(emap2, g, lg);
}

void postcon_try_bridging(char *bridge_seq, int seq_len, int bridge_type)
{
	test_bridge_result(bridge_seq, seq_len, bridge_type);
}
