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

void postcon_get_bridge(char *bridge_seq, int seq_len, int bridge_type)
{
	if (seq_len == 0)
		log_error("get_bridge function does not return a bridge");
	if (strlen(bridge_seq) != seq_len)
		log_error("The length of the sequence is not correct: %d vs %d",
				strlen(bridge_seq), seq_len);
	if (bridge_type != BRIDGE_LOCAL_NOT_FOUND && bridge_type != BRIDGE_TRIVIAL_BRIDGE
		&& bridge_type != BRIDGE_LOCAL_NOT_FOUND && bridge_type != BRIDGE_MULTIPLE_PATH)
		log_error("Bridge type is not in the expected list");
	// TODO: test if the seq contains any N character, then the number of them
	// must be either 0 or 100
}
