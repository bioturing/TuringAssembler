#ifndef __UNIT_TEST__
#define __UNIT_TEST__
#include "attribute.h"
#define pre_test(test_function, ...) do{\
	log_info("Testing pre-conditions for " #test_function);\
	precon_##test_function(__VA_ARGS__);\
	log_info("Test passed");\
	} while(0)\

#define post_test(test_function, ...) do{\
	log_info("Testing post-conditions for " #test_function);\
	postcon_##test_function(__VA_ARGS__);\
	log_info("Test passed");\
	} while(0)\

void precon_get_all_local_graphs(struct opt_proc_t *opt);
void precon_get_local_assembly(struct opt_proc_t *opt);
void postcon_get_local_assembly(int lksize, char *work_dir, int is_reads_exist);
void test_original_reads_sorted(struct opt_proc_t *opt);

void postcon_get_bridge(char *bridge_seq, int seq_len, int bridge_type);
#endif
