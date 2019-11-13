#include "scaffolding/global_params.h"
#include <assert.h>
#include "assembly_graph.h"
#include <stdlib.h>
#include "verbose.h"
#include "scaffolding/bin_hash.h"
#include "utils.h"
#include "log.h"

#define X(type, name, default_value) type name=default_value;
LIST_GLOBAL_PARAMS
#undef X

void check_global_params()
{
#define X(type, name, default_value) assert((name) != default_value);
	LIST_GLOBAL_PARAMS
#undef X
}

void init_global_params(struct asm_graph_t *g)
{
	log_info("------Init global params ---------");
	global_thres_length = 4000;
	global_thres_short_len = 100;
	global_molecule_length = 20000;
	global_avg_sum_bin_hash = get_avg_barcode(g);
	global_thres_coefficent = 0.20;
	global_genome_coverage = get_genome_coverage_h(g);
	log_info("Global genome coverrage %f", global_genome_coverage);
	global_filter_constant = 30;
	global_n_candidate = 11;
	global_count_bc_size = 3000;
	global_distance = 60000;
	global_number_n = 100;
	log_info("------Init global params done---------");
	check_global_params(g);
}

