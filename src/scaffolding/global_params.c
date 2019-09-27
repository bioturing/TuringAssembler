#include "scaffolding/global_params.h"
#include <assert.h>
#include "assembly_graph.h"
#include <stdlib.h>
#include "verbose.h"
#include "scaffolding/bin_hash.h"
#include "utils.h"

#define X(type, name, default_value) type name=default_value;
LIST_GLOBAL_PARAMS
#undef X

int log_level = 1;

void check_global_params()
{
#define X(type, name, default_value) assert((name) != default_value);
LIST_GLOBAL_PARAMS
#undef X
}

void init_global_params(struct asm_graph_t *g)
{
    VERBOSE_FLAG(1, "------Init global params ---------\n");
	global_thres_length = 4000;
	global_thres_short_len= 200;
	global_molecule_length = 20000;
	global_avg_sum_bin_hash = get_avg_barcode(g);
	global_thres_coefficent = 0.20;
	global_genome_coverage = get_genome_coverage_h(g);
	VERBOSE_FLAG(1, "Global genome coverrage %f\n", global_genome_coverage);
	global_filter_constant = 30;
	global_n_candidate = 11;
	global_count_bc_size = 3000;
	global_distance = 60000;
    VERBOSE_FLAG(1, "------Init global params done---------\n");
    check_global_params(g);
}

