#ifndef SCAFFOLDING_GLOBAL_PARAMS_H
#define SCAFFOLDING_GLOBAL_PARAMS_H

#include "assembly_graph.h"

#define LIST_GLOBAL_PARAMS \
	X(float, global_thres_bucks_score, -1)\
	X(int, global_thres_count_kmer , -1)\
	X(int, global_thres_length , -1)\
	X(int, global_thres_length_min , -1)\
	X(int, global_thres_n_buck_big_small , -1)\
	X(int, global_n_buck , -1)\
	X(float, global_genome_coverage, -1)\
	X(int, global_molecule_length, -1)\
	X(float, global_avg_sum_bin_hash, -1)\
	X(float, global_thres_coefficent, -1); 

#define X(type, name, default_value) extern type name;
LIST_GLOBAL_PARAMS
#undef X

extern int log_level;

void init_global_params(struct asm_graph_t *g);
void check_global_params();

#endif
