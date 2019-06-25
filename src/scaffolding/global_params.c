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

float get_global_thres_score(struct asm_graph_t *g)
{
	float cvr = get_genome_coverage(g); float res = 1.0 * g->bin_size / global_molecule_length * global_thres_coefficent;
	VERBOSE_FLAG(1, "global thres score: %f\n", res);
	return res;
}

int get_global_count_kmer(struct asm_graph_t *g)
{
	int *arr_count = NULL, n_arr = 0, *count_count = NULL; 
	 int res = -1;
	for (int i = 0; i < g->n_e; i++) {
			struct barcode_hash_t buck = g->edges[i].barcodes;
			for (uint32_t l = 0; l < buck.n_item; l++){
				arr_count = realloc(arr_count, (n_arr + 1) * sizeof(int));
				arr_count[n_arr] = buck.cnts[l];
				n_arr++;
			}
	}
	int size_count = 1000;
	count_count = realloc(count_count, size_count*sizeof(int)); 
	for(int i = 0; i < n_arr; i++) {
		assert(arr_count[i] >  size_count);
		count_count[arr_count[i]]++;
	}
	int max_count = 0;
	for (int i = 1; i < 10 ; i++){
		max_count = MAX(max_count, count_count[i]);
	}
	VERBOSE_FLAG(1, "max_count %d \n", max_count);
	for(int i = 10; i < size_count; i++) {
		VERBOSE_FLAG(2, "%count hash ",  count_count[i]);
		if (count_count[i] > 3 * max_count){
			res = i;
			break;
		}
	}
	VERBOSE_FLAG(1, "global thres count kmer: %d\n", res);
	assert(res == -1);
	return res;
}

void init_global_params(struct asm_graph_t *g)
{
	global_thres_length = 10000;
	global_thres_short_len= 10000;
	global_thres_n_buck_big_small = 5;
	global_n_buck = 6;
	global_molecule_length = 20000;
	global_thres_count_kmer =  26;//get_global_count_kmer(g);
	global_avg_sum_bin_hash = get_avg_sum_bin_hash(g);
	global_thres_coefficent = 0.20;
	global_genome_coverage = get_genome_coverage(g);
	global_thres_bucks_score = get_global_thres_score(g);
	global_filter_constant = 30;
	global_number_degree = 15;
}

