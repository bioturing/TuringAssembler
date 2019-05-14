#include "attribute.h"
#include "k31_count.h"
#include "assembly_graph.h"

#include "scaffolding/global_params.h"

float count_unique_bin_hash(struct asm_graph_t *g, struct barcode_hash_t *buck)
{
	extern int global_thres_count_kmer;
	const int thres_cnt = global_thres_count_kmer;
	int cnt = 0;
	
	for (uint32_t i = 0; i < buck->size; ++i) {
		if (buck->cnts[i] != (uint32_t)(-1)) {
			if (buck->cnts[i] >= (uint32_t)thres_cnt)
				cnt++;
		}
	}
	return cnt;
}

float count_sum_bin_hash(struct asm_graph_t *g, struct barcode_hash_t *buck)
{
	const int thres_cnt = global_thres_count_kmer;
	int cnt = 0;
	
	for (uint32_t i = 0; i < buck->size; ++i) {
		if (buck->cnts[i] != (uint32_t)(-1)) {
			cnt += buck->cnts[i];
		}
	}
	return cnt;
}

float get_avg_unique_bin_hash(struct asm_graph_t *g)
{
	int count = 0;
	long long sum = 0;
	for (int i = 0; i < g->n_e; ++i) {
		int n_bucks = (get_edge_len(&g->edges[i]) + g->bin_size-1) / g->bin_size;
		for (int j = 0; j < n_bucks - 1; j++){
			int tmp = count_unique_bin_hash(g, &g->edges[i].bucks[j]);
			if (tmp > 0) {
				count++;
				sum += tmp;
			}
		}
	}
	return 1.0*sum/count;
}

float get_avg_sum_bin_hash(struct asm_graph_t *g)
{
	int count = 0;
	long long sum = 0;
	for (int i = 0; i < g->n_e; ++i) {
		int n_bucks = (get_edge_len(&g->edges[i]) + g->bin_size-1) / g->bin_size;
		for (int j = 0; j < n_bucks - 1; j++){
			int tmp = count_sum_bin_hash(g, &g->edges[i].bucks[j]);
			if (tmp > 0) {
				count++;
				sum += tmp;
			}
		}
	}
	return 1.0*sum/count;
}

