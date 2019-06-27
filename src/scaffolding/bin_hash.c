#include "assembly_graph.h"
#include "barcode_hash.h"
#include "verbose.h"

#include "scaffolding/global_params.h"
#define K31_NULL		((k31key_t)-1)

int count_barcode(struct asm_graph_t *g, struct barcode_hash_t *buck)
{
	const int thres_cnt = global_thres_count_kmer;
	int cnt = 0;
	
	for (uint32_t i = 0; i < buck->size; ++i) {
		if (buck->keys[i] != (uint64_t)(-1)) {
			cnt++;
		}
	}
	return cnt;
}

float get_avg_barcode(struct asm_graph_t *g)
{
	int count = 0;
	uint64_t sum = 0;
	for (int i = 0; i < g->n_e; ++i) {
		int tmp = count_barcode(g, &g->edges[i].barcodes);
		sum += tmp;
		count++;
	}
	return 1.0*sum/count;
}

