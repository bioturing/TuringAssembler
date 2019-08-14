#include "assembly_graph.h"
#include "barcode_hash.h"
#include "verbose.h"

#include "scaffolding/global_params.h"
#define K31_NULL		((k31key_t)-1)

int count_barcode(struct asm_graph_t *g, struct barcode_hash_t *buck)
{
	return buck->n_item;
}

float get_avg_barcode(struct asm_graph_t *g)
{
	int count = 0;
	uint64_t sum = 0;
	for (int i = 0; i < g->n_e; ++i) {
		if (get_edge_len(&g->edges[i]) > global_count_bc_size) {
			int tmp = count_barcode(g, &g->edges[i].barcodes_scaf);
			VERBOSE_FLAG(0, "len and bc %d %d\n", get_edge_len(&g->edges[i]), tmp);
			sum += tmp;
			count++;
		}
	}
	return 1.0*sum/count;
}

