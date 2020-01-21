#include "assembly_graph.h"
#include "barcode_hash.h"

#include "scaffolding/global_params.h"

float get_avg_barcode(struct asm_graph_t *g)
{
	int count = 0;
	uint64_t sum = 0;
	for (int i = 0; i < g->n_e; ++i) {
		if (get_edge_len(&g->edges[i]) > (uint32_t) global_count_bc_size) {
			sum += g->edges[i].barcodes_scaf->n_item;
			count++;
		}
	}
	return 1.0*sum/count;
}

