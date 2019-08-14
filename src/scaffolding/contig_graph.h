#ifndef SCAFFOLDING_CONTIG_GRAPH_H
#define SCAFFOLDING_CONTIG_GRAPH_H
#include "assembly_graph.h"
struct contig_graph {
	int *head, *next, *edges, n_v, n_e, *edge_index;
};

struct contig_graph* build_graph_contig(struct asm_graph_t *g, int *mark_short);
void find_path_short(int src, int des, struct asm_graph_t *g, struct contig_graph *contig_g,
		int *n_insert_path, int **insert_path, int *mark_short, struct barcode_hash_t *hl,
		struct barcode_hash_t *hr);

#endif
