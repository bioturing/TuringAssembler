#ifndef __CLUSTER_MOLECULES__
#define __CLUSTER_MOLECULES__
#include "assembly_graph.h"
#include "complex_resolve.h"
#include "minimizers/minimizers.h"

struct dijkstra_node_t{
	int vertex;
	int len;
};

int get_shortest_path(struct asm_graph_t *g, int source, int target);
void get_sub_graph_of_molecules(struct asm_graph_t *g, struct mm_hits_t *hits);
#endif
