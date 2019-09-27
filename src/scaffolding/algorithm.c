#include "assembly_graph.h"

int get_rc_id(struct asm_graph_t *g, int i_e)
{
	return g->edges[i_e].rc_id;
}

