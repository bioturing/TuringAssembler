#include <stdlib.h>
#include <stdio.h>
#include "scaffolding/edge.h"
#include "assembly_graph.h"
#include "verbose.h"
#include "global_params.h"

void normalize_min_index(struct asm_graph_t *g, struct contig_edge *e)
{
	assert(e->src < g->n_e && e->des < g->n_e && e->src>=0 && e->des>=0);
	int rc_id_src = g->edges[e->src].rc_id;
	if (rc_id_src < e->src) {
		e->src = rc_id_src;
		e->rv_src ^= 1;
	}
	int rc_id_des = g->edges[e->des].rc_id;
	if (rc_id_des < e->des) {
		e->des = rc_id_des ;
		e->rv_des ^= 1;
	}
}

void normalize_one_dir(struct asm_graph_t *g, struct contig_edge *e)
{
	assert(e->src < g->n_e && e->des < g->n_e && e->src >=0 && e->des >=0);
	if (e->rv_src == 1) {
//		e->src = g->edges[e->src].rc_id;
		e->src ^= 1;
		e->rv_src = 0;
	}
	if (e->rv_des == 1) {
//		e->des = g->edges[e->des].rc_id;
		e->des ^= 1;
		e->rv_des = 0;
	}
}

