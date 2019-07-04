#include <stdint.h>
#include <assert.h>
#include "verbose.h"
#include <stdlib.h>
#include <math.h>
#include "assembly_graph.h"
#include "scaffolding/output.h"
#include "scaffolding/global_params.h"
#include "scaffolding/edge.h"
#include "scaffolding/contig_graph.h"

int get_relative_cov(struct asm_graph_t *g, int i_e)
{
	float cvr = global_genome_coverage;
	return lround(__get_edge_cov(&g->edges[i_e], g->ksize)/cvr);
}


int hround(float x)
{
	return lround(x);
}

int get_rc_id(struct asm_graph_t *g, int i_e)
{
	return g->edges[i_e].rc_id;
}

int get_rc_id_V(struct asm_graph_t *g, int n_v, int *listV, int i_e)
{
	int rc = g->edges[listV[i_e]].rc_id;
	if ((i_e^1) < n_v && listV[i_e^1] == rc) {
		return i_e^1;
	} else{
		for (int i = 0; i < n_v; i++) {
			if (listV[i] == rc) {
				VERBOSE_FLAG(0, "rcidV %d %d\n", i_e, i);
				return i;
			}
		}
	}
}

int get_new_size(int size)
{
	if (size == 0)
		return 1;
	return size*2;
}

float get_component_cov(struct asm_graph_t *g, int n_component_node, int *connected_component, int *listV)
{
	uint64_t sum_count = 0, sum_len = 0, sum_n_holes = 0;
	for (int i = 0 ; i < n_component_node; i++) {
		struct asm_edge_t *e = &g->edges[listV[connected_component[i]]];
		sum_count += e->count;
		sum_len += e->seq_len;
		sum_n_holes += e->n_holes;
	}
	return 1.0*sum_count/(sum_len - g->ksize * sum_n_holes - g->ksize * n_component_node);
}


