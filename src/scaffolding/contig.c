#include <stdlib.h>
#include <stdio.h>
#include "scaffolding/edge.h"
#include "scaffolding/score.h"
#include "assembly_graph.h"
#include "verbose.h"
#include "global_params.h"
#include <assert.h>

//void normalize_min_index(struct asm_graph_t *g, struct edges_score_type *edges_score)
//{
//
//	for (int i = 0; i < edges_score->n_edge; i++) {
//		struct scaffold_edge *e = &edges_score->list_edge[i];
//		assert(e->src < g->n_e && e->des < g->n_e && e->src>=0 && e->des>=0);
//		int rc_id_src = g->edges[e->src].rc_id;
//		if (rc_id_src < e->src) {
//			e->src = rc_id_src;
//			e->rv_src ^= 1;
//		}
//		int rc_id_des = g->edges[e->des].rc_id;
//		if (rc_id_des < e->des) {
//			e->des = rc_id_des ;
//			e->rv_des ^= 1;
//		}
//	}
//}
//
//void normalize_one_dir(struct asm_graph_t *g, struct edges_score_type *edges_score)
//{
//	for (int i = 0; i < edges_score->n_edge; i++) {
//		struct scaffold_edge *e = &edges_score->list_edge[i];
//		assert(e->src < g->n_e && e->des < g->n_e && e->src >=0 && e->des >=0);
//		if (e->rv_src == 1) {
//			e->src = g->edges[e->src].rc_id;
//			e->rv_src = 0;
//		}
//		if (e->rv_des == 1) {
//			e->des = g->edges[e->des].rc_id;
//			e->rv_des = 0;
//		}
//	}
//}

int is_long_contig(struct asm_edge_t *e)
{
	int len = e->seq_len;
	if (len >= global_thres_length){
		return 1;
	}
	return 0;
}

int is_short_contig(struct asm_edge_t *e)
{
	int len = e->seq_len;
	if (global_thres_length > len && len >= global_thres_short_len ){
		assert(len > 0);
		return 1;
	}
	return 0;
}

int is_very_short_contig(struct asm_edge_t *e)
{
	int len = e->seq_len;
	if (global_thres_short_len > len){
		return 1;
	}
	return 0;
}

