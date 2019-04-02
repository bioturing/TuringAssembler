#include <stdlib.h>
#include "assembly_graph.h"
#include "k31hash.h"
#include "pthread.h"
#include "verbose.h"

uint32_t min(uint32_t a, uint32_t b)
{
	if (a < b) return a;
	else return b;
}

float getScoreBucks(struct barcode_hash_t *buck0,struct barcode_hash_t *buck1) 
{
	const uint32_t thres_cnt = 30;
	uint32_t cnt0 = 0, cnt1 = 0, res2 = 0;

	for (uint32_t i = 0; i < buck0->size; ++i) {
		if (buck0->cnts[i] != (uint32_t)(-1) && buck0->cnts[i] >= thres_cnt) {
			cnt0++;
		}
	}
	for (uint32_t i = 0; i < buck1->size; ++i) {
		if (buck1->cnts[i] != (uint32_t)(-1) && buck1->cnts[i] >= thres_cnt) {
			cnt1++;
		}
	}

	for (uint32_t i = 0; i < buck0->size; ++i) {
		if ((buck0->keys[i]) != (uint64_t)(-1) && buck0->cnts[i] >= thres_cnt) {
			uint32_t tmp = barcode_hash_get(buck1, buck0->keys[i]);
			if (tmp != BARCODE_HASH_END(buck1) && buck1->cnts[tmp] >= thres_cnt) {
				res2++;
			}
		}
	}
	return 1.0 * res2 / (cnt0 + cnt1);
}

uint32_t abssub(uint32_t a, uint32_t b) {
	if (a>b) 
		return a-b;
	else 
		return b-a;
}

int32_t getScore(struct asm_edge_t *e0, struct asm_edge_t *e1, struct asm_graph_t *g) {
	uint32_t n_bucks = 5;
	uint32_t n0_bucks = (get_edge_len(e0) + g->bin_size-1) / g->bin_size;
	uint32_t n1_bucks = (get_edge_len(e1) + g->bin_size-1) / g->bin_size;
	const float thres_score = 0.01;
	uint32_t score2=0;
	if (n0_bucks < n_bucks || n1_bucks < n_bucks) 
		return -1;
	float res = 0;
	float maxtmp =0 ;
	for (uint32_t i = 0; i < n_bucks; ++i) {
		for (uint32_t j = 0; j < n_bucks; ++j) {
			float tmp = getScoreBucks(&e0->bucks[i], &e1->bucks[i]) / (n_bucks*n_bucks);
			__VERBOSE("%f ", tmp);
			if (tmp > maxtmp) {
				maxtmp = tmp;
				score2 = i+j; 
			}
			res += tmp;
		}
		__VERBOSE("\n");
	}
	__VERBOSE("\n");
	if (res < thres_score) 
		return -1;
	return score2;
}

void listContig(struct asm_graph_t *g) {
	const uint32_t thres_len_e = 10000; 
	uint32_t *listE = NULL;
	uint32_t n_e=0;
	for (uint32_t e = 0; e < g->n_e; ++e) {
		uint32_t len = get_edge_len(&g->edges[e]);
		if (len > thres_len_e) {
			++n_e;
			listE = realloc(listE, n_e*sizeof(uint32_t));
			listE[n_e-1] = e; 
		}
	}
	for (uint32_t i = 0; i < n_e; i++) {
		uint32_t e = listE[i];
		float max_score = 0;
		uint32_t rr = -1;
		for (uint32_t i1 = i+1; i1 < n_e; i1++) {
			uint32_t e1 = listE[i1];
			int32_t score = getScore(&g->edges[e], &g->edges[e1], g);
//			if (score > max_score && e != g->edges[e1].rc_id){
//				max_score = score;
//				rr = g->edges[e1].rc_id;
//			}
			if (score != -1) if (e != g->edges[e1].rc_id) {
				__VERBOSE("edge: %d %ld %d\n", score, g->edges[e1].rc_id, e); 
			}
		}
//		__VERBOSE("edge: %f %d %ld\n", max_score, rr, e); 

	}
}

