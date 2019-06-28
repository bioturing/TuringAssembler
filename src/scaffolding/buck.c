#include "assembly_graph.h"
#include "verbose.h"
#include "scaffolding/global_params.h"
#include "scaffolding/bin_hash.h"
#include "utils.h"

int get_amount_hole(struct asm_graph_t *g, struct asm_edge_t *e)
{
	int res = 0, l = 0, r = MIN_CONTIG_BARCODE, sum_holes = 0;
//	for (int i = 0; i < e->n_holes; ++i){
//		VERBOSE_FLAG(log_hole, "holeee %d %d %d\n" , e->seq_len, e->l_holes[i], e->p_holes[i]);
//	}
	for (uint32_t i = 0; i < e->n_holes; ++i){
		int pos = e->p_holes[i] + sum_holes;
		if (pos >= r){
			break;
		}
		if (pos >= l) {
			res += MIN(r-pos, (int)(e->l_holes[i]));
		}
		sum_holes += e->p_holes[i];
	}
	return res;
}

int check_qualify_buck(struct asm_graph_t *g, struct asm_edge_t *e, float avg_bin_hash, 
		struct opt_proc_t *opt)
{
	if (get_amount_hole(g, e)  > 0.7*MIN_CONTIG_BARCODE) {
		VERBOSE_FLAG(2, "NNNNN size is to big ");
		return 0;
	}
	return 1;
}

float get_share_barcode(struct barcode_hash_t *buck0, struct barcode_hash_t *buck1, float edge0_cov,
		float edge1_cov)
{
	int cnt0 = 0, cnt1 = 0, res2 = 0;
	float ratio0 = edge0_cov / global_genome_coverage, ratio1 = edge1_cov / global_genome_coverage;

	for (uint32_t i = 0; i < buck1->size; ++i) {
		if (buck1->keys[i] != (uint64_t)(-1)) {
			cnt1++;
		}
	}

	for (uint32_t i = 0; i < buck0->size; ++i) {
		if ((buck0->keys[i]) != (uint64_t)(-1)){
			cnt0++;
			uint32_t tmp = barcode_hash_get(buck1, buck0->keys[i]);
			if (tmp != BARCODE_HASH_END(buck1) && buck1->keys[tmp] != (uint64_t)(-1)) {
				res2++;
			}
		}
	}
	VERBOSE_FLAG(3, "res %d cnt0 %d cnt1 %d \n", res2, cnt0, cnt1 );
	if (MIN(cnt0, cnt1) == 0) 
		return 0;
	return 1.0 * res2 / global_avg_sum_bin_hash; 
}

float get_share_mate(struct asm_graph_t *g, int i0, int i1)
{
	int rev_i0 = g->edges[i0].rc_id;
	struct asm_edge_t *rev_e0 = &g->edges[rev_i0];
	int score = 0;
	for (int i = 0; i < rev_e0->n_mate_contigs; i++){
		if (rev_e0->mate_contigs[i] == i1) {
			score += rev_e0->mate_barcodes[i].n_item;
		}
	}
	return score;
}

float get_share_mate_2(struct asm_graph_t *g, int i0, int i1)
{
	int rev_i0 = g->edges[i0].rc_id;
	struct asm_edge_t *rev_e0 = &g->edges[rev_i0];
	int score = 0;
	for (int i = 0; i < rev_e0->n_mate_contigs_2; i++){
		if (rev_e0->mate_contigs_2[i] == i1) {
			score += rev_e0->mate_barcodes_2[i].n_item;
		}
	}
	return score;
}
