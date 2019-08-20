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

float get_bc_score(int count_share, int size0, int size1, float avg_bin_hash)
{
	if (MIN(size0, size1) < avg_bin_hash/15)
		return 0;
	return 1.0*count_share/MIN(size0 , size1);
}

float get_share_barcode(struct barcode_hash_t *buck0, struct barcode_hash_t *buck1, float edge0_cov,
		float edge1_cov, float avg_bin_hash)
{
	int cnt0 = 0, cnt1 = 0, res2 = 0;
	float ratio0 = edge0_cov / global_genome_coverage, ratio1 = edge1_cov / global_genome_coverage;

	for (uint32_t i = 0; i < buck0->size; ++i) {
		if ((buck0->keys[i]) != (uint64_t)(-1)){
			uint32_t tmp = barcode_hash_get(buck1, buck0->keys[i]);
			if (tmp != BARCODE_HASH_END(buck1) && buck1->keys[tmp] != (uint64_t)(-1)) {
				res2++;
			}
		}
	}
	VERBOSE_FLAG(3, "res %d cnt0 %d cnt1 %d \n", res2, cnt0, cnt1 );
	cnt0 = buck0->n_item;
	cnt1 = buck1->n_item;
	return get_bc_score(res2, cnt0, cnt1, avg_bin_hash);
}

float get_share_mate(struct asm_graph_t *g, int i0, int i1)
{
    int res = 0;
    struct barcode_hash_t *b0 = &g->edges[i0].barcodes_scaf2, *b1 = &g->edges[i1].barcodes_scaf2;
    for (uint32_t i = 0; i < b0->size; i++) {
        if (b0->keys[i] != (uint64_t)(-1)) {
            uint32_t  tmp = barcode_hash_get(b1, b0->keys[i]);
            if (tmp != BARCODE_HASH_END(b1) && b1->keys[tmp] != (uint64_t)(-1)) {
                res+=1;
            }
        }
    }
    return (float) res;
}

