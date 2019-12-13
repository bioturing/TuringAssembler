#include "assembly_graph.h"
#include "verbose.h"
#include "scaffolding/global_params.h"
#include "scaffolding/bin_hash.h"
#include "utils.h"

int get_amount_hole(struct asm_graph_t *g, struct asm_edge_t *e)
{
	int res = 0, l = 0, r = MIN_CONTIG_BARCODE, sum_holes = 0;
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
		log_debug("NNNNN size is to big ");
		return 0;
	}
	return 1;
}
#define MIN_SHARE_BARCODE 75
float get_bc_score(int count_share, int size0, int size1, float avg_bin_hash, int src, int des)
{
//	if (MIN(size0, size1) < avg_bin_hash/15)
//		return 0;
	if (size0 < MIN_SHARE_BARCODE|| size1 < MIN_SHARE_BARCODE){
		log_warn("Two edges %d %d have low no. of barcode %d %d, threshold %d", src, des, size0, size1, MIN_SHARE_BARCODE);

		return 0;
	}
	return 1.0*count_share/MIN(size0 , size1);
}


