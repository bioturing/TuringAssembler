#include "assembly_graph.h"
#include "verbose.h"
#include "scaffolding/global_params.h"
#include "scaffolding/bin_hash.h"

int get_amount_hole(struct asm_graph_t *g, struct asm_edge_t *e, int b)
{
	int res = 0, l = b * g->bin_size, r = (b+1) * g->bin_size, sum_holes = 0;
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

int check_qualify_buck(struct asm_graph_t *g, struct asm_edge_t *e, int b, float avg_bin_hash)
{
	if (get_amount_hole(g, e, b)  > 0.7*g->bin_size) {
		VERBOSE_FLAG(2, "NNNNN size is to big ");
		return 0;
	}
	int cnt = 0, cov = global_genome_coverage, normal_count = (g->bin_size - g->ksize +1) * cov;

	struct barcode_hash_t *buck = &e->bucks[b];
	for (uint32_t i = 0; i < buck->size; ++i) {
		if (buck->cnts[i] != (uint32_t)(-1)) {
			cnt += buck->cnts[i];
		}
	}
	assert(normal_count != 0);
	if  (cnt > 2 * normal_count || cnt < normal_count * 0.5) 
	{
		VERBOSE_FLAG(2, "count hash is abnormal: %d %d\n", cnt, normal_count);
		return 0;
	}
	return 1;
}

int check_count_hash_buck(struct asm_graph_t *g, struct asm_edge_t *e, struct barcode_hash_t *buck,  float avg_bin_hash)
{
	int cnt = 0, cov = global_genome_coverage, normal_count = (g->bin_size - g->ksize +1) * cov;
	for (uint32_t i = 0; i < buck->size; ++i) {
		if (buck->cnts[i] != (uint32_t)(-1)) {
			cnt += buck->cnts[i];
		}
	}
	assert(normal_count != 0);
	if  (cnt > 2 * normal_count || cnt < normal_count * 0.5) 
	{
		VERBOSE_FLAG(2, "count hash is abnormal: %d %d\n", cnt, normal_count);
		return 0;
	}
	return 1;
}

float get_score_bucks(struct barcode_hash_t *buck0, struct barcode_hash_t *buck1, float edge0_cov,
		float edge1_cov)
{
	const int thres_cnt = global_thres_count_kmer;
	int cnt0 = 0, cnt1 = 0, res2 = 0;
	float ratio0 = edge0_cov / global_genome_coverage, ratio1 = edge1_cov / global_genome_coverage;

//	VERBOSE_FLAG(0, "ratio %f %f\n" , ratio0, ratio1);
	for (uint32_t i = 0; i < buck1->size; ++i) {
		if (buck1->cnts[i] != (uint32_t)(-1)) {
			if (buck1->cnts[i] >= (uint32_t)thres_cnt) {
				cnt1+= buck1->cnts[i];
			}
		}
	}

	for (uint32_t i = 0; i < buck0->size; ++i) {
		if ((buck0->keys[i]) != (uint64_t)(-1)){
			if (buck0->cnts[i] >= (uint32_t)thres_cnt) {
				cnt0 += buck0->cnts[i];
				uint32_t tmp = barcode_hash_get(buck1, buck0->keys[i]);
				if (tmp != BARCODE_HASH_END(buck1) && buck1->cnts[tmp] >= (uint32_t)thres_cnt) {
					res2+= MIN(ratio0 * buck0->cnts[i] , ratio1 * buck1->cnts[tmp]);
//					VERBOSE_FLAG(3, "MIN %d\n", MIN(buck0->cnts[i], buck1->cnts[tmp]));
				}
			}
		}
	}
	VERBOSE_FLAG(3, "res %d cnt0 %d cnt1 %d \n", res2, cnt0, cnt1 );
	if (MIN(cnt0, cnt1) == 0) 
		return 0;
	return 1.0 * res2 / global_avg_sum_bin_hash; 
}

float get_score_multiple_buck(struct asm_graph_t *g, struct asm_edge_t *e, struct barcode_hash_t *b_left, struct barcode_hash_t *b_right)
{
	float avg_bin_hash = get_avg_unique_bin_hash(g);
	float res = 0;
	int count = 0;
	float cov_e = __get_edge_cov(e, g->ksize);
	for (int i = 0; i < 3; i++) if (check_count_hash_buck(g, e, &b_left[i], avg_bin_hash)) {
		for(int j = 0; j < 3; j++) if (check_count_hash_buck(g, e, &b_right[j], avg_bin_hash)) {
			count ++;
			float t = get_score_bucks(&b_left[i], &b_right[j], cov_e, cov_e);
//			VERBOSE_FLAG(log_check_contig, "score %f\n", t);
			res += t;
		}
	}
	if (count == 0)
		return 0;
	return res/count;
}

