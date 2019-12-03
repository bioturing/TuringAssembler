#ifndef SCAFFOLDS_CHECK_BUCK_H
#define SCAFFOLDS_CHECK_BUCK_H

#include "assembly_graph.h"

float get_share_barcode(struct barcode_hash_t *buck0, struct barcode_hash_t *buck1, float edge0_cov,
		float edge1_cov, float avg_bin_hash);
float get_share_mate(struct asm_graph_t *g, int i0, int i1);
float get_bc_score(int count_share, int size0, int size1, float avg_bin_hash, int src, int des);
#endif
