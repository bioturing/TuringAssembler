#ifndef SCAFFOLDS_CHECK_BUCK_H
#define SCAFFOLDS_CHECK_BUCK_H

#include "assembly_graph.h"

float get_share_barcode(struct barcode_hash_t *buck0, struct barcode_hash_t *buck1, float edge0_cov,
		float edge1_cov);

float get_share_mate(struct asm_graph_t *g, int i0, int i1);
float get_share_mate_2(struct asm_graph_t *g, int i0, int i1);
#endif
