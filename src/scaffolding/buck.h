#ifndef SCAFFOLDS_CHECK_BUCK_H
#define SCAFFOLDS_CHECK_BUCK_H

#include "assembly_graph.h"
#include "k31hash.h"

float get_score_multiple_buck(struct asm_graph_t *g, struct asm_edge_t *e, struct barcode_hash_t *b_left, struct barcode_hash_t *b_right);
float get_score_bucks(struct barcode_hash_t *buck0, struct barcode_hash_t *buck1, float edge0_cov,
		float edge1_cov);
int check_qualify_buck(struct asm_graph_t *g, struct asm_edge_t *e, int b, float avg_bin_hash);
#endif
