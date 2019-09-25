#ifndef __BARCODE_RESOLVE2_H__
#define __BARCODE_RESOLVE2_H__
#include "assembly_graph.h"

void get_local_reads_intersect(struct read_path_t *reads, struct read_path_t *rpath,
			khash_t(bcpos) *dict, struct asm_graph_t *g,
			gint_t e1, gint_t e2, const char *prefix);
int check_large_pair_superior(struct asm_graph_t *g, gint_t e1,
							gint_t e2, gint_t e2a);
void construct_read_index(struct read_path_t *rpath, khash_t(bcpos) *h);
#endif
