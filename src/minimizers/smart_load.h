//
// Created by BioTuring on 2019-09-25.
//

#ifndef SRC_SMART_LOAD_H
#define SRC_SMART_LOAD_H
#include "khash.h"
#include "../radix_sort.h"
#include "../attribute.h"
#include "../assembly_graph.h"
#include "../barcode_resolve2.h"

#define read_index_get_key(p) ((p).r1_offset)

//#ifndef RS_IMPL_READ_INDEX
//#define RS_IMPL_READ_INDEX
//RS_IMPL(read_index, struct read_index_t, 64, 8, read_index_get_key);
//#endif

struct bc_hit_bundle_t{
	struct asm_graph_t *g;
	struct read_path_t *read_sorted_path;
	khash_t(bcpos) *bc_pos_dict;
	struct mm_db_edge_t *mm_edges;
};

void rs_sort_read_index(struct read_index_t *beg, struct read_index_t *end);

void smart_construct_read_index(struct read_path_t *rpath, khash_t(bcpos) *h);

void smart_load_barcode(struct opt_proc_t *opt);

void stream_filter_read(struct read_path_t *ref, khash_t(bcpos) *dict,
                        uint64_t *shared, int n_shared, char **buf1, char **buf2,
                        uint64_t *size1, uint64_t *size2);

struct mm_hits_t *get_hits_from_barcode(char *bc, struct bc_hit_bundle_t *bc_hit_bundle);
void get_bc_hit_bundle(struct opt_proc_t *opt, struct bc_hit_bundle_t *bc_hit_bundle);
void bc_hit_bundle_destroy(struct bc_hit_bundle_t *bc_hit_bundle);
#endif //SRC_SMART_LOAD_H
