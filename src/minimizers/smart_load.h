//
// Created by BioTuring on 2019-09-25.
//

#ifndef SRC_SMART_LOAD_H
#define SRC_SMART_LOAD_H
#include "khash.h"
#include "radix_sort.h"
#include "attribute.h"
#include "assembly_graph.h"
#include "barcode_resolve2.h"

#define read_index_get_key(p) ((p).r1_offset)

//#ifndef RS_IMPL_READ_INDEX
//#define RS_IMPL_READ_INDEX
//RS_IMPL(read_index, struct read_index_t, 64, 8, read_index_get_key);
//#endif

void rs_sort_read_index(struct read_index_t *beg, struct read_index_t *end);

void smart_construct_read_index(struct read_path_t *rpath, khash_t(bcpos) *h);

void smart_load_barcode(struct opt_proc_t *opt);

void stream_filter_read(struct read_path_t *ref, khash_t(bcpos) *dict,
                        uint64_t *shared, int n_shared, char **buf1, char **buf2,
                        uint64_t *size1, uint64_t *size2);

#endif //SRC_SMART_LOAD_H
