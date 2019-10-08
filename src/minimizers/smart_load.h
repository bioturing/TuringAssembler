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

void smart_construct_read_index(struct read_path_t *rpath, khash_t(bcpos) *h);

void smart_load_barcode(struct opt_proc_t *opt);

void stream_filter_read(struct read_path_t *ref, khash_t(bcpos) *dict,
                        uint64_t *shared, int n_shared, char **buf1, char **buf2,
                        uint64_t *size1, uint64_t *size2);

#endif //SRC_SMART_LOAD_H
