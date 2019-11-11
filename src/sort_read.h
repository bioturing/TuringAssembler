#ifndef __SORT_READ_H__
#define __SORT_READ_H__

#include "attribute.h"
#include "khash.h"

struct readbc_t {
    uint64_t barcode;
    int64_t offset;
    int len1;
    int len2;
};
struct read_index_t {
    int64_t r1_offset;
    int64_t r2_offset;
    int64_t r1_len;
    int64_t r2_len;
};

KHASH_MAP_INIT_INT64(bcpos, struct read_index_t);

void sort_read(struct opt_proc_t *opt, struct read_path_t *rpath);
RS_PROTO(read_index, struct read_index_t);

#endif /* __SORT_READ_H__ */
