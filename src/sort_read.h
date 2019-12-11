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
struct readsort_bundle_t {
	struct dqueue_t *q;
	char prefix[MAX_PATH];
	int64_t sm;
};

void sort_read(struct opt_proc_t *opt, struct read_path_t *rpath);
void *biot_buffer_iterator(void *data);
void *ust_buffer_iterator(void *data);
void *x10_buffer_iterator(void *data);
RS_PROTO(read_index, struct read_index_t);

#endif /* __SORT_READ_H__ */
