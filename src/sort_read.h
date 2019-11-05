#ifndef __SORT_READ_H__
#define __SORT_READ_H__

#include "attribute.h"

struct readbc_t {
	uint64_t barcode;
	int64_t offset;
	int len1;
	int len2;
};

void sort_read(struct opt_proc_t *opt, struct read_path_t *rpath);

#endif /* __SORT_READ_H__ */
