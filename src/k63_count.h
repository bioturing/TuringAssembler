#ifndef __K63_COUNT_H__
#define __K63_COUNT_H__

#include "attribute.h"
#include "k63hash.h"

void k63_test_process(struct opt_count_t *opt);

void build_k63_table_lazy(struct opt_count_t *opt, struct k63hash_t *h, int ksize);

#endif
