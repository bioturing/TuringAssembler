#ifndef __K31_COUNT_H__
#define __K31_COUNT_H__

#include "attribute.h"
#include "k31hash.h"

void k31_test_process(struct opt_count_t *opt);

void build_k31_table_lazy(struct opt_count_t *opt, struct k31hash_t *h, int ksize);

#endif
