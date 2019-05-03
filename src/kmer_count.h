#ifndef __K31_COUNT_H__
#define __K31_COUNT_H__

#include "attribute.h"
#include "k31hash.h"
#include "k63hash.h"

void build_k31_table_from_scratch(struct opt_count_t *opt, struct k31hash_t *h, int ksize);
void build_k31_table_from_k31_table(struct opt_count_t *opt, struct k31hash_t *dst,
		struct k31hash_t *src, int ksize_dst, int ksize_src);

void build_k63_table_from_scratch(struct opt_count_t *opt, struct k63hash_t *h, int ksize);
void build_k63_table_from_k31_table(struct opt_count_t *opt, struct k63hash_t *dst,
		struct k31hash_t *src, int ksize_dst, int ksize_src);
void build_k63_table_from_k63_table(struct opt_count_t *opt, struct k63hash_t *dst,
		struct k63hash_t *src, int ksize_dst, int ksize_src);

#endif
