#ifndef _KMER_COUNT_H_
#define _KMER_COUNT_H_

#include "attribute.h"
#include "kmhash.h"

void kmer_test_process(struct opt_count_t *opt);

void kmer_fastq_count(struct opt_count_t *opt);

struct kmhash_t *count_kmer(struct opt_count_t *opt);

struct kmhash_t *build_kmer_table_lazy(struct opt_count_t *opt);

#endif
