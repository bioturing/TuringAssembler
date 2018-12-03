#ifndef __FASTQ_READER_H__
#define __FASTQ_READER_H__

#include <pthread.h>

#include "dqueue.h"
#include "get_buffer.h"

struct pair_buffer_t {
	char *buf1;
	char *buf2;
	int input_format;
};

struct producer_bundle_t *init_fastq_PE(struct opt_count_t *opt);

void *fastq_PE_producer(void *data);

void free_fastq_PE(struct producer_bundle_t *bundles, int n);

#endif
