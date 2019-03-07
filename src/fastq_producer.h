#ifndef __FASTQ_READER_H__
#define __FASTQ_READER_H__

#include <pthread.h>

#include "attribute.h"
#include "dqueue.h"
#include "get_buffer.h"

struct pair_buffer_t {
	char *buf1;
	char *buf2;
	int input_format;
};

struct pair_buffer_t *init_pair_buffer();

void free_pair_buffer(struct pair_buffer_t *p);

// struct producer_bundle_t *init_fastq_PE(struct opt_count_t *opt);
struct producer_bundle_t *init_fastq_PE(int n_thread, int n_files,
						char **files_1, char **files_2);

void *fastq_PE_producer(void *data);

void free_fastq_PE(struct producer_bundle_t *bundles, int n);

#endif
