#ifndef __FASTQ_READER_H__
#define __FASTQ_READER_H__

#include <pthread.h>

#include "attribute.h"
#include "dqueue.h"
#include "get_buffer.h"

struct producer_bundle_t {
	struct dqueue_t *q;
	int *n_consumer;
	void *stream;
	pthread_barrier_t *barrier;
	pthread_mutex_t *lock;

	void *(*buffer_init)();
	void (*buffer_free)(void *);
	int64_t (*gb_get_data)(void *, void *);

	int64_t total_size;
	int64_t *processed_size;
};

void *init_single_buffer();
void *init_trip_buffer();
void *init_pair_buffer();

void free_single_buffer(void *vp);
void free_pair_buffer(void *vp);
void free_trip_buffer(void *vp);

struct producer_bundle_t *init_fastq_single(int n_thread, int n_file, char **files);
struct producer_bundle_t *init_fastq_pair(int n_thread, int n_file,
					char **R1_files, char **R2_files);
struct producer_bundle_t *init_fastq_triple(int n_thread, int n_file,
			char **R1_files, char **R2_files, char **I_files);

void free_fastq_single(struct producer_bundle_t *bundles, int n);
void free_fastq_pair(struct producer_bundle_t *bundles, int n);
void free_fastq_triple(struct producer_bundle_t *bundles, int n);

void *fastq_producer(void *data);

#endif
