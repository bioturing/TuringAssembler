#ifndef _ATTRIBUTE_H_
#define _ATTRIBUTE_H_

#include <stdint.h>

#define SIZE_1MB		1048576
#define SIZE_2MB		2097152
#define SIZE_4MB		4194304
#define SIZE_16MB		16777216
#define SIZE_128MB		134217728

#include "khash.h"

KHASH_DECLARE(kvert, uint64_t, struct kvert_info_t)

struct kvert_info_t {
	int idx;
	uint32_t cnt;
};

struct read_t {
	char *seq;
	char *name;
	char *qual;
	char *note;
	char *info;
	int len;
};

struct opt_count_t {
	int n_threads;
	int hash_size;
	int kmer_master;
	int kmer_slave;
	int filter_thres;
	int n_files;
	char **files_1, **files_2;
	char *out_dir;
};

struct producer_bundle_t {
	struct dqueue_t *q;
	int *n_consumer;
	void *stream;
	pthread_barrier_t *barrier;
	pthread_mutex_t *lock;
};

#endif
