#ifndef __FASTQ_REDUCER__
#define __FASTQ_REDUCER__
#include "dqueue.h"
#include "../include/bwa.h"
#include "../include/bwamem.h"
#include "attribute.h"
#define STRICT_HEAD_LEN 1000
#define BUF_LEN (1 << 26)

struct fastq_reducer_bundle_t{
	struct dqueue_t *q;
	bwaidx_t *bwa_idx;
	mem_opt_t *bwa_opt;
	pthread_mutex_t *buf_lock[2];
	pthread_mutex_t *file_lock[2];
	char *buf[2];
	int *buf_pos[2];
	FILE *f[2];
};

void fastq_reducer(struct opt_proc_t *opt, struct read_path_t *org_rpath,
		struct read_path_t *reduced_path);
void *fastq_reducer_iterator(void *data);
void reduce_read(struct read_t *r1, struct read_t *r2,
		struct fastq_reducer_bundle_t *bundle);
#endif
