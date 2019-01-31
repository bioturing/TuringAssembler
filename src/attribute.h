#ifndef _ATTRIBUTE_H_
#define _ATTRIBUTE_H_

#include <stdint.h>

#define SIZE_1MB		1048576
#define SIZE_2MB		2097152
#define SIZE_4MB		4194304
#define SIZE_16MB		16777216
#define SIZE_128MB		134217728

#define HM_COFF_1			UINT64_C(0x87c37b91114253d5)
#define HM_COFF_2			UINT64_C(0x4cf5ad432745937f)

#define HM_MIX_1			UINT64_C(0xbf58476d1ce4e5b9)
#define HM_MIX_2			UINT64_C(0x94d049bb133111eb)

#define KMFLAG_EMPTY			0
#define KMFLAG_OLD			1
#define KMFLAG_NEW			2
#define KMFLAG_LOADING			3

#define KMHASH_MAX_SIZE			UINT64_C(0x400000000)
#define KMHASH_SINGLE_RESIZE		UINT64_C(0x100000)

typedef uint64_t kmint_t;

#define __round_up_kmint(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16,	       \
				 (x) |= (x) >> 32, ++(x))

#define __rotl64(x, r) (((x) << (r)) | ((x) >> (64 - (r))))

static inline uint64_t __hash_int(uint64_t k)
{
	uint64_t x = k;
	x = (x ^ (x >> 30)) * HM_MIX_1;
	x = (x ^ (x >> 27)) * HM_MIX_2;
	x ^= (x >> 31);
	return x;
}

static inline kmint_t estimate_probe_3(kmint_t size)
{
	kmint_t s, i;
	i = s = 0;
	while (s < size) {
		++i;
		s += i * i * i * 64;
	}
	return i;
}

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
