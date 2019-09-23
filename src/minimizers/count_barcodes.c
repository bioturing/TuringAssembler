/*
 * =====================================================================================
 *
 *       Filename:  count_barcodes.c
 *
 *    Description:  Count the frequency of each barcode
 *
 *        Version:  1.0
 *        Created:  09/22/2019 10:55:15 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Phan Thai Duy Tan (ptdtan), tan@bioturing.com
 *   Organization:  BioTuring Dev Team
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <stdbool.h>

#include "fastq_producer.h"
#include "attribute.h"
#include "dqueue.h"
#include "get_buffer.h"
#include "utils.h"
#include "verbose.h"
#include "atomic.h"

pthread_mutex_t lock_key;

#define HASH_END LONG_MAX

#define BARCODES100M 100663320
#define BIG_CONSTANT(x) (x##LLU)

#if defined(_MSC_VER)

#define FORCE_INLINE	__forceinline

#include <stdlib.h>

#define ROTL32(x,y)	_rotl(x,y)
#define ROTL64(x,y)	_rotl64(x,y)

#define BIG_CONSTANT(x) (x)

// Other compilers

#else	// defined(_MSC_VER)

#define	FORCE_INLINE inline __attribute__((always_inline))

inline uint32_t rotl32(uint32_t x, int8_t r)
{
	return (x << r) | (x >> (32 - r));
}

inline uint64_t rotl64(uint64_t x, int8_t r)
{
	return (x << r) | (x >> (64 - r));
}

#define	ROTL32(x,y)	rotl32(x,y)
#define ROTL64(x,y)	rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)

#endif // !defined(_MSC_VER)

static inline uint64_t fmix64(uint64_t k)
{
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xff51afd7ed558ccd);
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
	k ^= k >> 33;

	return k;
}

static struct readbc_t {
    uint64_t barcode;
    int64_t offset;
    int len1;
    int len2;
};

static struct readsort_bundle_t {
    struct dqueue_t *q;
    char prefix[MAX_PATH];
    int64_t sm;
};

struct mini_hash_t {
    uint64_t *h;
    uint64_t *key;
    unsigned int size;
};

struct mini_hash_t *init_mini_hash()
{
	struct mini_hash_t *h_table = malloc(sizeof(struct mini_hash_t));
	h_table->h = calloc(BARCODES100M, sizeof(uint64_t));
	h_table->key = calloc(BARCODES100M, sizeof(uint64_t));
	h_table->size = BARCODES100M;
	memset(h_table->h, 0, h_table->size);
	memset(h_table->key, HASH_END, h_table->size);
	return h_table;
}

static struct mini_hash_t *h_table;

static inline uint64_t get_barcode_biot(char *s, struct read_t *r)
{
	if (s == NULL) {
		r->seq = r->qual = NULL;
		return (uint64_t)-1;
	}
	int i, len = 0;
	uint64_t ret = 0;
	char *p;
	p = strstr(s, "BX:Z:");
	if (p == NULL) {
		r->seq = r->qual = NULL;
		return (uint64_t)-1;
	}
	r->seq = p += 5;
	for (i = 0; p[i] && !__is_sep(p[i]); ++i) {
		ret = ret * 5 + nt4_table[(int)p[i]];
		++len;
	}
	p = strstr(s, "QB:Z:");
	if (p == NULL)
		r->qual = NULL;
	else
		r->qual = p + 5;
	r->len = len;
	return ret;
}

static inline uint64_t getblock64(const uint64_t * p, int i)
{
	return p[i];
}

static inline uint64_t MurmurHash3_x64_64(const uint8_t *data, const int len)
{
	int n_blocks = len >> 4;
	uint64_t h1 = BIG_CONSTANT(0x13097);
	uint64_t h2 = BIG_CONSTANT(0x13097);

	const uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
	const uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

	const uint64_t *blocks = (const uint64_t *)(data);

	int i;
	for (i = 0; i < n_blocks; ++i) {
		uint64_t k1 = getblock64(blocks, i << 1);
		uint64_t k2 = getblock64(blocks, (i << 1) + 1);
		k1 *= c1; k1 = ROTL64(k1, 31); k1 *= c2; h1 ^= k1;
		h1 = ROTL64(h1, 27); h1 += h2; h1 = h1 * 5 + 0x52dce729;
		k2 *= c2; k2 = ROTL64(k2, 33); k2 *= c1; h2 ^= k2;
		h2 = ROTL64(h2, 31); h2 += h1; h2 = h2 * 5 + 0x38495ab5;
	}

	const uint8_t *tail = (const uint8_t *)(data + (n_blocks << 4));
	uint64_t k1 = 0;
	uint64_t k2 = 0;

	switch (len & 15) {
		case 15: k2 ^= ((uint64_t)tail[14]) << 48;
		case 14: k2 ^= ((uint64_t)tail[13]) << 40;
		case 13: k2 ^= ((uint64_t)tail[12]) << 32;
		case 12: k2 ^= ((uint64_t)tail[11]) << 24;
		case 11: k2 ^= ((uint64_t)tail[10]) << 16;
		case 10: k2 ^= ((uint64_t)tail[ 9]) << 8;
		case  9: k2 ^= ((uint64_t)tail[ 8]) << 0;
			k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

		case  8: k1 ^= ((uint64_t)tail[ 7]) << 56;
		case  7: k1 ^= ((uint64_t)tail[ 6]) << 48;
		case  6: k1 ^= ((uint64_t)tail[ 5]) << 40;
		case  5: k1 ^= ((uint64_t)tail[ 4]) << 32;
		case  4: k1 ^= ((uint64_t)tail[ 3]) << 24;
		case  3: k1 ^= ((uint64_t)tail[ 2]) << 16;
		case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;
		case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;
			k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
	}

	h1 ^= len; h2 ^= len;

	h1 += h2;
	h2 += h1;

	h1 = fmix64(h1);
	h2 = fmix64(h2);

	h1 += h2;
	h2 += h1;

	return h2;
}

void mini_inc(uint8_t *data, int len)
{
	uint64_t key = MurmurHash3_x64_64(data, len);
	uint64_t mask = h_table->size - 1;
	uint64_t slot = key % mask;
	uint64_t d = *((uint64_t *)data);
	uint64_t prev = atomic_val_CAS64(h_table->h + slot, 0, 1);
	if (!prev) { // slot is empty -> fill in
		h_table->key[slot] = d;
	} else {
		uint64_t probe = slot + 1;
		while (!(prev = atomic_val_CAS64(h_table->h + probe, 0, 1)) && !atomic_bool_CAS64(h_table->key + probe, d, d)) { //TODO: not this condition
			probe = (++probe) % mask;
		}
		if (probe == slot)
			__ERROR("No more slot in the hash table!");
		if (!prev) { //room at probe is empty -> fill in
			h_table->key[probe] = d;
		} else{
			atomic_add_and_fetch64(h_table->h + probe, 1);
		}
	}
}

void mini_print()
{
	FILE *fp = fopen("barcode_frequencies.txt", "w");
	int i;
	for (i = 0; i < h_table->size; ++i) {
		if (h_table->h[i])
		fprintf(fp, "%lld\t%d\n", h_table->key[i], h_table->h[i]);
	}
}

static inline void *biot_buffer_iterator_simple(void *data);

void count_bx_freq(struct opt_proc_t *opt, struct read_path_t *r_path)
{
	pthread_mutex_init(&lock_key, NULL);
	h_table = init_mini_hash();

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	void *(*buffer_iterator)(void *) = biot_buffer_iterator_simple;
	struct producer_bundle_t *producer_bundles = NULL;
	if (opt->lib_type != LIB_TYPE_BIOT) {
		__ERROR("Only accept 'bioturing' library at this stage, Sorry!");
	}
	producer_bundles = init_fastq_pair(opt->n_threads, opt->n_files,
		                                   opt->files_1, opt->files_2);
	buffer_iterator = biot_buffer_iterator_simple;
	struct readsort_bundle_t *worker_bundles; //use an arbitrary structure for worker bundle
	worker_bundles = malloc(opt->n_threads * sizeof(struct readsort_bundle_t));
	int i;
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_producer,
		               producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, buffer_iterator,
		               worker_bundles + i);

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	free_fastq_pair(producer_bundles, opt->n_files);
	free(worker_bundles);

	mini_print();

}

static inline void *biot_buffer_iterator_simple(void *data)
{
	struct readsort_bundle_t *bundle = (struct readsort_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2, readbc;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_trip_buffer();

	char *R1_buf, *R2_buf;
	int pos1, pos2, rc1, rc2, input_format;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf)
			break;
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		R1_buf = ext_buf->R1_buf;
		R2_buf = ext_buf->R2_buf;
		input_format = ext_buf->input_format;

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, R1_buf, &pos1) :
				get_read_from_fa(&read1, R1_buf, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, R2_buf, &pos2) :
				get_read_from_fa(&read2, R2_buf, &pos2);

			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file\n");

			/* read_name + \t + BX:Z: + barcode + \t + QB:Z: + barcode_quality + \n */
			uint64_t barcode = get_barcode_biot(read1.info, &readbc);
			int record_len, len1, len2;
			if (barcode != (uint64_t)-1) { //read doesn't have barcode
			} else {
				// any main stuff goes here
				mini_inc(&barcode, sizeof(uint64_t) / sizeof(uint8_t));
			}
			if (rc1 == READ_END)
				break;
		}
	}
	return NULL;
}
