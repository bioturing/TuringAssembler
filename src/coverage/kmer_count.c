//
// Created by BioTuring on 2019-12-08.
//
#include "../utils.h"

#include "../assembly_graph.h"
#include "../minimizers/minimizers.h"
#include "../minimizers/count_barcodes.h"
#include "../fastq_producer.h"
#include "../atomic.h"

#include "../log.h"

#define KMER_SIZE_COVERAGE 31
#define MAX_KMER_COUNT 999

#define __reverse_bit(a, b) \
	b = ((a & 0x5555555555555555) << 1)  | ((a >> 1)  & 0x5555555555555555); \
	b = ((b & 0x3333333333333333) << 2)  | ((b >> 2)  & 0x3333333333333333); \
	b = ((b & 0x0f0f0f0f0f0f0f0f) << 4)  | ((b >> 4)  & 0x0f0f0f0f0f0f0f0f); \
	b = ((b & 0x00ff00ff00ff00ff) << 8)  | ((b >> 8)  & 0x00ff00ff00ff00ff); \
	b = ((b & 0x0000ffff0000ffff) << 16) | ((b >> 16) & 0x0000ffff0000ffff); \
	b = ((b & 0x00000000ffffffff) << 32) | ((b >> 32) & 0x00000000ffffffff);

struct cov_bundle_t {
    struct dqueue_t *q;
    struct mini_hash_t *table;
};

/**
 * Get the i kmer of the encoded DNA sequence s (from minimizers/minimizers.c)
 * @param s DNA sequence s
 * @param i pos
 * @param k k size
 * @return encoded kmer
 */
static inline uint64_t get_km_i_bin(uint32_t *s, int i, int k)
{
	uint64_t km, c;
	int j;
	int pad = (32 - k - 1)*2;
	km = 0;
	for (j = 0; j < k; ++j) {
		c = (uint64_t)__binseq_get(s, i + j);
		km |= c;
		km <<= 2;
	}
	km <<= pad;
	return km;
}

static inline uint64_t get_km_i_str(char *s, int i, int k)
{
	uint64_t km, c;
	int j;
	int pad = (32 - k - 1)*2;
	km = 0;
	for (j = 0; j < k; ++j) {
		c = (uint64_t)nt4_table[s[i + j]];
		km |= c;
		km <<= 2;
	}
	km <<= pad;
	return km;
}

//Only for consecutive edge
int index_bin_edge(struct mini_hash_t **table, struct asm_edge_t e)
{
	uint32_t i;
	uint64_t c, cnt = 0;
	uint64_t *slot;
	int pad = (32 - KMER_SIZE_COVERAGE - 1)*2;
	uint64_t km = get_km_i_bin(e.seq, 0, KMER_SIZE_COVERAGE);
	for (i = 0 ; i < e.seq_len - KMER_SIZE_COVERAGE + 1; ++i) {
		c = (uint64_t)__binseq_get(e.seq, i + KMER_SIZE_COVERAGE - 1);
		km |= ((uint64_t) c << (pad + 2));
		slot = mini_put(table, km);
		if (*slot == 0)
			cnt++;
		km <<= 2;
	}
	return cnt;
}

int get_and_add_kmer(struct mini_hash_t *table, struct read_t read)
{
	uint32_t i;
	uint64_t c, cnt = 0;
	uint64_t *slot;
	int pad = (32 - KMER_SIZE_COVERAGE - 1)*2;
	uint64_t km = get_km_i_str(read.seq, 0, KMER_SIZE_COVERAGE);
	uint64_t rev;
	for (i = 0 ; i < read.len - KMER_SIZE_COVERAGE + 1; ++i) {
		c = (uint64_t) nt4_table[read.seq[i + KMER_SIZE_COVERAGE - 1]];
		km |= ((uint64_t) c << (pad + 2));
		__reverse_bit(km, rev);
		rev ^= 0xf;
		rev <<= (pad + 2);
		slot = mini_get(table, km);
		if (slot != (uint64_t *)EMPTY_SLOT)
			atomic_add_and_fetch64(slot, 1);
		slot = mini_get(table, rev);
		if (slot != (uint64_t *)EMPTY_SLOT)
			atomic_add_and_fetch64(slot, 1);
		km <<= 2;
	}
	return cnt;
}

void add_cnt_to_graph(struct asm_graph_t *g, struct mini_hash_t *kmer_table)
{
	int i, j;
	uint64_t *slot;
	uint64_t km, c;
	int pad = (32 - KMER_SIZE_COVERAGE - 1)*2;
	for (i = 0; i < g->n_e; ++i) {
		g->edges[i].count = 0;
		km = get_km_i_bin(g->edges[i].seq, 0, KMER_SIZE_COVERAGE);
		for (j = 0 ; j < g->edges[i].seq_len - KMER_SIZE_COVERAGE + 1; ++j) {
			c = (uint64_t)__binseq_get(g->edges[i].seq, j + KMER_SIZE_COVERAGE - 1);
			km |= ((uint64_t) c << (pad + 2));
			slot = mini_get(kmer_table, km);
			if (slot != (uint64_t *)EMPTY_SLOT) // TODO:why?
				g->edges[i].count += __min(*slot, MAX_KMER_COUNT); 
			km <<= 2;
		}
	}
}

struct mini_hash_t * construct_edges_hash(struct asm_graph_t *g)
{
	int64_t i, kmer_cnt = 0;
	struct mini_hash_t *table;
	init_mini_hash(&table, INIT_PRIME_INDEX);
	for (i = 0; i < g->n_e; ++i) {
		if (g->edges[i].seq_len < (KMER_SIZE_COVERAGE + 1))
			continue;
		kmer_cnt += index_bin_edge(&table, g->edges[i]);
	}
	log_info("Total number of kmer %d", kmer_cnt);
	return table;
}

static inline void *kmer_count_iterator(void *data)
{
	struct cov_bundle_t *bundle = (struct cov_bundle_t *) data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2, readbc;
	struct pair_buffer_t *own_buf, *ext_buf;
	uint64_t *slot;
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
				log_error("Wrong format file");

			if (read1.len > KMER_SIZE_COVERAGE)
				get_and_add_kmer(bundle->table, read1);
			if (read2.len > KMER_SIZE_COVERAGE)
				get_and_add_kmer(bundle->table, read2);

			if (rc1 == READ_END)
				break;
		}
	}
	return NULL;
}
struct mini_hash_t *kmer_count_on_edges(struct opt_proc_t *opt)
{
	uint32_t i;
	struct asm_graph_t *g = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g, opt->in_file);
	struct mini_hash_t *kmer_table = construct_edges_hash(g);
	asm_graph_destroy(g);

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	void *(*buffer_iterator)(void *) = kmer_count_iterator;
	struct producer_bundle_t *producer_bundles = NULL;
	producer_bundles = init_fastq_pair(opt->n_threads, opt->n_files,
	                                   opt->files_1, opt->files_2);

	struct cov_bundle_t *worker_bundles; //use an arbitrary structure for worker bundle
	worker_bundles = malloc(opt->n_threads * sizeof(struct cov_bundle_t));

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].table = kmer_table;
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

	return kmer_table;
}

