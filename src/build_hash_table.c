//
// Created by che on 11/11/2019.
//
#include "attribute.h"
#include "utils.h"
#include "fastq_producer.h"
#include "verbose.h"
#include "build_hash_table.h"
#include "assembly_graph.h"
#include "kmhash.h"
#include "atomic.h"

void get_seq(char *seq, int start, int len, uint8_t *res)
{
	for (int i = start, i_res = 0; i < start + len; i++, i_res++) {
		res[i_res >> 2] |= ((seq[i >> 2] >> ((i & 3) << 1)) & 3) << ((i_res & 3) << 1);
	}
}

void copy_seq32_seq8(uint32_t *seq, int start, uint8_t *res, int start_res, int len)
{
	for (int i = start, i_res = start_res; i < start + len; i++, i_res++) {
		res[i_res >> 2] |= ((seq[i >> 4] >> ((i & 15) << 1)) & 3) << ((i_res & 3) << 1);
	}
}

void print_u8_seq(uint8_t *a, int len)
{
	for (int i = 0; i < len; i++) {
		printf("%c", nt4_char[(a[i >> 2] >> ((i & 3) << 1)) & 3]);
	}
	printf("\n");
}

uint8_t *compress_seq(char *a)
{
	int len = strlen(a);
	uint8_t *res = calloc(len, 1);
	for (int i = 0; i < len; i++) {
		res[i >> 2] |= nt4_table[a[i]] << ((i & 3) << 1);
	}
	return res;
}

void ust_add_big_kmer(struct read_t *r, khash_t(pair_kmer_count) *table, int ksize, pthread_mutex_t *lock)
{
	int8_t *seq = compress_seq(r->seq);
	int big_ksize = BIG_KSIZE;

	uint8_t *left = calloc((big_ksize + 3) >> 2, sizeof(uint8_t));
	for (int i = 0; i < r->len - DISTANCE_KMER; i++) {
		memset(left, 0, (big_ksize + 3) >> 2);
		int64_t res;
		get_seq(seq, i, big_ksize, left);
		res = MurmurHash3_x64_64(left, (big_ksize + 3) >> 2);
		pthread_mutex_lock(lock);
		khint_t k = kh_get(pair_kmer_count, table, res);
		if (k == kh_end(table)) {
			int tmp;
			k = kh_put(pair_kmer_count, table, res, &tmp);
			kh_value(table, k) = 0;
			assert(tmp == 1);
		}
		pthread_mutex_unlock(lock);
		atomic_add_and_fetch32(&kh_value(table, k), 1);
	}
	free(left);
}


void *get_pair_kmer_ust_iterator(void *data)
{
	struct kmer_pair_iterator_bundle_t *bundle = (struct kmer_pair_iterator_bundle_t *) data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	khash_t(pair_kmer_count) *table = bundle->table;
	own_buf = init_trip_buffer();
	int64_t sm = bundle->sm;

	int m_read;
	m_read = sm / 600;
	__round_up_32(m_read);
	char *R1_buf, *R2_buf, *I_buf;
	int pos1, pos2, posI, rc1, rc2, rcI, input_format;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf)
			break;
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = posI = 0;
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

			if (rc1 == READ_FAIL || rc2 == READ_FAIL )
				__ERROR("\nWrong format file ust\n");

			ust_add_big_kmer(&read1, table, bundle->ksize, bundle->table_lock);
			ust_add_big_kmer(&read2, table, bundle->ksize, bundle->table_lock);
			if (rc1 == READ_END)
				break;
		}
	}
}

void build_pair_kmer_table(struct opt_proc_t *opt, khash_t(pair_kmer_count) *table)
{
	struct read_path_t read_path;
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	void *(*buffer_iterator)(void *) = NULL;
	struct producer_bundle_t *producer_bundles = NULL;
	if (opt->lib_type == LIB_TYPE_BIOT || opt->lib_type == LIB_TYPE_10X) {
//		producer_bundles = init_fastq_pair(opt->n_threads, opt->n_files,
//										   opt->files_1, opt->files_2);
//		if (opt->lib_type == LIB_TYPE_BIOT)
//			buffer_iterator = biot_buffer_iterator;
//		else
//			buffer_iterator = x10_buffer_iterator;
		__ERROR("not handle yet");
	} else if (opt->lib_type == LIB_TYPE_UST || opt->lib_type == LIB_TYPE_SORTED) {
		producer_bundles = init_fastq_pair(opt->n_threads, opt->n_files,
		                                   opt->files_1, opt->files_2);
		buffer_iterator = get_pair_kmer_ust_iterator;
	} else {
		__ERROR("Wrong library format\n");
	}
	struct kmer_pair_iterator_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct kmer_pair_iterator_bundle_t));
	int i;
	int64_t sm_in_byte, sm_per_thread;
	sm_in_byte = (int64_t) opt->mmem * 1024 * 1024 * 1024;
	sm_per_thread = sm_in_byte / opt->n_threads;
	__round_up_64(sm_per_thread);
	if (sm_per_thread > sm_in_byte / opt->n_threads)
		sm_per_thread >>= 1;
	char *read_dir = malloc(MAX_PATH);
	sprintf(read_dir, "%s/reads", opt->out_dir);
	mkdir(read_dir, 0755);

	pthread_mutex_t lock;
	pthread_mutex_init(&lock, NULL);
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].sm = sm_per_thread;
		worker_bundles[i].table = table;
		worker_bundles[i].ksize = KMER_PAIR_SIZE;
		worker_bundles[i].table_lock = &lock;

		sprintf(worker_bundles[i].prefix, "%s/thread_%d", read_dir, i);
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

	if (opt->lib_type == LIB_TYPE_BIOT || opt->lib_type == LIB_TYPE_10X)
		free_fastq_pair(producer_bundles, opt->n_files);
	else if (opt->lib_type == LIB_TYPE_UST)
		free_fastq_pair(producer_bundles, opt->n_files);
	free(worker_bundles);
}

