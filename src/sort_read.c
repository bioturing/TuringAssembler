#include "sort_read.h"

int estimate_ust_record_len(struct read_t *r1, struct read_t *r2, struct read_t *rI)
{
	int ret = (r1->len + 1) * 2 + 2 + (r2->len + 1) * 2 + 2;
	ret += strlen(r1->name) + strlen(r2->name) + (1 + rI->len + 6 + rI->len + 6) * 2;
	return ret;
}

void *ust_buffer_iterator(void *data)
{
	struct readsort_bundle_t *bundle = (struct readsort_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2, readI;
	struct trip_buffer_t *own_buf, *ext_buf;
	own_buf = init_trip_buffer();

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
		I_buf = ext_buf->I_buf;
		input_format = ext_buf->input_format;

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, buf1, &pos1) :
				get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, buf2, &pos2) :
				get_read_from_fa(&read2, buf2, &pos2);

			rcI = input_format = TYPE_FASTQ ?
				get_read_from_fq(&readI, I_buf, &posI) :
				get_read_from_fa(&readI, I_buf, &posI);

			if (rc1 == READ_FAIL || rc2 == READ_FAIL || rcI == READ_FAIL)
				__ERROR("\nWrong format file\n");

			/* read_name + \t + BX:Z: + barcode + \t + QB:Z: + barcode_quality + \n */
			barcode = get_barcode_ust_raw(&readI);
			record_len = strlen(read1.name) + strlen(read2.name) +
				2 * (1 + (readI.len + 6) * 2) +
				(read1.len + 1) * 2 + (read2.len + 1) * 2 + 2 * 2;
			if (record_len > buf_size)
				/* extend buf size */;
			if (record_len > buf_size - buf_len)
				/* sort and flush buffer */;
			/* append record */
			len1 = ust_add_record(&read1, &readI, buf + buf_len, input_format);
			len2 = ust_add_record(&read2, &readI, buf + buf_len + len1, input_format);

			if (rc1 == READ_END)
				break;
		}
	}
}

static inline int ust_add_record(struct read_t *r, struct read_t *rI,
						char *buf, int input_format)
{
	int len = 0, tlen, i;
	tlen = strlen(r->name);
	memcpy(buf, r->name, tlen);
	len = tlen;
	buf[len++] = '\t';
	memcpy(buf + len, "BX:Z:", 5);
	len += 5;
	memcpy(buf + len, rI->seq, rI->len);
	len += rI->len;
	buf[len++] = '\t';
	memcpy(buf + len, "QB:Z:", 5);
	len += 5;
	if (input_format == TYPE_FASTQ) {
		memcpy(buf + len, rI->qual, rI->len);
		len += rI->len;
	} else {
		for (i = 0; i < rI->len; ++i) {
			buf[len++] = nt4_table[rI->seq[i]] < 4 ? 'F' : '#';
		}
	}
	buf[len++] = '\n';
	memcpy(buf + len, r->seq, r->len);
	len += r->len;
	buf[len++] = '\n';
	memcpy(buf + len, "+\n", 2);
	len += 2;
	if (input_format == TYPE_FASTQ) {
		memcpy(buf + len, r->qual, r->len);
		len += r->len;
	} else {
		for (i = 0; i < r->len; ++i) {
			buf[len++] = nt4_table[r->seq[i]] < 4 ? 'F' : '#';
		}
	}
	buf[len++] = '\n';
	return len;
}

void sort_read(struct opt_proc_t *opt)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	struct producer_bundle_t *producer_bundles;
	if (opt->lib_type == LIB_TYPE_BIOT || opt->lib_type == LIB_TYPE_10X)
		producer_bundles = init_fastq_pair(opt->n_threads, opt->n_files,
						opt->files_1, opt->files_2);
	else if (opt->lib_type == LIB_TYPE_UST)
		producer_bundles = init_fastq_triple(opt->n_threads, opt->n_files,
				opt->files_1, opt->files_2, opt->files_I);
	struct readsort_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct readsort_bundle_t));
	pthread_mutex_t lock;
	pthread_mutex_init(&lock, NULL);
	int i;
	sm_in_byte = opt->sm * 1024 * 1024 * 1024;
	sm_per_thread = sm_in_byte / opt->n_threads;
	__round_up_32(sm_per_thread);
	if (sm_per_thread > sm_in_byte / opt->n_threads)
		sm_per_thread >>= 1;
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].sm = sm_per_thread;
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

	free(producer_threads);
	free(worker_threads);

}

void barcode_start_count(struct opt_proc_t *opt, struct bccount_bundle_t *ske)
{
	int i;

	struct bccount_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct bccount_bundle_t));
	int64_t n_reads = 0;

	pthread_mutex_t lock;
	pthread_mutex_init(&lock, NULL);
	uint64_t hash_sum;
	hash_sum = 0;
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].n_reads = &n_reads;
		worker_bundles[i].g = ske->g;
		worker_bundles[i].aux_build = ske->aux_build;
		worker_bundles[i].bwa_idx = ske->bwa_idx;
		worker_bundles[i].bwa_opt = ske->bwa_opt;
		worker_bundles[i].barcode_calculator = ske->barcode_calculator;
		worker_bundles[i].lock = &lock;
		worker_bundles[i].hash_sum = &hash_sum;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	// __VERBOSE("hash sum = %lu\n", hash_sum);

}

