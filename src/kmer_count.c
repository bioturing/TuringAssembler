#include <stdlib.h>
#include <string.h>

#include "dqueue.h"
#include "get_buffer.h"
#include "kmer_count.h"
#include "kmhash.h"
#include "utils.h"
#include "verbose.h"

struct pair_buffer_t {
	char *buf1;
	char *buf2;
	int input_format;
};

struct producer_bundle_t {
	int *n_consumer;
	void *stream;
	pthread_barrier_t *barrier;
	pthread_mutex_t *lock;
	struct dqueue_t *q;
	struct kmhash_t *h;
};

struct worker_stat_t {
	int64_t r1_sum;
	int64_t r2_sum;
	int nread;
};

struct worker_bundle_t {
	struct dqueue_t *q;
	struct kmhash_t *h;
	pthread_mutex_t *lock_count;
	// pthread_barrier_t *barrier_hash;
	struct worker_stat_t *global_stat;
	int n_threads;
	pthread_mutex_t *lock_hash;
	int ksize;
	int count_reverse;
	// int thread_no;
};

struct pair_buffer_t *init_pair_buffer()
{
	struct pair_buffer_t *ret = malloc(sizeof(struct pair_buffer_t));
	ret->buf1 = malloc(BUF_SIZE + 1);
	ret->buf2 = malloc(BUF_SIZE + 1);
	return ret;
}

void free_pair_buffer(struct pair_buffer_t *p)
{
	if (!p) return;
	free(p->buf1);
	free(p->buf2);
	free(p);
}

struct dqueue_t *init_dqueue_PE(int cap)
{
	struct dqueue_t *ret = init_dqueue(cap);
	struct pair_buffer_t *p;
	int i;
	for (i = 0; i < cap; ++i) {
		p = init_pair_buffer();
		d_enqueue_out(ret, p);
	}
	return ret;
}

void *producer_worker(void *data)
{
	struct producer_bundle_t *bundle = (struct producer_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct gb_pair_data *input_stream = bundle->stream;
	struct pair_buffer_t *own_buf = init_pair_buffer();
	struct pair_buffer_t *external_buf;
	int64_t offset;
	int64_t n_frag = 0;

	while ((offset = gb_get_pair(input_stream, &own_buf->buf1, &own_buf->buf2)) != -1) {
		own_buf->input_format = input_stream->type;
		external_buf = d_dequeue_out(q);
		d_enqueue_in(q, own_buf);
		own_buf = external_buf;
		++n_frag;
	}
	free_pair_buffer(own_buf);

	int cur;
	pthread_barrier_wait(bundle->barrier);
	while (1) {
		pthread_mutex_lock(bundle->lock);
		cur = *(bundle->n_consumer);
		if (*(bundle->n_consumer) > 0)
			--*(bundle->n_consumer);
		pthread_mutex_unlock(bundle->lock);
		if (cur == 0)
			break;
		external_buf = d_dequeue_out(q);
		free_pair_buffer(external_buf);
		d_enqueue_in(q, NULL);
	}

	pthread_exit(NULL);
}

void add_sum_read(struct read_t *r1, struct read_t *r2, struct worker_stat_t *stat)
{
	++stat->nread;
	int i;
	for (i = 0; i < r1->len; ++i) {
		stat->r1_sum += r1->seq[i];
		stat->r2_sum += r2->seq[i];
	}
}

void count_kmer_read(struct read_t *r, struct worker_bundle_t *bundle)
{
	struct kmhash_t *h = bundle->h;
	pthread_mutex_t *lock_hash = bundle->lock_hash;
	int i, k, last, c, len, lmc;
	char *seq;
	len = r->len;
	seq = r->seq;
	k = bundle->ksize;
	kmkey_t knum, krev, kmask;
	kmask = ((kmkey_t)1 << (k << 1)) - 1;
	knum = krev = 0;
	last = 0;
	lmc = (k - 1) << 1;
	for (i = 0; i < len; ++i) {
		c = nt4_table[(int)seq[i]];
		knum = (knum << 2) & kmask;
		krev = (krev >> 2) & kmask;
		if (c < 4) {
			knum |= c;
			krev |= (kmkey_t)(c ^ 3) << lmc;
			++last;
		} else {
			last = 0;
		}
		if (last >= k) {
			if (knum < krev)
				kmhash_inc_val_wrap(h, knum, lock_hash);
			else
				kmhash_inc_val_wrap(h, krev, lock_hash);
			//	kmhash_inc_val_wrap(h, knum, lock_hash);
		}
	}
	// if (bundle->count_reverse) {
	// 	knum = 0;
	// 	last = 0;
	// 	for (i = 0; i < len; ++i) {
	// 		c = nt4_table[seq[len - i - 1]];
	// 		knum = (knum << 2) & kmask;
	// 		if (c < 4) {
	// 			knum |= (c ^ 3);
	// 			++last;
	// 		} else {
	// 			last = 0;
	// 		}
	// 		if (last >= k)
	// 			kmhash_inc_val_wrap(h, knum, lock_hash);
	// 	}
	// }
}

void *count_worker(void *data)
{
	struct worker_bundle_t *bundle = (struct worker_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	// struct kmhash_t *h = bundle->h;
	// pthread_mutex_t *lock_hash = bundle->lock_hash;

	// struct kmthread_bundle_t kmhash_bundle;
	// kmhash_bundle.thread_no = bundle->thread_no;
	// kmhash_bundle.n_threads = bundle->n_threads;
	// kmhash_bundle.barrier = bundle->barrier_hash;
	// kmhash_init(h, bundle->init_hash_size, 0, &kmhash_bundle);

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	struct worker_stat_t stat;
	memset(&stat, 0, sizeof(struct worker_stat_t));

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf)
			break;
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		buf1 = ext_buf->buf1;
		buf2 = ext_buf->buf2;
		input_format = ext_buf->input_format;
		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, buf1, &pos1) :
				get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, buf2, &pos2) :
				get_read_from_fa(&read2, buf2, &pos2);


			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file\n");

			add_sum_read(&read1, &read2, &stat);
			// add_sum_reads(&read2, &(own_result.sum_reads));
			// align_chromium_read(&read1, &read2, &kmhash_bundle, bundle);
			count_kmer_read(&read1, bundle);
			count_kmer_read(&read2, bundle);

			if (rc1 == READ_END)
				break;
		}
	}
	pthread_mutex_lock(bundle->lock_count);
	bundle->global_stat->r1_sum += stat.r1_sum;
	bundle->global_stat->r2_sum += stat.r2_sum;
	bundle->global_stat->nread += stat.nread;
	pthread_mutex_unlock(bundle->lock_count);

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

struct kmhash_t *count_kmer(struct opt_count_t *opt)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	struct dqueue_t *q;
	q = init_dqueue_PE(opt->n_threads * 2);
	int n_consumer;
	n_consumer = opt->n_threads;

	struct producer_bundle_t *producer_bundles;
	pthread_t *producer_threads;

	producer_bundles = malloc(opt->n_files * sizeof(struct producer_bundle_t));
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));

	pthread_mutex_t producer_lock;
	pthread_barrier_t producer_barrier;

	pthread_mutex_init(&producer_lock, NULL);
	pthread_barrier_init(&producer_barrier, NULL, opt->n_files);

	int i;
	for (i = 0; i < opt->n_files; ++i) {
		struct gb_pair_data *data = calloc(1, sizeof(struct gb_pair_data));
		gb_pair_init(data, opt->files_1[i], opt->files_2[i]);

		producer_bundles[i].n_consumer = &n_consumer;
		producer_bundles[i].stream = (void *)data;
		producer_bundles[i].q = q;
		producer_bundles[i].barrier = &producer_barrier;
		producer_bundles[i].lock = &producer_lock;
		pthread_create(producer_threads + i, &attr, producer_worker, producer_bundles + i);
	}

	struct worker_bundle_t *worker_bundles;
	pthread_t *worker_threads;

	struct worker_stat_t result;
	memset(&result, 0, sizeof(struct worker_stat_t));

	worker_bundles = malloc(opt->n_threads * sizeof(struct worker_bundle_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	pthread_mutex_t lock_count;
	pthread_mutex_init(&lock_count, NULL);

	struct kmhash_t *h;
	h = init_kmhash(opt->hash_size, opt->n_threads);
	// pthread_barrier_t barrier_hash;
	// pthread_barrier_init(&barrier_hash, NULL, opt->n_threads);

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = q;
		worker_bundles[i].h = h;
		worker_bundles[i].n_threads = opt->n_threads;
		// worker_bundles[i].thread_no = i;
		// worker_bundles[i].barrier_hash = &barrier_hash;
		worker_bundles[i].lock_count = &lock_count;
		worker_bundles[i].lock_hash = h->locks + i;
		worker_bundles[i].global_stat = &result;

		worker_bundles[i].ksize = opt->kmer_size;
		worker_bundles[i].count_reverse = 0;
		pthread_create(worker_threads + i, &attr, count_worker, worker_bundles + i);
	}

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	__VERBOSE_LOG("Result", "Number of read                 : %20d\n", result.nread);
	__VERBOSE_LOG("Result", "Number of kmer                 : %20llu\n", (unsigned long long)h->n_items);
	return h;
}

