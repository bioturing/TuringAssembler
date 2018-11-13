
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "attribute.h"
#include "dqueue.h"
#include "get_buffer.h"
#include "graph.h"
#include "kmer_count.h"
#include "kmhash.h"
#include "kseq.h"
#include "utils.h"
#include "verbose.h"


struct e_bundle_t {
	struct dqueue_t *q;
	khash_t(kvert) *h;
	int16_t *e;
	int ksize;
};

static struct pair_buffer_t *init_pair_buffer()
{
	struct pair_buffer_t *ret = malloc(sizeof(struct pair_buffer_t));
	ret->buf1 = malloc(BUF_SIZE + 1);
	ret->buf2 = malloc(BUF_SIZE + 1);
	return ret;
}

static void free_pair_buffer(struct pair_buffer_t *p)
{
	if (!p) return;
	free(p->buf1);
	free(p->buf2);
	free(p);
}

static struct dqueue_t *init_dqueue_PE(int cap)
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

static void *producer_worker(void *data)
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

void add_edge(struct read_t *r, struct e_bundle_t *bundle)
{
	khash_t(kvert) *h = bundle->h;
	int16_t *e = bundle->e;

	int i, k, last, ci, ck, len, lmc;
	char *seq;
	len = r->len;
	seq = r->seq;
	kmkey_t knum, krev, pknum, pkrev, kmask;
	khint_t ki, kk;


	k = bundle->ksize;
	kmask = ((kmkey_t)1 << (k << 1)) - 1;
	knum = krev = pknum = pkrev = 0;
	last = 0;
	lmc = (k - 1) << 1;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		knum = (knum << 2) & kmask;
		krev = (krev >> 2) & kmask;
		if (ci < 4) {
			knum |= ci;
			krev |= (kmkey_t)(ci ^ 3) << lmc;
			++last;
		} else {
			last = 0;
		}
		if (last >= k + 1) { // k + 1 for an edge
			ck = nt4_table[(int)seq[i - k]];
			// insert forward edge
			if (pknum < pkrev) {
				ki = kh_get(kvert, h, pknum);
			} else {
				ki = kh_get(kvert, h, pkrev);
				ci += 4;
			}
			if (knum < krev) {
				kk = kh_get(kvert, h, knum);
			} else {
				kk = kh_get(kvert, h, krev);
				ck += 4;
			}
			if (ki != kh_end(h) && kk != kh_end(h)) {
				__sync_fetch_and_add(e + (ki * 8 + ci), 1);
				__sync_fetch_and_add(e + (kk * 8 + ck), 1);
			}
		}
	}

}

void *edge_worker(void *data)
{
	struct e_bundle_t *bundle = (struct e_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	// khash_t(kvert) *h = bundle->h;
	// int16_t *e = bundle->e;

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

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

			add_edge(&read1, bundle);
			add_edge(&read2, bundle);
			// count_kmer_read(&read1, bundle);
			// count_kmer_read(&read2, bundle);

			if (rc1 == READ_END)
				break;
		}
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);

}

int16_t *get_edges(struct opt_count_t *opt, khash_t(kvert) *h)
{
	int nvert = kh_size(h);
	int16_t *edges = calloc(nvert * 8, sizeof(int16_t));

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

	struct e_bundle_t *worker_bundles;
	pthread_t *worker_threads;

	worker_bundles = malloc(opt->n_threads * sizeof(struct e_bundle_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = q;
		worker_bundles[i].h = h;
		worker_bundles[i].e = edges;
		// worker_bundles[i].n_threads = opt->n_threads;
		// worker_bundles[i].thread_no = i;
		// worker_bundles[i].barrier_hash = &barrier_hash;
		// worker_bundles[i].lock_count = &lock_count;
		// worker_bundles[i].lock_hash = h->locks + i;
		// worker_bundles[i].global_stat = &result;

		worker_bundles[i].ksize = opt->kmer_size;
		pthread_create(worker_threads + i, &attr, edge_worker, worker_bundles + i);
	}

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	return edges;
}
