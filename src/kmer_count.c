#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "dqueue.h"
#include "get_buffer.h"
#include "kmer_count.h"
#include "kmhash.h"
#include "kseq.h"
#include "utils.h"
#include "verbose.h"

struct precount_bundle_t {
	struct dqueue_t *q;
	struct kmhash_t *h;
	int ksize;
	int64_t *n_reads;
};

struct maincount_bundle_t {
	struct dqueue_t *q;
	struct kmhash_t *h;
	struct kmhash_t *dict; // read-only small kmer dictionary
	int ksize;
};

void count_kmer_read(struct read_t *r, struct kmhash_t *h, int ksize)
{
	int i, last, c, len, lmc;
	char *seq;
	len = r->len;
	seq = r->seq;

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
		if (last >= ksize) {
			if (knum < krev)
				kmhash_inc_val(h, knum);
			else
				kmhash_inc_val(h, krev);
		}
	}
}

void *PE_minor_worker(void *data)
{
	struct precount_bundle_t *bundle = (struct precount_bundlet *)data;
	struct dqueue_t *q = bundle->q;
	struct kmhash_t *h = bundle->h;

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	int64_t n_reads;
	int64_t *gcnt_reads;
	gcnt_reads = bundle->n_reads;

	int ksize;
	ksize = bundle->ksize;

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

		n_reads = 0;
		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, buf1, &pos1) :
				get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, buf2, &pos2) :
				get_read_from_fa(&read2, buf2, &pos2);


			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file\n");

			++n_reads;
			count_kmer_read(&read1, h, ksize);
			count_kmer_read(&read2, h, ksize);

			if (rc1 == READ_END)
				break;
		}
		n_reads = __sync_add_and_fetch(gcnt_reads, n_reads);
		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

struct kmhash_t *count_kmer_minor(struct opt_count_t *opt, int kmer_size)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	int i;

	struct kmhash_t *kmer_hash;
	kmer_hash = init_kmhash((kmint_t)opt->hash_size, opt->n_threads);

	struct fastq_bundle_t *producer_bundles;
	producer_bundles = init_fastq_PE(opt);

	struct precount_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct precount_bundle_t));

	int64_t n_reads;
	n_reads = 0;

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].h = kmer_hash;
		worker_bundles[i].ksize = kmer_size;
		worker_bundles[i].n_reads = &n_reads;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, PE_minor_worker,
				worker_bundles + i);

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i]);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i])

	free_fastq_PE(producer_bundles, opt->n_files);
	free(worker_bundles);

	// get some stat
	return kmer_hash;
}

struct kmhash_t *count_kmer_master(struct opt_count_t *opt, struct kmhash_t *dict)
{
	struct kmhash_t *kmer_hash;
	kmer_hash = int_kmhash((kmint_t)opt->hash_size, opt->n_threads);
	return kmer_hash;
}

struct kmhash_t *count_kmer(struct opt_count_t *opt)
{
	struct kmhash_t *hs, *hm;
	hs = count_kmer_minor(opt, opt->kmer_slave);
	hm = count_kmer_master(opt, hs);
	kmhash_destroy(hs);
	return hm;
}

void kmer_test_process(struct opt_count_t *opt)
{
	struct kmhash_t *hs, *hl, *hm;
	hs = count_kmer_minor(opt, opt->kmer_slave);
	hl = count_kmer_minor(opt, opt->kmer_master);
	kmhash_destroy(hl);
	hm = count_kmer_master(opt, hs);
	kmhash_destroy(hs);
}

void kmer_fastq_count(struct opt_count_t *opt)
{
}

/**************************************************************/

// void *count_worker_fa(void *data)
// {
// 	struct worker_bundle_t *bundle = (struct worker_bundle_t *)data;
// 	// struct dqueue_t *q = bundle->q;
// 	// struct kmhash_t *h = bundle->h;
// 	// pthread_mutex_t *lock_hash = bundle->lock_hash;

// 	// struct kmthread_bundle_t kmhash_bundle;
// 	// kmhash_bundle.thread_no = bundle->thread_no;
// 	// kmhash_bundle.n_threads = bundle->n_threads;
// 	// kmhash_bundle.barrier = bundle->barrier_hash;
// 	// kmhash_init(h, bundle->init_hash_size, 0, &kmhash_bundle);

// 	struct read_t read;
// 	// struct read_t read1, read2;
// 	// struct pair_buffer_t *own_buf, *ext_buf;
// 	// own_buf = init_pair_buffer();

// 	// char *buf1, *buf2;
// 	// int pos1, pos2, rc1, rc2, input_format;
	
// 	kseq_t *seq;
// 	seq = bundle->kseq;

// 	struct worker_stat_t stat;
// 	memset(&stat, 0, sizeof(struct worker_stat_t));
// 	pthread_mutex_t *lock_input;
// 	lock_input = bundle->lock_input;

// 	while (1) {
// 		pthread_mutex_lock(lock_input);
// 		if (kseq_read(seq) >= 0) {
// 			__VERBOSE("Read sequence: %s\n", seq->name.s);
// 			// read.name = strdup(seq->name.s);
// 			read.len = seq->seq.l;
// 			read.seq = strdup(seq->seq.s);
// 		} else {
// 			read.seq = NULL;
// 			read.len = 0;
// 		}
// 		pthread_mutex_unlock(lock_input);
// 		if (read.seq == NULL)
// 			break;
// 		count_kmer_read(&read, bundle);
// 		free(read.seq);
// 	}
// 	// pthread_mutex_lock(bundle->lock_count);
// 	// bundle->global_stat->r1_sum += stat.r1_sum;
// 	// bundle->global_stat->r2_sum += stat.r2_sum;
// 	__sync_add_and_fetch(&(bundle->global_stat->nread), stat.nread);
// 	// bundle->global_stat->nread += stat.nread;
// 	// pthread_mutex_unlock(bundle->lock_count);

// 	pthread_exit(NULL);
// }

// struct kmhash_t *count_kmer_fasta(struct opt_count_t *opt)
// {
// 	pthread_attr_t attr;
// 	pthread_attr_init(&attr);
// 	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
// 	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

// 	struct worker_bundle_t *worker_bundles;
// 	pthread_t *worker_threads;

// 	struct worker_stat_t result;
// 	memset(&result, 0, sizeof(struct worker_stat_t));

// 	worker_bundles = malloc(opt->n_threads * sizeof(struct worker_bundle_t));
// 	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

// 	pthread_mutex_t lock_count, lock_input;
// 	pthread_mutex_init(&lock_count, NULL);
// 	pthread_mutex_init(&lock_input, NULL);

// 	struct kmhash_t *h;
// 	h = init_kmhash(opt->hash_size, opt->n_threads);

// 	kseq_t *kseq;
// 	gzFile fp = gzopen(opt->files_1[0], "r");
// 	kseq = kseq_init(fp);

// 	int i;
// 	for (i = 0; i < opt->n_threads; ++i) {
// 		// worker_bundles[i].q = q;
// 		worker_bundles[i].h = h;
// 		worker_bundles[i].n_threads = opt->n_threads;
// 		worker_bundles[i].kseq = kseq;
// 		// worker_bundles[i].thread_no = i;
// 		// worker_bundles[i].barrier_hash = &barrier_hash;
// 		worker_bundles[i].lock_count = &lock_count;
// 		worker_bundles[i].lock_input = &lock_input;
// 		worker_bundles[i].lock_hash = h->locks + i;
// 		worker_bundles[i].global_stat = &result;

// 		worker_bundles[i].ksize = opt->kmer_size;
// 		pthread_create(worker_threads + i, &attr, count_worker_fa, worker_bundles + i);
// 	}

// 	for (i = 0; i < opt->n_threads; ++i)
// 		pthread_join(worker_threads[i], NULL);

// 	kseq_destroy(kseq);
// 	gzclose(fp);

// 	__VERBOSE_LOG("Result", "Number of read                 : %20d\n", result.nread);
// 	__VERBOSE_LOG("Result", "Number of kmer                 : %20llu\n", (unsigned long long)h->n_items);
// 	return h;
// }

