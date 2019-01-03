#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "attribute.h"
#include "dqueue.h"
#include "fastq_producer.h"
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
	pthread_mutex_t *lock_hash;
};

struct maincount_bundle_t {
	struct dqueue_t *q;
	struct kmhash_t *h;
	struct kmhash_t *dict; // read-only small kmer dictionary
	int ksmall;
	int klarge;
	int64_t *n_reads;
	pthread_mutex_t *lock_hash;
};

void count_kmer_read(struct read_t *r, struct kmhash_t *h, int ksize, pthread_mutex_t *lock_hash)
{
	int i, last, c, len, lmc;
	char *seq;
	len = r->len;
	seq = r->seq;

	kmkey_t knum, krev, kmask;
	kmask = ((kmkey_t)1 << (ksize << 1)) - 1;
	knum = krev = 0;
	last = 0;
	lmc = (ksize - 1) << 1;
	for (i = 0; i < len; ++i) {
		c = nt4_table[(int)seq[i]];
		knum = (knum << 2) & kmask;
		krev = krev >> 2;
		if (c < 4) {
			knum |= c;
			krev |= (kmkey_t)(c ^ 3) << lmc;
			++last;
		} else {
			last = 0;
		}
		if (last >= ksize) {
			if (knum < krev)
				// kmhash_inc_val(h, knum, lock_hash);
				kmhash_put(h, knum, lock_hash);
			else
				// kmhash_inc_val(h, krev, lock_hash);
				kmhash_put(h, krev, lock_hash);
		}
	}
}

void count_kmer_filter(struct read_t *r, struct kmhash_t *h,
			int ksmall, int klarge, struct kmhash_t *dict,
			pthread_mutex_t *lock_hash)
{
	int i, last, cnt_small, n_small, c, len, lmc_small, lmc_large;
	char *seq;
	kmint_t k;
	len = r->len;
	seq = r->seq;

	kmkey_t knum_small, knum_large, krev_small, krev_large, kmask_small, kmask_large;
	kmask_small = ((kmkey_t)1 << (ksmall << 1)) - 1;
	kmask_large = ((kmkey_t)1 << (klarge << 1)) - 1;
	knum_small = knum_large = krev_small = krev_large = 0;
	last = cnt_small = 0;
	n_small = klarge - ksmall + 1;
	lmc_small = (ksmall - 1) << 1;
	lmc_large = (klarge - 1) << 1;
	for (i = 0; i < len; ++i) {
		c = nt4_table[(int)seq[i]];
		knum_small = (knum_small << 2) & kmask_small;
		krev_small = krev_small >> 2;
		knum_large = (knum_large << 2) & kmask_large;
		krev_large = krev_large >> 2;
		if (c < 4) {
			knum_small |= c;
			knum_large |= c;
			krev_small |= (kmkey_t)(c ^ 3) << lmc_small;
			krev_large |= (kmkey_t)(c ^ 3) << lmc_large;
			++last;
			if (knum_small < krev_small)
				k = kmphash_get(dict, knum_small);
			else
				k = kmphash_get(dict, krev_small);
			if (k != KMHASH_MAX_SIZE)
				++cnt_small;
			else
				cnt_small = 0;
		} else {
			last = 0;
			cnt_small = 0;
		}
		if (last >= klarge && cnt_small >= n_small) {
			if (knum_large < krev_large)
				// kmhash_inc_val(h, knum_large, lock_hash);
				kmhash_put(h, knum_large, lock_hash);
			else
				// kmhash_inc_val(h, krev_large, lock_hash);
				kmhash_put(h, krev_large, lock_hash);
		}
	}
}

void *PE_slave_worker(void *data)
{
	struct precount_bundle_t *bundle = (struct precount_bundle_t *)data;
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

	pthread_mutex_t *lock_hash;
	lock_hash = bundle->lock_hash;

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
			count_kmer_read(&read1, h, ksize, lock_hash);
			count_kmer_read(&read2, h, ksize, lock_hash);

			if (rc1 == READ_END)
				break;
		}
		n_reads = __sync_add_and_fetch(gcnt_reads, n_reads);
		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

void *PE_master_worker(void *data)
{
	struct maincount_bundle_t *bundle = (struct maincount_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct kmhash_t *h = bundle->h;
	struct kmhash_t *dict = bundle->dict;

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	int64_t n_reads;
	int64_t *gcnt_reads;
	gcnt_reads = bundle->n_reads;

	int ksmall, klarge;
	ksmall = bundle->ksmall;
	klarge = bundle->klarge;

	pthread_mutex_t *lock_hash;
	lock_hash = bundle->lock_hash;

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
			count_kmer_filter(&read1, h, ksmall, klarge, dict, lock_hash);
			count_kmer_filter(&read2, h, ksmall, klarge, dict, lock_hash);

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
	kmer_hash = init_filter_kmhash((kmint_t)opt->hash_size - 1, opt->n_threads);

	struct producer_bundle_t *producer_bundles;
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
		worker_bundles[i].lock_hash = kmer_hash->locks + i;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, PE_slave_worker,
				worker_bundles + i);

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	free_fastq_PE(producer_bundles, opt->n_files);
	free(worker_bundles);

	free(producer_threads);
	free(worker_threads);

	// get some stat
	return kmer_hash;
}

struct kmhash_t *count_kmer_master(struct opt_count_t *opt, struct kmhash_t *dict)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	int i;

	struct kmhash_t *kmer_hash;
	kmer_hash = init_filter_kmhash((kmint_t)opt->hash_size - 1, opt->n_threads);

	struct producer_bundle_t *producer_bundles;
	producer_bundles = init_fastq_PE(opt);

	struct maincount_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct maincount_bundle_t));

	int64_t n_reads;
	n_reads = 0;

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].h = kmer_hash;
		worker_bundles[i].dict = dict;
		worker_bundles[i].klarge = opt->kmer_master;
		worker_bundles[i].ksmall = opt->kmer_slave;
		worker_bundles[i].n_reads = &n_reads;
		worker_bundles[i].lock_hash = kmer_hash->locks + i;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, PE_master_worker,
				worker_bundles + i);

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	free_fastq_PE(producer_bundles, opt->n_files);
	free(worker_bundles);

	free(worker_threads);
	free(producer_threads);

	return kmer_hash;
}

struct kmhash_t *count_kmer(struct opt_count_t *opt)
{
	struct kmhash_t *hs, *hm;
	__VERBOSE("\nCounting %d-mer\n", opt->kmer_slave);
	hs = count_kmer_minor(opt, opt->kmer_slave);
	__VERBOSE("\n");
	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer: %llu\n", opt->kmer_slave,
		(long long unsigned)hs->n_items);

	__VERBOSE("Filtering %d-mer\n", opt->kmer_slave);
	kmhash_filter_singleton(hs, 0);
	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %llu\n",
			opt->kmer_slave, (long long unsigned)hs->n_items);
	// uint64_t cnt_good;
	// kmint_t i;
	// cnt_good = 0;
	// for (i = 0; i < hs->size; ++i) {
	// 	if (hs->keys[i] != TOMB_STONE && hs->vals[i] > opt->filter_thres)
	// 		++cnt_good;
	// 	else
	// 		hs->vals[i] = 0;
	// }
	// __VERBOSE_LOG("KMER COUNT", "Number of %d-mer with count greater than (%d): %llu\n",
	// 	opt->kmer_slave, opt->filter_thres, (long long unsigned)cnt_good);
	__VERBOSE("\nCounting %d-mer (early singleton-filtered with %d-mer)\n", opt->kmer_master, opt->kmer_slave);
	hm = count_kmer_master(opt, hs);
	__VERBOSE("\n");
	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer (early singleton-filtered with %d-mer): %llu\n", opt->kmer_master, opt->kmer_slave, (long long unsigned)hm->n_items);

	kmhash_destroy(hs);
	kmhash_filter_singleton(hm, 0);
	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %llu\n", opt->kmer_master, (long long unsigned)hm->n_items);
	return hm;
}

void kmer_test_process(struct opt_count_t *opt)
{
	struct kmhash_t *hs, *hl, *hm;
	__VERBOSE("\nCounting %d-mer\n", opt->kmer_slave);
	hs = count_kmer_minor(opt, opt->kmer_slave);

	__VERBOSE("\n");
	__VERBOSE("Number of %d-mer: %llu\n", opt->kmer_slave, (long long unsigned)hs->n_items);
	// __VERBOSE("Number of non-singleton %d-mer: %llu\n", opt->kmer_slave, (long long unsigned)cnt_slave);

	// kmint_t cnt, k, ret;
	// cnt = 0;
	// for (k = 0; k < hs->size; ++k)
	// 	if ((hs->vals[k >> 5] >> (k & 31)) & 1)
	// 		++cnt;
	// kmkey_t *tmp;
	// tmp = malloc(cnt * sizeof(kmkey_t));
	// cnt = 0;
	// for (k = 0; k < hs->size; ++k)
	// 	if ((hs->vals[k >> 5] >> (k & 31)) & 1)
	// 		tmp[cnt++] = hs->keys[k];

	kmhash_filter_singleton(hs, 0);
	__VERBOSE("Number of non-singleton %d-mer: %llu\n", opt->kmer_slave, (long long unsigned)hs->n_items);

	// test 
	
	// fprintf(stderr, "pre-count = %llu\n", cnt);

	// for (k = 0; k < cnt; ++k) {
	// 	ret = kmphash_get(hs, tmp[k]);
	// 	if (ret == KMHASH_MAX_SIZE)
	// 		fprintf(stderr, "Fail test\n");
	// }
	// fprintf(stderr, "Pass test\n");

	__VERBOSE("\nCounting %d-mer\n", opt->kmer_master);
	hl = count_kmer_minor(opt, opt->kmer_master);
	// uint64_t cnt_slave, cnt_master;
	// kmint_t k;
	// cnt_slave = cnt_master = 0;
	// for (k = 0; k < hs->size; ++k)
	// 	if (hs->keys[k] != TOMB_STONE && hs->vals[k] > opt->filter_thres)
	// 		++cnt_slave;
	// 	else
	// 		hs->vals[k] = 0;
	// for (k = 0; k < hl->size; ++k)
	// 	if (hl->keys[k] != TOMB_STONE && hl->vals[k] > opt->filter_thres)
	// 		++cnt_master;

	__VERBOSE("\n");
	__VERBOSE("Number of %d-mer: %llu\n", opt->kmer_master, (long long unsigned)hl->n_items);
	// __VERBOSE("Number of non-singleton %d-mer: %llu\n", opt->kmer_master, (long long unsigned)cnt_master);

	kmhash_filter_singleton(hl, 0);
	__VERBOSE("Number of non-singleton %d-mer: %llu\n", opt->kmer_master, (long long unsigned)hl->n_items);

	// kmhash_destroy(hl);
	hm = count_kmer_master(opt, hs);
	// uint64_t cnt_master_2;
	// cnt_master_2 = 0;
	// for (k = 0; k < hm->size; ++k)
	// 	if (hm->keys[k] != TOMB_STONE && hm->vals[k] > opt->filter_thres)
	// 		++cnt_master_2;
	// 	else
	// 		hm->vals[k] = 0;
	__VERBOSE("\n");
	__VERBOSE("Number of %d-mer (2-step counting): %llu\n", opt->kmer_master, (long long unsigned)hm->n_items);
	// __VERBOSE("Number of non-singleton %d-mer (2-step counting): %llu\n", opt->kmer_master, (long long unsigned)cnt_master_2);

	kmhash_filter_singleton(hm, 0);
	__VERBOSE("Number of non-singleton %d-mer (2-step counting): %llu\n", opt->kmer_master, (long long unsigned)hm->n_items);

	// for (k = 0; k < hl->size; ++k) {
	// 	if (hl->keys[k] != TOMB_STONE) {
	// 		ret = kmphash_get(hm, hl->keys[k]);
	// 		assert(ret != KMHASH_MAX_SIZE && "Fail test #1");
	// 	}
	// }

	// for (k = 0; k < hm->size; ++k) {
	// 	if (hm->keys[k] != TOMB_STONE) {
	// 		ret = kmphash_get(hl, hm->keys[k]);
	// 		assert(ret != KMHASH_MAX_SIZE && "Fail test #2");
	// 	}
	// }

	kmhash_destroy(hs);
}

// struct kmhash_t *count_kmer(struct opt_count_t *opt)
// {
// 	struct kmhash_t *hs, *hm;
// 	__VERBOSE("\nCounting %d-mer\n", opt->kmer_slave);
// 	hs = count_kmer_minor(opt, opt->kmer_slave);
// 	__VERBOSE("\n");
// 	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer: %llu\n", opt->kmer_slave,
// 		(long long unsigned)hs->n_items);
// 	uint64_t cnt_good;
// 	kmint_t i;
// 	cnt_good = 0;
// 	for (i = 0; i < hs->size; ++i) {
// 		if (hs->keys[i] != TOMB_STONE && hs->vals[i] > opt->filter_thres)
// 			++cnt_good;
// 		else
// 			hs->vals[i] = 0;
// 	}
// 	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer with count greater than (%d): %llu\n",
// 		opt->kmer_slave, opt->filter_thres, (long long unsigned)cnt_good);
// 	__VERBOSE("\nCounting %d-mer (early filtered with %d-mer)\n", opt->kmer_master, opt->kmer_slave);
// 	hm = count_kmer_master(opt, hs);
// 	kmhash_destroy(hs);
// 	return hm;
// }

// void kmer_test_process(struct opt_count_t *opt)
// {
// 	struct kmhash_t *hs, *hl, *hm;
// 	__VERBOSE("\nCounting %d-mer\n", opt->kmer_slave);
// 	hs = count_kmer_minor(opt, opt->kmer_slave);
// 	__VERBOSE("\nCounting %d-mer\n", opt->kmer_master);
// 	hl = count_kmer_minor(opt, opt->kmer_master);
// 	uint64_t cnt_slave, cnt_master;
// 	kmint_t k;
// 	cnt_slave = cnt_master = 0;
// 	for (k = 0; k < hs->size; ++k)
// 		if (hs->keys[k] != TOMB_STONE && hs->vals[k] > opt->filter_thres)
// 			++cnt_slave;
// 		else
// 			hs->vals[k] = 0;
// 	for (k = 0; k < hl->size; ++k)
// 		if (hl->keys[k] != TOMB_STONE && hl->vals[k] > opt->filter_thres)
// 			++cnt_master;
// 	__VERBOSE("\n");
// 	__VERBOSE("Number of %d-mer: %llu\n", opt->kmer_slave, (long long unsigned)hs->n_items);
// 	__VERBOSE("Number of non-singleton %d-mer: %llu\n", opt->kmer_slave, (long long unsigned)cnt_slave);

// 	__VERBOSE("Number of %d-mer: %llu\n", opt->kmer_master, (long long unsigned)hl->n_items);
// 	__VERBOSE("Number of non-singleton %d-mer: %llu\n", opt->kmer_master, (long long unsigned)cnt_master);

// 	kmhash_destroy(hl);
// 	hm = count_kmer_master(opt, hs);
// 	uint64_t cnt_master_2;
// 	cnt_master_2 = 0;
// 	for (k = 0; k < hm->size; ++k)
// 		if (hm->keys[k] != TOMB_STONE && hm->vals[k] > opt->filter_thres)
// 			++cnt_master_2;
// 		else
// 			hm->vals[k] = 0;
// 	__VERBOSE("\n");
// 	__VERBOSE("Number of %d-mer (2-step counting): %llu\n", opt->kmer_master, (long long unsigned)hm->n_items);
// 	__VERBOSE("Number of non-singleton %d-mer (2-step counting): %llu\n", opt->kmer_master, (long long unsigned)cnt_master_2);

// 	kmhash_destroy(hs);
// }

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

