#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "attribute.h"
#include "dqueue.h"
#include "fastq_producer.h"
#include "get_buffer.h"
#include "kmer_count.h"
#include "k63hash.h"
#include "kmhash.h"
#include "kseq.h"
#include "utils.h"
#include "verbose.h"

struct precount_bundle_t {
	struct dqueue_t *q;
	struct k63hash_t *h;
	int ksize;
	int64_t *n_reads;
	pthread_mutex_t *lock_hash;
};

struct maincount_bundle_t {
	struct dqueue_t *q;
	struct k63hash_t *h;
	struct k63hash_t *dict; // read-only small kmer dictionary
	int ksmall;
	int klarge;
	int64_t *n_reads;
	pthread_mutex_t *lock_hash;
};

#define __k63_lt(x, y) ((x).bin[1] < (y).bin[1] || ((x).bin[1] == (y).bin[1] && (x).bin[0] < (y).bin[0]))

#define __k63_lshift2(k) (((k).bin[1] = ((k).bin[1] << 2) | ((k).bin[0] >> 62)), \
				((k).bin[0] <<= 2))
#define __k63_rshift2(k) (((k).bin[0] = ((k).bin[0] >> 2) | (((k).bin[1] & 0x3ull) << 62)), \
				((k).bin[1] >>= 2))

#define __k63_lshift(k, l) (((k).bin[1] = ((k).bin[1] << (l)) | ((k).bin[0] >> (64 - (l)))), \
				((k).bin[0] <<= (l)))
#define __k63_rshift(k, l) (((k).bin[0] = ((k).bin[0] >> (l)) | (((k).bin[1] & ((1ull << (l)) - 1)) << (64 - (l))), \
				((k).bin[1] >>= (l))))

#define __k63_and(k, v) ((k).bin[0] &= (v).bin[0], (k).bin[1] &= (v).bin[1])

#define __reverse_bit_order64(x)					       \
(									       \
	(x) = (((x) & 0xffffffff00000000ull) >> 32) | (((x) & 0x00000000ffffffffull) << 32), \
	(x) = (((x) & 0xffff0000ffff0000ull) >> 16) | (((x) & 0x0000ffff0000ffffull) << 16), \
	(x) = (((x) & 0xff00ff00ff00ff00ull) >>  8) | (((x) & 0x00ff00ff00ff00ffull) <<  8), \
	(x) = (((x) & 0xf0f0f0f0f0f0f0f0ull) >>  4) | (((x) & 0x0f0f0f0f0f0f0f0full) <<  4), \
	(x) = (((x) & 0xccccccccccccccccull) >>  2) | (((x) & 0x3333333333333333ull) <<  2)  \
)

#define __k63_revc_num(y, x, l, mask)					       \
(									       \
	(x) = (y), __k63_lshift(x, 128 - ((l) << 1)),			       \
	__reverse_bit_order64((x).bin[0]), __reverse_bit_order64((x).bin[1]),  \
	(x).bin[0] ^= 0xffffffffffffffffull, (x).bin[0] &= (mask).bin[0],      \
	(x).bin[1] ^= 0xffffffffffffffffull, (x).bin[1] &= (mask).bin[1]       \
)

#define __get_revc_num(y, x, l, mask)					       \
(									       \
	(x) = (y) << (64 - ((l) << 1)),					       \
	__reverse_bit_order64(x), (x) ^= 0xffffffffffffffffull, (x) &= (mask)  \
)

/*
 * Move the kmer window along reads and add kmer to hash table
 */
static void count_lazy_from_read(struct read_t *r, struct k63hash_t *h,
					int ksize, pthread_mutex_t *lock_hash)
{
	int i, last, last_i, ci, ck, len, lmc, kedge;
	char *seq;
	len = r->len;
	seq = r->seq;

	k63key_t knum, krev, pknum, pkrev, kmask;
	kmask.bin[0] = (1ull << (ksize << 1)) - 1;
	kmask.bin[1] = (1ull << ((ksize << 1) - 64)) - 1;
	knum = krev = (k63key_t){0ull, 0ull};
	last = 0;
	lmc = (ksize - 1) << 1;
	kedge = ksize + 1;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		__k63_lshift2(knum); __k63_and(knum, kmask);
		__k63_rshift2(krev);
		if (ci < 4) {
			knum.bin[0] |= ci;
			/* please make sure 128 < lmc <= 64 */
			krev.bin[1] |= (uint64_t)(ci ^ 3) << (lmc - 64);
			++last;
		} else {
			last = 0;
		}
		if (last >= ksize) {
			if (__k63_lt(knum, krev))
				k63hash_put_adj(h, knum, lock_hash);
			else
				k63hash_put_adj(h, krev, lock_hash);
		}
		if (last >= kedge) {
			ck = nt4_table[(int)seq[i - ksize]] ^ 3;

			if (__k63_lt(pknum, pkrev))
				k63hash_add_edge(h, pknum, ci, lock_hash);
			else
				k63hash_add_edge(h, pkrev, ci + 4, lock_hash);

			if (__k63_lt(knum, krev))
				k63hash_add_edge(h, knum, ck + 4, lock_hash);
			else
				k63hash_add_edge(h, krev, ck, lock_hash);
		}
		pknum = knum;
		pkrev = krev;
	}
}

// static void count_small_from_read(struct read_t *r, struct k63hash_t *h, int ksize,
// 						pthread_mutex_t *lock_hash)
// {
// 	int i, last, c, len, lmc;
// 	char *seq;
// 	len = r->len;
// 	seq = r->seq;

// 	kmkey_t knum, krev, kmask;
// 	kmask = ((kmkey_t)1 << (ksize << 1)) - 1;
// 	knum = krev = 0;
// 	last = 0;
// 	lmc = (ksize - 1) << 1;
// 	for (i = 0; i < len; ++i) {
// 		c = nt4_table[(int)seq[i]];
// 		knum = (knum << 2) & kmask;
// 		krev = krev >> 2;
// 		if (c < 4) {
// 			knum |= c;
// 			krev |= (kmkey_t)(c ^ 3) << lmc;
// 			++last;
// 		} else {
// 			last = 0;
// 		}
// 		if (last >= ksize) {
// 			if (knum < krev)
// 				kmhash_put(h, knum, lock_hash);
// 			else
// 				kmhash_put(h, krev, lock_hash);
// 		}
// 	}
// }

// static void count_large_from_read(struct read_t *r, struct k63hash_t *h,
// 			int ksmall, int klarge, struct k63hash_t *dict,
// 			pthread_mutex_t *lock_hash)
// {
// 	kmkey_t knum_small, knum_large, krev_small, krev_large, pknum, pkrev;
// 	kmkey_t kmask_small, kmask_large;
// 	int i, last, cnt_small, pcnt_small, nk_edge, nk_large, ci, ck, len;
// 	int lmc_small, lmc_large, kedge, last_i;
// 	char *seq;
// 	kmint_t k;
// 	len = r->len;
// 	seq = r->seq;

// 	kmask_small = ((kmkey_t)1 << (ksmall << 1)) - 1;
// 	kmask_large = ((kmkey_t)1 << (klarge << 1)) - 1;
// 	kedge = klarge + 1;
// 	knum_small = knum_large = krev_small = krev_large = 0;
// 	last = cnt_small = pcnt_small = 0;
// 	nk_edge = kedge - ksmall + 1;
// 	nk_large = klarge - ksmall + 1;
// 	lmc_small = (ksmall - 1) << 1;
// 	lmc_large = (klarge - 1) << 1;
// 	for (i = 0; i < len; ++i) {
// 		ci = nt4_table[(int)seq[i]];
// 		knum_small = (knum_small << 2) & kmask_small;
// 		krev_small = krev_small >> 2;
// 		knum_large = (knum_large << 2) & kmask_large;
// 		krev_large = krev_large >> 2;
// 		if (ci < 4) {
// 			knum_small |= ci;
// 			knum_large |= ci;
// 			krev_small |= ((kmkey_t)(ci ^ 3) << lmc_small);
// 			krev_large |= ((kmkey_t)(ci ^ 3) << lmc_large);
// 			++last;
// 			if (last >= ksmall) {
// 				if (knum_small < krev_small)
// 					k = kmhash_get(dict, knum_small);
// 				else
// 					k = kmhash_get(dict, krev_small);
// 				if (k != KMHASH_MAX_SIZE)
// 					++cnt_small;
// 				else
// 					cnt_small = 0;
// 			} else {
// 				cnt_small = 0;
// 			}
// 		} else {
// 			last = 0;
// 			cnt_small = 0;
// 		}
// 		if (cnt_small >= nk_large) {
// 			if (knum_large < krev_large)
// 				kmhash_put_adj(h, knum_large, lock_hash);
// 			else
// 				kmhash_put_adj(h, krev_large, lock_hash);
// 		}
// 		if (last >= kedge && cnt_small >= nk_edge) {
// 			ck = nt4_table[(int)seq[i - klarge]] ^ 3;
// 			if (pknum < pkrev)
// 				kmhash_add_edge(h, pknum, ci, lock_hash);
// 			else
// 				kmhash_add_edge(h, pkrev, ci + 4, lock_hash);

// 			if (knum_large < krev_large)
// 				kmhash_add_edge(h, knum_large, ck + 4, lock_hash);
// 			else
// 				kmhash_add_edge(h, krev_large, ck, lock_hash);
// 		}
// 		pknum = knum_large;
// 		pkrev = krev_large;
// 	}
// }

/*
 * Iterate through all reads in buffer
 * Give read to counter process
 */
static void *PE_count_lazy_worker(void *data)
{
	struct maincount_bundle_t *bundle = (struct maincount_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct k63hash_t *h = bundle->h;

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	int64_t n_reads;
	int64_t *gcnt_reads;
	gcnt_reads = bundle->n_reads;

	int ksize;
	ksize = bundle->klarge;

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
			count_lazy_from_read(&read1, h, ksize, lock_hash);
			count_lazy_from_read(&read2, h, ksize, lock_hash);

			if (rc1 == READ_END)
				break;
		}
		n_reads = __sync_add_and_fetch(gcnt_reads, n_reads);
		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

// static void *PE_small_count_worker(void *data)
// {
// 	struct precount_bundle_t *bundle = (struct precount_bundle_t *)data;
// 	struct dqueue_t *q = bundle->q;
// 	struct k63hash_t *h = bundle->h;

// 	struct read_t read1, read2;
// 	struct pair_buffer_t *own_buf, *ext_buf;
// 	own_buf = init_pair_buffer();

// 	char *buf1, *buf2;
// 	int pos1, pos2, rc1, rc2, input_format;

// 	int64_t n_reads;
// 	int64_t *gcnt_reads;
// 	gcnt_reads = bundle->n_reads;

// 	int ksize;
// 	ksize = bundle->ksize;

// 	pthread_mutex_t *lock_hash;
// 	lock_hash = bundle->lock_hash;

// 	while (1) {
// 		ext_buf = d_dequeue_in(q);
// 		if (!ext_buf)
// 			break;
// 		d_enqueue_out(q, own_buf);
// 		own_buf = ext_buf;
// 		pos1 = pos2 = 0;
// 		buf1 = ext_buf->buf1;
// 		buf2 = ext_buf->buf2;
// 		input_format = ext_buf->input_format;

// 		n_reads = 0;
// 		while (1) {
// 			rc1 = input_format == TYPE_FASTQ ?
// 				get_read_from_fq(&read1, buf1, &pos1) :
// 				get_read_from_fa(&read1, buf1, &pos1);

// 			rc2 = input_format == TYPE_FASTQ ?
// 				get_read_from_fq(&read2, buf2, &pos2) :
// 				get_read_from_fa(&read2, buf2, &pos2);


// 			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
// 				__ERROR("\nWrong format file\n");

// 			++n_reads;
// 			count_small_from_read(&read1, h, ksize, lock_hash);
// 			count_small_from_read(&read2, h, ksize, lock_hash);

// 			if (rc1 == READ_END)
// 				break;
// 		}
// 		n_reads = __sync_add_and_fetch(gcnt_reads, n_reads);
// 		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
// 	}

// 	free_pair_buffer(own_buf);
// 	pthread_exit(NULL);
// }

// static void *PE_large_count_worker(void *data)
// {
// 	struct maincount_bundle_t *bundle = (struct maincount_bundle_t *)data;
// 	struct dqueue_t *q = bundle->q;
// 	struct k63hash_t *h = bundle->h;
// 	struct k63hash_t *dict = bundle->dict;

// 	struct read_t read1, read2;
// 	struct pair_buffer_t *own_buf, *ext_buf;
// 	own_buf = init_pair_buffer();

// 	char *buf1, *buf2;
// 	int pos1, pos2, rc1, rc2, input_format;

// 	int64_t n_reads;
// 	int64_t *gcnt_reads;
// 	gcnt_reads = bundle->n_reads;

// 	int ksmall, klarge;
// 	ksmall = bundle->ksmall;
// 	klarge = bundle->klarge;

// 	pthread_mutex_t *lock_hash;
// 	lock_hash = bundle->lock_hash;

// 	while (1) {
// 		ext_buf = d_dequeue_in(q);
// 		if (!ext_buf)
// 			break;
// 		d_enqueue_out(q, own_buf);
// 		own_buf = ext_buf;
// 		pos1 = pos2 = 0;
// 		buf1 = ext_buf->buf1;
// 		buf2 = ext_buf->buf2;
// 		input_format = ext_buf->input_format;

// 		n_reads = 0;
// 		while (1) {
// 			rc1 = input_format == TYPE_FASTQ ?
// 				get_read_from_fq(&read1, buf1, &pos1) :
// 				get_read_from_fa(&read1, buf1, &pos1);

// 			rc2 = input_format == TYPE_FASTQ ?
// 				get_read_from_fq(&read2, buf2, &pos2) :
// 				get_read_from_fa(&read2, buf2, &pos2);

// 			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
// 				__ERROR("\nWrong format file\n");

// 			++n_reads;
// 			count_large_from_read(&read1, h, ksmall, klarge,
// 							dict, lock_hash);
// 			count_large_from_read(&read2, h, ksmall, klarge,
// 							dict, lock_hash);

// 			if (rc1 == READ_END)
// 				break;
// 		}
// 		n_reads = __sync_add_and_fetch(gcnt_reads, n_reads);
// 		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
// 	}

// 	free_pair_buffer(own_buf);
// 	pthread_exit(NULL);
// }

/*
 * Start thread for counting processes
 */
static void count_kmer_lazy(struct opt_count_t *opt, struct k63hash_t *h, int ksize)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	int64_t n_reads;
	int i;

	struct producer_bundle_t *producer_bundles;
	producer_bundles = init_fastq_PE(opt);

	struct maincount_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct maincount_bundle_t));

	n_reads = 0;
	k63hash_init(h, (kmint_t)opt->hash_size - 1, opt->n_threads, 1);

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].h = h;
		worker_bundles[i].klarge = ksize;
		worker_bundles[i].n_reads = &n_reads;
		worker_bundles[i].lock_hash = h->locks + i;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, PE_count_lazy_worker,
				worker_bundles + i);

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	free_fastq_PE(producer_bundles, opt->n_files);
	free(worker_bundles);

	free(worker_threads);
	free(producer_threads);
}

// static void precount_small_kmer(struct opt_count_t *opt, struct k63hash_t *h,
// 							int kmer_size)
// {
// 	pthread_attr_t attr;
// 	pthread_attr_init(&attr);
// 	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
// 	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

// 	int64_t n_reads;
// 	int i;

// 	struct producer_bundle_t *producer_bundles;
// 	producer_bundles = init_fastq_PE(opt);

// 	struct precount_bundle_t *worker_bundles;
// 	worker_bundles = malloc(opt->n_threads * sizeof(struct precount_bundle_t));

// 	n_reads = 0;
// 	kmhash_init(h, (kmint_t)opt->hash_size - 1, opt->n_threads, 0);

// 	for (i = 0; i < opt->n_threads; ++i) {
// 		worker_bundles[i].q = producer_bundles->q;
// 		worker_bundles[i].h = h;
// 		worker_bundles[i].ksize = kmer_size;
// 		worker_bundles[i].n_reads = &n_reads;
// 		worker_bundles[i].lock_hash = h->locks + i;
// 	}

// 	pthread_t *producer_threads, *worker_threads;
// 	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
// 	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

// 	for (i = 0; i < opt->n_files; ++i)
// 		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
// 				producer_bundles + i);

// 	for (i = 0; i < opt->n_threads; ++i)
// 		pthread_create(worker_threads + i, &attr, PE_small_count_worker,
// 				worker_bundles + i);

// 	for (i = 0; i < opt->n_files; ++i)
// 		pthread_join(producer_threads[i], NULL);

// 	for (i = 0; i < opt->n_threads; ++i)
// 		pthread_join(worker_threads[i], NULL);

// 	free_fastq_PE(producer_bundles, opt->n_files);
// 	free(worker_bundles);

// 	free(producer_threads);
// 	free(worker_threads);
// }

// static void count_kmer_and_adj(struct opt_count_t *opt, struct k63hash_t *h,
// 				struct k63hash_t *dict, int ksmall, int klarge)
// {
// 	pthread_attr_t attr;
// 	pthread_attr_init(&attr);
// 	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
// 	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

// 	int64_t n_reads;
// 	int i;

// 	struct producer_bundle_t *producer_bundles;
// 	producer_bundles = init_fastq_PE(opt);

// 	struct maincount_bundle_t *worker_bundles;
// 	worker_bundles = malloc(opt->n_threads * sizeof(struct maincount_bundle_t));

// 	n_reads = 0;
// 	kmhash_init(h, (kmint_t)opt->hash_size - 1, opt->n_threads, 1);

// 	for (i = 0; i < opt->n_threads; ++i) {
// 		worker_bundles[i].q = producer_bundles->q;
// 		worker_bundles[i].h = h;
// 		worker_bundles[i].dict = dict;
// 		worker_bundles[i].klarge = klarge;
// 		worker_bundles[i].ksmall = ksmall;
// 		worker_bundles[i].n_reads = &n_reads;
// 		worker_bundles[i].lock_hash = h->locks + i;
// 	}

// 	pthread_t *producer_threads, *worker_threads;
// 	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
// 	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

// 	for (i = 0; i < opt->n_files; ++i)
// 		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
// 				producer_bundles + i);

// 	for (i = 0; i < opt->n_threads; ++i)
// 		pthread_create(worker_threads + i, &attr, PE_large_count_worker,
// 				worker_bundles + i);

// 	for (i = 0; i < opt->n_files; ++i)
// 		pthread_join(producer_threads[i], NULL);

// 	for (i = 0; i < opt->n_threads; ++i)
// 		pthread_join(worker_threads[i], NULL);

// 	free_fastq_PE(producer_bundles, opt->n_files);
// 	free(worker_bundles);

// 	free(worker_threads);
// 	free(producer_threads);
// }

void k63_correct_edge(struct k63hash_t *h, int ksize)
{
	k63key_t knum, krev, nknum, nkrev, kmask;
	kmint_t i, k;
	int c, lmc;
	kmask.bin[0] = (1ull << (ksize << 1)) - 1;
	kmask.bin[1] = (1ull << ((ksize << 1) - 64)) - 1;
	knum = krev = (k63key_t){0ull, 0ull};
	lmc = (ksize - 1) << 1;

	for (i = 0; i < h->size; ++i) {
		if (h->flag[i] == KMFLAG_EMPTY)
			continue;
		knum = h->keys[i];
		__k63_revc_num(knum, krev, ksize, kmask);
		for (c = 0; c < 4; ++c) {
			if ((h->adjs[i] >> c) & 1) {
				nknum = knum; __k63_lshift2(nknum);
				__k63_and(nknum, kmask); nknum.bin[0] |= c;
				nkrev = krev; __k63_rshift2(nkrev);
				nkrev.bin[1] |= (uint64_t)(c ^ 3) << (lmc - 64);
				if (__k63_lt(nknum, nkrev))
					k = k63hash_get(h, nknum);
				else
					k = k63hash_get(h, nkrev);
				if (k == KMHASH_MAX_SIZE)
					h->adjs[i] &= ~((uint8_t)1 << c);
			}
			if ((h->adjs[i] >> (c + 4)) & 1) {
				nknum = krev; __k63_lshift2(nknum);
				__k63_and(nknum, kmask); nknum.bin[0] |= c;
				nkrev = knum; __k63_rshift2(nkrev);
				nkrev.bin[1] |= (uint64_t)(c ^ 3) << (lmc - 64);
				if (__k63_lt(nknum, nkrev))
					k = k63hash_get(h, nknum);
				else
					k = k63hash_get(h, nkrev);
				if (k == KMHASH_MAX_SIZE)
					h->adjs[i] &= ~((uint8_t)1 << (c + 4));
			}
		}
	}
}

// static inline void dump_seq(const char *label, kmkey_t kmer, int size)
// {
// 	char *seq = malloc(size + 1);
// 	int i;
// 	for (i = 0; i < size; ++i) {
// 		seq[size - i - 1] = nt4_char[kmer & 3];
// 		kmer >>= 2;
// 	}
// 	seq[size] = '\0';
// 	fprintf(stderr, "%s%s\n", label, seq);
// 	free(seq);
// }

// void check_edge(struct k63hash_t *h, int ksize)
// {
// 	kmkey_t *keys;
// 	uint8_t *adjs;
// 	kmint_t i, k;
// 	kmkey_t kmer, rkmer, nkmer, rnkmer, mask;
// 	int c, lmc, rc;

// 	keys = h->keys;
// 	adjs = h->adjs;
// 	mask = ((kmkey_t)1 << (ksize << 1)) - 1;
// 	lmc = (ksize - 1) << 1;
// 	for (i = 0; i < h->size; ++i) {
// 		if (keys[i] == TOMB_STONE)
// 			continue;
// 		kmer = keys[i];
// 		__get_revc_num(kmer, rkmer, ksize, mask);
// 		for (c = 0; c < 4; ++c) {
// 			if ((adjs[i] >> c) & 1) {
// 				nkmer = ((kmer << 2) & mask) | c;
// 				rnkmer = (rkmer >> 2) | ((kmkey_t)(c ^ 3) << lmc);
// 				if (nkmer < rnkmer)
// 					k = kmhash_get(h, nkmer);
// 				else
// 					k = kmhash_get(h, rnkmer);
// 				assert(k != KMHASH_MAX_SIZE);
// 				rc = kmer >> lmc;
// 				assert(rc >= 0 && rc < 4);
// 				if (nkmer < rnkmer)
// 					assert((adjs[k] >> ((rc ^ 3) + 4)) & 1);
// 				else
// 					assert((adjs[k] >> (rc ^ 3)) & 1);
// 			}

// 			if ((adjs[i] >> (c + 4)) & 1) {
// 				nkmer = ((rkmer << 2) & mask) | c;
// 				rnkmer = (kmer >> 2) | ((kmkey_t)(c ^ 3) << lmc);
// 				if (nkmer < rnkmer)
// 					k = kmhash_get(h, nkmer);
// 				else
// 					k = kmhash_get(h, rnkmer);
// 				assert(k != KMHASH_MAX_SIZE);
// 				rc = rkmer >> lmc;
// 				assert(rc >= 0 && rc < 4);
// 				if (nkmer < rnkmer)
// 					assert((adjs[k] >> ((rc ^ 3) + 4)) & 1);
// 				else
// 					assert((adjs[k] >> (rc ^ 3)) & 1);
// 			}
// 		}
// 	}

// 	__VERBOSE("Check edge kmhash done\n");
// }

// void recount_edge(struct k63hash_t *h)
// {
// 	kmint_t ret, i;
// 	int c;
// 	ret = 0;
// 	for (i = 0; i < h->size; ++i) {
// 		if (h->keys[i] == TOMB_STONE)
// 			continue;
// 		for (c = 0; c < 8; ++c)
// 			ret += ((h->adjs[i] >> c) & 1);
// 	}
// 	__VERBOSE("Edge count: %llu\n", (long long unsigned)ret);
// 	kmkey_t sum;
// 	sum = 0;
// 	for (i = 0; i < h->size; ++i) {
// 		if (h->keys[i] == TOMB_STONE)
// 			continue;
// 		sum += (h->keys[i] ^ (kmkey_t)h->adjs[i]);
// 	}
// 	__VERBOSE("hash sum: %llu\n", (long long unsigned)sum);
// }

void build_k63_table_lazy(struct opt_count_t *opt, struct k63hash_t *h)
{
	__VERBOSE("\nCounting %d-mer\n", opt->kmer_master);
	count_kmer_lazy(opt, h, opt->kmer_master);
	__VERBOSE("\n");
	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer: %llu\n",
			opt->kmer_master, (long long unsigned)h->n_item);
	// check_edge(h, opt->kmer_master);
	// recount_edge(h);

	/* Filter singleton kmer */
	k63hash_filter(h, 1);
	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %llu\n",
			opt->kmer_master, (long long unsigned)h->n_item);

	/* Correct edges */
	k63_correct_edge(h, opt->kmer_master);
	// recount_edge(h);
	// check_edge(h, opt->kmer_master);
}

// void build_kmer_table(struct opt_count_t *opt, struct k63hash_t *h)
// {
// 	struct k63hash_t *hs;
// 	hs = calloc(1, sizeof(struct k63hash_t));
// 	__VERBOSE("\nCounting %d-mer\n", opt->kmer_slave);
// 	precount_small_kmer(opt, h, opt->kmer_slave);
// 	__VERBOSE("\n");
// 	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer: %llu\n", opt->kmer_slave,
// 					(long long unsigned)hs->n_item);

// 	__VERBOSE("Filtering %d-mer\n", opt->kmer_slave);
// 	kmhash_filter(hs, 0);
// 	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %llu\n",
// 			opt->kmer_slave, (long long unsigned)hs->n_item);

// 	__VERBOSE("\nCounting %d-mer\n", opt->kmer_master);
// 	count_kmer_and_adj(opt, h, hs, opt->kmer_slave, opt->kmer_master);
// 	__VERBOSE("\n");
// 	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer : %llu\n",
// 			opt->kmer_master, (long long unsigned)h->n_item);

// 	kmhash_destroy(hs);
// 	free(hs);
// 	__VERBOSE("Filtering %d-mer\n", opt->kmer_master);
// 	kmhash_filter(h, 1);
// 	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %llu\n",
// 			opt->kmer_master, (long long unsigned)h->n_item);

// 	correct_edge(h, opt->kmer_master);
// }

// void test_relocate(struct k63hash_t *h)
// {
// 	kmint_t i, cnt, k;
// 	kmkey_t hash;
// 	cnt = 0;
// 	hash = 0;
// 	for (i = 0; i < h->size; ++i) {
// 		if (h->keys[i] != TOMB_STONE) {
// 			k = kmhash_get(h, h->keys[i]);
// 			if (k == KMHASH_MAX_SIZE)
// 				__ERROR("test relocate failed");
// 			++cnt;
// 			hash += h->keys[i];
// 		}
// 	}
// 	__VERBOSE("Recount kmer = %llu\n", (long long unsigned)cnt);
// 	__VERBOSE("Hash = %llu\n", (long long unsigned)hash);
// }

void k63_test_process(struct opt_count_t *opt)
{
	struct k63hash_t *hm;
	hm = calloc(1, sizeof(struct k63hash_t));
	__VERBOSE("Lazy count %d-mer\n", opt->kmer_master);
	build_k63_table_lazy(opt, hm);

	// struct k63hash_t *hs, *hl, *hm;
	// __VERBOSE("\nCounting %d-mer\n", opt->kmer_slave);
	// hs = calloc(1, sizeof(struct k63hash_t));
	// precount_small_kmer(opt, hs, opt->kmer_slave);

	// __VERBOSE("\n");
	// __VERBOSE("Number of %d-mer: %llu\n", opt->kmer_slave,
	// 				(long long unsigned)hs->n_item);

	// k63hash_filter(hs, 0);
	// __VERBOSE("Number of non-singleton %d-mer: %llu\n", opt->kmer_slave,
	// 				(long long unsigned)hs->n_item);

	// test_relocate(hs);
	// // fprintf(stderr, "pre-count = %llu\n", cnt);

	// // for (k = 0; k < cnt; ++k) {
	// // 	ret = kmphash_get(hs, tmp[k]);
	// // 	if (ret == KMHASH_MAX_SIZE)
	// // 		fprintf(stderr, "Fail test\n");
	// // }
	// // fprintf(stderr, "Pass test\n");

	// /* Lazy counting kmer */
	// __VERBOSE("\nLazy count %d-mer\n", opt->kmer_master);
	// hl = calloc(1, sizeof(struct k63hash_t));
	// build_kmer_table_lazy(opt, hl);

	// recount_edge(hl);

	// __VERBOSE("\nCount %d-mer from %d-mer\n", opt->kmer_master,
	// 						opt->kmer_slave);
	// hm = calloc(1, sizeof(struct k63hash_t));
	// count_kmer_and_adj(opt, hm, hs, opt->kmer_slave, opt->kmer_master);
	// __VERBOSE("\n");
	// __VERBOSE("Number of %d-mer (2-step counting): %llu\n",
	// 		opt->kmer_master, (long long unsigned)hm->n_item);

	// kmhash_filter(hm, 1);
	// __VERBOSE("Number of non-singleton %d-mer (2-step counting): %llu\n",
	// 		opt->kmer_master, (long long unsigned)hm->n_item);

	// test_relocate(hm);

	// correct_edge(hm, opt->kmer_master);

	// recount_edge(hm);
}

