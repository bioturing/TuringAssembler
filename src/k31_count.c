#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "attribute.h"
#include "dqueue.h"
#include "fastq_producer.h"
#include "get_buffer.h"
#include "k31hash.h"
#include "k31_count.h"
#include "kseq.h"
#include "test.h"
#include "utils.h"
#include "verbose.h"

struct precount_bundle_t {
	struct dqueue_t *q;
	struct k31hash_t *h;
	int ksize;
	int64_t *n_reads;
	pthread_mutex_t *lock_hash;
};

struct maincount_bundle_t {
	struct dqueue_t *q;
	struct k31hash_t *h;
	struct k31hash_t *dict; // read-only small kmer dictionary
	int ksmall;
	int klarge;
	int64_t *n_reads;
	pthread_mutex_t *lock_hash;
};

#define __reverse_bit_order64(x)					       \
(									       \
	(x) = (((x) & 0xffffffff00000000ull) >> 32) | (((x) & 0x00000000ffffffffull) << 32), \
	(x) = (((x) & 0xffff0000ffff0000ull) >> 16) | (((x) & 0x0000ffff0000ffffull) << 16), \
	(x) = (((x) & 0xff00ff00ff00ff00ull) >>  8) | (((x) & 0x00ff00ff00ff00ffull) <<  8), \
	(x) = (((x) & 0xf0f0f0f0f0f0f0f0ull) >>  4) | (((x) & 0x0f0f0f0f0f0f0f0full) <<  4), \
	(x) = (((x) & 0xccccccccccccccccull) >>  2) | (((x) & 0x3333333333333333ull) <<  2)  \
)

#define __k31_revc_num(y, x, l, mask)					       \
(									       \
	(x) = (y) << (64 - ((l) << 1)),					       \
	__reverse_bit_order64(x), (x) ^= 0xffffffffffffffffull, (x) &= (mask)  \
)

/*
 * Move the kmer window along reads and add kmer to hash table
 */
static void count_lazy_from_read(struct read_t *r, struct k31hash_t *h,
					int ksize, pthread_mutex_t *lock_hash)
{
	int i, last, ci, ck, len, lmc, kedge;
	char *seq;
	len = r->len;
	seq = r->seq;

	k31key_t knum, krev, pknum, pkrev, kmask;
	kmask = ((k31key_t)1 << (ksize << 1)) - 1;
	knum = krev = 0;
	last = 0;
	lmc = (ksize - 1) << 1;
	kedge = ksize + 1;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		knum = (knum << 2) & kmask;
		krev >>= 2;
		if (ci < 4) {
			knum |= ci;
			krev |= (k31key_t)(ci ^ 3) << lmc;
			++last;
		} else {
			last = 0;
		}
		if (last >= ksize) {
			if (knum < krev)
				k31hash_put_adj(h, knum, lock_hash);
			else
				k31hash_put_adj(h, krev, lock_hash);
		}
		if (last >= kedge) {
			ck = nt4_table[(int)seq[i - ksize]] ^ 3;

			if (pknum < pkrev)
				k31hash_add_edge(h, pknum, ci, lock_hash);
			else
				k31hash_add_edge(h, pkrev, ci + 4, lock_hash);

			if (knum < krev)
				k31hash_add_edge(h, knum, ck + 4, lock_hash);
			else
				k31hash_add_edge(h, krev, ck, lock_hash);
		}
		pknum = knum;
		pkrev = krev;
	}
}

/*
 * Iterate through all reads in buffer
 * Give read to counter process
 */
static void *PE_count_lazy_worker(void *data)
{
	struct maincount_bundle_t *bundle = (struct maincount_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct k31hash_t *h = bundle->h;

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

/*
 * Start thread for counting processes
 */
static void count_kmer_lazy(struct opt_count_t *opt, struct k31hash_t *h, int ksize)
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
	k31hash_init(h, (kmint_t)opt->hash_size - 1, opt->n_threads, 1);

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

void k31_correct_edge(struct k31hash_t *h, int ksize)
{
	k31key_t knum, krev, nknum, nkrev, kmask;
	kmint_t i, k;
	int c, lmc;
	kmask = ((k31key_t)1 << (ksize << 1)) - 1;
	lmc = (ksize - 1) << 1;

	for (i = 0; i < h->size; ++i) {
		if (h->keys[i] == K31_NULL)
			continue;
		knum = h->keys[i];
		__k31_revc_num(knum, krev, ksize, kmask);
		for (c = 0; c < 4; ++c) {
			if ((h->adjs[i] >> c) & 1) {
				nknum = ((knum << 2) & kmask) | c;
				nkrev = (krev >> 2) | ((k31key_t)(c ^ 3) << lmc);
				if (nknum < nkrev)
					k = k31hash_get(h, nknum);
				else
					k = k31hash_get(h, nkrev);
				if (k == KMHASH_MAX_SIZE)
					h->adjs[i] &= ~((uint8_t)1 << c);
			}
			if ((h->adjs[i] >> (c + 4)) & 1) {
				nknum = ((krev << 2) & kmask) | c;
				nkrev = (knum >> 2) | ((k31key_t)(c ^ 3) << lmc);
				if (nknum < nkrev)
					k = k31hash_get(h, nknum);
				else
					k = k31hash_get(h, nkrev);
				if (k == KMHASH_MAX_SIZE)
					h->adjs[i] &= ~((uint8_t)1 << (c + 4));
			}
		}
	}
}

static void k31_test_hash(struct k31hash_t *h)
{
	kmint_t cnt, i, k;
	cnt = 0;
	for (i = 0; i < h->size; ++i) {
		if (h->keys[i] != K31_NULL) {
			++cnt;
			k = k31hash_get(h, h->keys[i]);
			if (k == KMHASH_END(h)) {
				fprintf(stderr, "i = %lu; key = %lu\n", i, h->keys[i]);
				assert(0);
			}
		}
	}
	assert(cnt == h->n_item);
	fprintf(stderr, "test hash relocate done\n");
}

static void k31_check_sum(struct k31hash_t *h)
{
	kmint_t i;
	uint64_t sum = 0;
	for (i = 0; i < h->size; ++i) {
		if (h->keys[i] != K31_NULL) {
			sum += h->keys[i];
		}
	}
	fprintf(stderr, "check sum = %lu\n", sum);
}

static void k31_check_edge(struct k31hash_t *h, int ksize)
{
	kmint_t i, k;
	k31key_t knum, krev, nknum, nkrev, kmask;
	int c, lmc, rc;

	kmask = ((k31key_t)1 << (ksize << 1)) - 1;
	lmc = (ksize - 1) << 1;
	kmint_t sum_deg = 0;
	for (i = 0; i < h->size; ++i) {
		if (h->keys[i] == K31_NULL)
			continue;
		knum = h->keys[i];
		__k31_revc_num(knum, krev, ksize, kmask);
		for (c = 0; c < 4; ++c) {
			if ((h->adjs[i] >> c) & 1) {
				nknum = ((knum << 2) & kmask) | c;
				nkrev = (krev >> 2) | ((k31key_t)(c ^ 3) << lmc);
				if (nknum < nkrev)
					k = k31hash_get(h, nknum);
				else
					k = k31hash_get(h, nkrev);
				assert(k != KMHASH_END(h));
				if (nknum == krev) {
					assert(nkrev == knum);
				}
				rc = knum >> lmc;
				assert(rc >= 0 && rc < 4);
				if (nknum < nkrev)
					assert((h->adjs[k] >> ((rc ^ 3) + 4)) & 1);
				else
					assert((h->adjs[k] >> (rc ^ 3)) & 1);
				++sum_deg;
			}

			if ((h->adjs[i] >> (c + 4)) & 1) {
				nknum = ((krev << 2) & kmask) | c;
				nkrev = (knum >> 2) | ((k31key_t)(c ^ 3) << lmc);
				if (nknum < nkrev)
					k = k31hash_get(h, nknum);
				else
					k = k31hash_get(h, nkrev);
				assert(k != KMHASH_END(h));
				if (nknum == knum) {
					assert(nkrev == krev);
				}
				rc = krev >> lmc;
				assert(rc >= 0 && rc < 4);
				if (nknum < nkrev)
					assert((h->adjs[k] >> ((rc ^ 3) + 4)) & 1);
				else
					assert((h->adjs[k] >> (rc ^ 3)) & 1);
				++sum_deg;
			}
		}
	}

	__VERBOSE("Check edge kmhash done. Sum degree = %lu\n", sum_deg);
}

void build_k31_table_lazy(struct opt_count_t *opt, struct k31hash_t *h, int ksize)
{
	__VERBOSE("\nCounting %d-mer\n", ksize);
	count_kmer_lazy(opt, h, ksize);
	__VERBOSE("\n");
	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer: %llu\n", ksize,
					(long long unsigned)h->n_item);
	k31_check_edge(h, ksize);
	k31_test_hash(h);

	/* Filter singleton kmer */
	__VERBOSE("Filtering singleton %d-mer\n", ksize);
	k31hash_filter(h, 1);
	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %llu\n",
					ksize, (long long unsigned)h->n_item);

	/* Correct edges */
	k31_correct_edge(h, ksize);

	k31_test_hash(h);
	k31_check_sum(h);
	k31_check_edge(h, ksize);
}

void k31_test_process(struct opt_count_t *opt)
{
	struct k31hash_t *hm;
	hm = calloc(1, sizeof(struct k31hash_t));
	__VERBOSE("Lazy count %d-mer with max word size %lu\n",
		opt->kmer_master, (long unsigned)sizeof(k31key_t));
	build_k31_table_lazy(opt, hm, opt->kmer_master);
	test_kmer_count(opt, opt->kmer_master);
}

