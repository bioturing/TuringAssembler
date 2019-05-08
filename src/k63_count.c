#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "attribute.h"
#include "dqueue.h"
#include "fastq_producer.h"
#include "get_buffer.h"
#include "k63hash.h"
#include "k63_count.h"
#include "kseq.h"
#include "test.h"
#include "utils.h"
#include "verbose.h"

struct maincount_bundle_t {
	struct dqueue_t *q;
	struct k63hash_t *h;
	struct k63hash_t *dict; // read-only small kmer dictionary
	int ksmall;
	int klarge;
	int64_t *n_reads;
	pthread_mutex_t *lock_hash;
	int is10X;
};

static void k63_dump_kmer(k63key_t key, char *seq, int l)
{
	seq[l] = '\0';
	while (l) {
		seq[--l] = nt4_char[key.bin[0] & 3];
		__k63_rshift2(key);
	}
}


/*
 * Move the kmer window along reads and add kmer to hash table
 */
static void count_lazy_from_read(struct read_t *r, struct k63hash_t *h,
					int ksize, pthread_mutex_t *lock_hash)
{
}

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
	// producer_bundles = init_fastq_PE(opt);
	producer_bundles = init_fastq_PE(opt->n_threads, opt->n_files,
						opt->files_1, opt->files_2);

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
		worker_bundles[i].is10X = opt->is10X;
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

static inline void k63_test_revc(k63key_t knum, k63key_t krev, int l)
{
	char *snum, *srev;
	snum = alloca(l + 1);
	srev = alloca(l + 1);
	int i;
	for (i = 0; i < l; ++i) {
		snum[i] = nt4_char[knum.bin[0] & 3];
		__k63_rshift2(knum);
		srev[i] = nt4_char[krev.bin[0] & 3];
		__k63_rshift2(krev);
	}
	snum[l] = srev[l] = '\0';
	for (i = 0; i < l; ++i) {
		if (nt4_table[(int)snum[i]] != (nt4_table[(int)srev[l - i -1]] ^ 3)) {
			fprintf(stderr, "snum = %s\nsrev = %s\n", snum, srev);
			assert(0);
		}
	}
}

static void k63_check_edge(struct k63hash_t *h, int ksize)
{
	kmint_t i, k;
	k63key_t knum, krev, nknum, nkrev, kmask;
	int c, lmc, rc;

	kmask.bin[0] =(uint64_t)-1;
	kmask.bin[1] = (1ull << ((ksize << 1) - 64)) - 1;
	lmc = (ksize - 1) << 1;
	kmint_t sum_deg = 0;
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
				assert(k != KMHASH_END(h));
				if (__k63_equal(nknum, krev))
					assert(__k63_equal(nkrev, knum));
				rc = knum.bin[1] >> (lmc - 64);
				assert(rc >= 0 && rc < 4);
				if (__k63_lt(nknum, nkrev))
					assert((h->adjs[k] >> ((rc ^ 3) + 4)) & 1);
				else
					assert((h->adjs[k] >> (rc ^ 3)) & 1);
				++sum_deg;
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
				assert(k != KMHASH_END(h));
				if (__k63_equal(nknum, knum))
					assert(__k63_equal(nkrev, krev));
				rc = krev.bin[1] >> (lmc - 64);
				assert(rc >= 0 && rc < 4);
				if (__k63_lt(nknum, nkrev))
					assert((h->adjs[k] >> ((rc ^ 3) + 4)) & 1);
				else
					assert((h->adjs[k] >> (rc ^ 3)) & 1);
				++sum_deg;
			}
		}
	}

	__VERBOSE("Check edge kmhash done. Sum degree = %ld\n", sum_deg);
}

void k63_correct_edge(struct k63hash_t *h, int ksize)
{
	k63key_t knum, krev, nknum, nkrev, kmask;
	kmint_t i, k;
	int c, lmc;
	kmask.bin[0] = (uint64_t)-1;
	kmask.bin[1] = (1ull << ((ksize << 1) - 64)) - 1;
	knum = krev = (k63key_t){{0ull, 0ull}};
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

void build_k63_table_lazy(struct opt_count_t *opt, struct k63hash_t *h, int ksize)
{
	__VERBOSE("\nCounting %d-mer\n", ksize);
	count_kmer_lazy(opt, h, ksize);
	__VERBOSE("\n");
	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer: %llu\n", ksize,
					(long long unsigned)h->n_item);
	// k63_check_edge(h, ksize);

	/* Filter singleton kmer */
	k63hash_filter(h, 1);
	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %llu\n",
					ksize, (long long unsigned)h->n_item);

	/* Correct edges */
	k63_correct_edge(h, ksize);
	// k63_check_edge(h, ksize);
	// recount_edge(h);
	// check_edge(h, opt->k1);
}

void k63_test_process(struct opt_count_t *opt)
{
	struct k63hash_t *hm;
	hm = calloc(1, sizeof(struct k63hash_t));
	__VERBOSE("Lazy count %d-mer with max word size %lu\n",
		opt->k1, (long unsigned)sizeof(k63key_t));
	build_k63_table_lazy(opt, hm, opt->k1);
	test_kmer_count(opt, opt->k1);
}

