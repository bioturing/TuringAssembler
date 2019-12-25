#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "attribute.h"
#include "dqueue.h"
#include "fastq_producer.h"
#include "get_buffer.h"
#include "khash.h"
#include "kseq.h"
#include "utils.h"
#include "verbose.h"

KHASH_MAP_INIT_STR(kmap_t, uint8_t)

struct count_bundle_t {
	struct dqueue_t *q;
	khash_t(kmap_t) *h;
	int ksize;
	int64_t *n_reads;
	pthread_mutex_t *lock_hash;
};

static void count_from_read(struct read_t *r, khash_t(kmap_t) *h,
					int ksize, pthread_mutex_t *lock_hash)
{
	int len, i, c, last, j, k, ret;
	char *seq = r->seq;
	len = r->len;
	last = 0;
	for (i = 0; i < len; ++i) {
		c = nt4_table[(int)seq[i]];
		if (c < 4)
			++last;
		else
			last = 0;
		if (last >= ksize) {
			char *snum, *srev;
			snum = malloc(ksize + 1);
			srev = malloc(ksize + 1);
			j = i - ksize + 1;
			for (k = 0; k < ksize; ++k) {
				snum[k] = nt4_char[nt4_table[(int)seq[j + k]]];
				srev[k] = rev_nt4_char[nt4_table[(int)seq[i - k]]];
			}
			snum[ksize] = srev[ksize] = '\0';
			pthread_mutex_lock(lock_hash);
			if (strcmp(snum, srev) < 0) {
				khiter_t k_it = kh_put(kmap_t, h, snum, &ret);
				assert(ret != -1 && ret != 2);
				if (ret == 1)
					kh_val(h, k_it) = 0;
				else if (ret == 0) {
					kh_val(h, k_it) = 1;
					free(snum);
				}
				free(srev);
			} else {
				khiter_t k_it = kh_put(kmap_t, h, srev, &ret);
				assert(ret != -1 && ret != 2);
				if (ret == 1)
					kh_val(h, k_it) = 0;
				else if (ret == 0) {
					kh_val(h, k_it) = 1;
					free(srev);
				}
				free(snum);
			}
			pthread_mutex_unlock(lock_hash);
		}
	}
}

static void *PE_count(void *data)
{
	struct count_bundle_t *bundle = (struct count_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	khash_t(kmap_t) *h = bundle->h;

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	int64_t *gcnt_reads = bundle->n_reads;
	int ksize = bundle->ksize;
	pthread_mutex_t *lock_hash = bundle->lock_hash;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf)
			break;
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		buf1 = ext_buf->R1_buf;
		buf2 = ext_buf->R2_buf;
		input_format = ext_buf->input_format;

		int64_t n_reads = 0;
		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, buf1, &pos1) :
				get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, buf2, &pos2) :
				get_read_from_fa(&read2, buf2, &pos2);


			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				log_error("\nWrong format file pe count\n");

			++n_reads;
			count_from_read(&read1, h, ksize, lock_hash);
			count_from_read(&read2, h, ksize, lock_hash);

			if (rc1 == READ_END)
				break;
		}
		n_reads = __sync_add_and_fetch(gcnt_reads, n_reads);
		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);

}

static void count_kmer(struct opt_count_t *opt, khash_t(kmap_t) *h, int ksize)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	int64_t n_reads;
	int i;

	struct producer_bundle_t *producer_bundles;
	producer_bundles = init_fastq_pair(opt->n_threads, opt->n_files,
						opt->files_1, opt->files_2);

	struct count_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct count_bundle_t));

	n_reads = 0;
	pthread_mutex_t lock_hash;
	pthread_mutex_init(&lock_hash, NULL);

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].h = h;
		worker_bundles[i].ksize = ksize;
		worker_bundles[i].n_reads = &n_reads;
		worker_bundles[i].lock_hash = &lock_hash;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, PE_count,
				worker_bundles + i);

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	free_fastq_pair(producer_bundles, opt->n_files);
	free(worker_bundles);

	free(worker_threads);
	free(producer_threads);
}

void test_kmer_count(struct opt_count_t *opt, int ksize)
{
	khash_t(kmap_t) *h = kh_init(kmap_t);
	__VERBOSE("Naiive counting %d-mer\n", ksize);
	count_kmer(opt, h, ksize);
	__VERBOSE("\n");
	__VERBOSE_LOG("TEST", "Number of %d-mer: %u\n", ksize,
						(unsigned)kh_size(h));
	uint32_t cnt = 0;
	khint_t k;
	for (k = kh_begin(h); k != kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			if (kh_val(h, k) == 1)
				++cnt;
			free((char *)kh_key(h, k));
		}
	}
	__VERBOSE_LOG("TEST", "Number of non-singleton %d-mer: %u\n", ksize,
								(unsigned)cnt);
	kh_destroy(kmap_t, h);
}
