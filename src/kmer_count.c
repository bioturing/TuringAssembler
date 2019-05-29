#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "attribute.h"
#include "dqueue.h"
#include "fastq_producer.h"
#include "get_buffer.h"
#include "k31hash.h"
#include "k63hash.h"
#include "kseq.h"
#include "test.h"
#include "utils.h"
#include "verbose.h"

struct kmer_count_bundle_t {
	struct dqueue_t *q;
	void *dst;
	void *src;
	int ksize_dst;
	int ksize_src;
	int64_t *n_reads;
	pthread_mutex_t *lock_hash;
	void (*read_process_func)(struct read_t *, struct kmer_count_bundle_t *);
	uint64_t (*barcode_calculator)(struct read_t *, struct read_t *);
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

uint64_t ust_get_barcode(struct read_t *r1, struct read_t *r2)
{
	char *s = r1->info;
	if (s == NULL)
		return (uint64_t)(-1);
	int i, k, len = 0;
	uint64_t ret = 0;
	for (i = 0; s[i]; ++i) {
		if (strncmp(s + i, "BX:Z:", 5) == 0) {
			for (k = i + 5; s[k] && !__is_sep(s[k]); ++k) {
				ret = ret * 5 + nt4_table[(int)s[k]];
				++len;
			}
			break;
		}
	}
	return ret;
}

uint64_t x10_get_barcode(struct read_t *r1, struct read_t *r2)
{
	char *s = r1->seq;
	assert(r1->len >= 23); /* 16 bp barcode and 7 bp UMI */
	uint64_t ret = 0;
	int i;
	for (i = 0; i < 16; ++i)
		ret = ret * 5 + nt4_table[(int)s[i]];
	r1->seq += 23;
	return ret;
}

uint64_t (*barcode_calculators[])(struct read_t *, struct read_t *) = {ust_get_barcode, x10_get_barcode};

/*
 * Move the kmer window along reads and add kmer to hash table
 */

void k63_add_edge_str(void *p, const char *seq, int ksize)
{
	struct k63hash_t *h = (struct k63hash_t *)p;
	int lmc, i, ci, ck;
	k63key_t kmask, knum1, knum2, krev1, krev2;
	kmask.bin[0] = (uint64_t)-1;
	kmask.bin[1] = (1ull << ((ksize << 1) - 64)) - 1;
	knum1 = krev1 = (k63key_t){{0ull, 0ull}};
	lmc = (ksize - 1) << 1;
	for (i = 0; i < ksize; ++i) {
		ci = nt4_table[(int)seq[i]];
		__k63_lshift2(knum1);
		knum1.bin[0] |= ci;
		__k63_rshift2(krev1);
		krev1.bin[1] |= (uint64_t)(ci ^ 3) << (lmc - 64);
	}
	ci = nt4_table[(int)seq[ksize]];
	ck = nt4_table[(int)seq[0]] ^ 3;
	knum2 = knum1;
	__k63_lshift2(knum2); __k63_and(knum2, kmask);
	knum2.bin[0] |= ci;
	krev2 = krev1;
	__k63_rshift2(krev2);
	krev2.bin[1] |= (uint64_t)(ci ^ 3) << (lmc - 64);
	if (__k63_lt(knum1, krev1))
		k63hash_add_edge_simple(h, knum1, ci);
	else
		k63hash_add_edge_simple(h, krev1, ci + 4);
	if (__k63_lt(knum2, krev2))
		k63hash_add_edge_simple(h, knum2, ck + 4);
	else
		k63hash_add_edge_simple(h, krev2, ck);
}

void k31_add_edge_str(void *p, const char *seq, int ksize)
{
	struct k31hash_t *h = (struct k31hash_t *)p;
	int lmc, i, ci, ck;
	k31key_t kmask, knum1, knum2, krev1, krev2;
	lmc = (ksize - 1) << 1;
	kmask = ((k31key_t)1 << (ksize << 1)) - 1;
	knum1 = krev1 = 0;
	for (i = 0; i < ksize; ++i) {
		ci = nt4_table[(int)seq[i]];
		knum1 <<= 2;
		knum1 |= ci;
		krev1 >>= 2;
		krev1 |= ((k31key_t)(ci ^ 3) << lmc);
	}
	ci = nt4_table[(int)seq[ksize]];
	ck = nt4_table[(int)seq[0]] ^ 3;
	knum2 = ((knum1 << 2) & kmask) | ci;
	krev2 = (krev1 >> 2) | ((k31key_t)(ci ^ 3) << lmc);
	if (knum1 < krev1)
		k31hash_add_edge_simple(h, knum1, ci);
	else
		k31hash_add_edge_simple(h, krev1, ci + 4);
	if (knum2 < krev2)
		k31hash_add_edge_simple(h, knum2, ck + 4);
	else
		k31hash_add_edge_simple(h, krev2, ck);
}

void k31_count_from_k31(struct read_t *r, struct kmer_count_bundle_t *bundle)
{
	struct k31hash_t *dst, *src;
	int ksize_src, ksize_dst, ksize_edge;
	pthread_mutex_t *lock_hash = bundle->lock_hash;
	dst = (struct k31hash_t *)bundle->dst;
	src = (struct k31hash_t *)bundle->src;
	ksize_src = bundle->ksize_src;
	ksize_dst = bundle->ksize_dst;
	ksize_edge = ksize_dst + 1;

	kmint_t k;
	int i, last, cnt_small, ci, ck, len;
	int lmc_src, lmc_dst, num_kdst_diff, num_kedge_diff;
	char *seq;
	k31key_t knum_src, knum_dst, krev_src, krev_dst, pknum_dst, pkrev_dst;
	k31key_t kmask_src, kmask_dst;
	len = r->len;
	seq = r->seq;

	kmask_src = ((k31key_t)1 << (ksize_src << 1)) - 1;
	kmask_dst = ((k31key_t)1 << (ksize_dst << 1)) - 1;
	knum_src = krev_src = knum_dst = krev_dst = 0;
	pknum_dst = pkrev_dst = 0;
	last = cnt_small = 0;
	num_kdst_diff = ksize_dst - ksize_src + 1;
	num_kedge_diff = ksize_edge - ksize_src + 1;
	lmc_src = (ksize_src - 1) << 1;
	lmc_dst = (ksize_dst - 1) << 1;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		knum_src = (knum_src << 2) & kmask_src;
		krev_src >>= 2;
		knum_dst = (knum_dst << 2) & kmask_dst;
		krev_dst >>= 2;
		if (ci < 4) {
			knum_src |= ci;
			knum_dst |= ci;
			krev_src |= ((k31key_t)(ci ^ 3) << lmc_src);
			krev_dst |= ((k31key_t)(ci ^ 3) << lmc_dst);
			++last;
			if (last >= ksize_src) {
				if (knum_src < krev_src)
					k = k31hash_get(src, knum_src);
				else
					k = k31hash_get(src, krev_src);
				if (k != KMHASH_END(src))
					++cnt_small;
				else
					cnt_small = 0;
			} else {
				cnt_small = 0;
			}
		} else {
			last = 0;
			cnt_small = 0;
		}
		if (cnt_small >= num_kdst_diff) {
			if (knum_dst < krev_dst)
				k31hash_put_adj(dst, knum_dst, lock_hash);
			else
				k31hash_put_adj(dst, krev_dst, lock_hash);
		}
		if (cnt_small >= num_kedge_diff) {
			ck = nt4_table[(int)seq[i - ksize_dst]] ^ 3;
			if (pknum_dst < pkrev_dst)
				k31hash_add_edge(dst, pknum_dst, ci, lock_hash);
			else
				k31hash_add_edge(dst, pkrev_dst, ci + 4, lock_hash);
			if (knum_dst < krev_dst)
				k31hash_add_edge(dst, knum_dst, ck + 4, lock_hash);
			else
				k31hash_add_edge(dst, krev_dst, ck, lock_hash);
		}
		pknum_dst = knum_dst;
		pkrev_dst = krev_dst;
	}
}

void k31_count_from_scratch(struct read_t *r, struct kmer_count_bundle_t *bundle)
{
	struct k31hash_t *h = (struct k31hash_t *)bundle->dst;
	int ksize = bundle->ksize_dst;
	pthread_mutex_t *lock_hash = bundle->lock_hash;
	int i, last, ci, ck, len, lmc, kedge;
	char *seq;
	len = r->len;
	seq = r->seq;

	k31key_t knum, krev, pknum, pkrev, kmask;
	kmask = ((k31key_t)1 << (ksize << 1)) - 1;
	knum = krev = pknum = pkrev = 0;
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

void k63_count_from_scratch(struct read_t *r, struct kmer_count_bundle_t *bundle)
{
	struct k63hash_t *h = (struct k63hash_t *)bundle->dst;
	int ksize = bundle->ksize_dst;
	pthread_mutex_t *lock_hash = bundle->lock_hash;
	int i, last, ci, ck, len, lmc, kedge;
	char *seq;
	len = r->len;
	seq = r->seq;

	k63key_t knum, krev, pknum, pkrev, kmask;
	kmask.bin[0] = (uint64_t)-1;
	kmask.bin[1] = (1ull << ((ksize << 1) - 64)) - 1;
	knum = krev = pknum = pkrev = (k63key_t){{0ull, 0ull}};
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

void k63_count_from_k31(struct read_t *r, struct kmer_count_bundle_t *bundle)
{
	struct k31hash_t *src;
	struct k63hash_t *dst;
	int ksize_src, ksize_dst, ksize_edge;
	pthread_mutex_t *lock_hash = bundle->lock_hash;
	dst = (struct k63hash_t *)bundle->dst;
	src = (struct k31hash_t *)bundle->src;
	ksize_src = bundle->ksize_src;
	ksize_dst = bundle->ksize_dst;
	ksize_edge = ksize_dst + 1;

	kmint_t k;
	int i, last, cnt_small, ci, ck, len;
	int lmc_src, lmc_dst, num_kdst_diff, num_kedge_diff;
	char *seq;
	k31key_t knum_src, krev_src, kmask_src;
	k63key_t knum_dst, krev_dst, kmask_dst, pknum_dst, pkrev_dst;
	len = r->len;
	seq = r->seq;

	kmask_src = ((k31key_t)1 << (ksize_src << 1)) - 1;
	knum_src = krev_src = 0;
	kmask_dst.bin[0] = (uint64_t)-1;
	kmask_dst.bin[1] = (1ull << ((ksize_dst << 1) - 64)) - 1;
	knum_dst = krev_dst = pknum_dst = pkrev_dst = (k63key_t){{0ull, 0ull}};
	last = cnt_small = 0;
	num_kdst_diff = ksize_dst - ksize_src + 1;
	num_kedge_diff = ksize_edge - ksize_src + 1;
	lmc_src = (ksize_src - 1) << 1;
	lmc_dst = (ksize_dst - 1) << 1;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		knum_src = (knum_src << 2) & kmask_src;
		krev_src >>= 2;
		__k63_lshift2(knum_dst); __k63_and(knum_dst, kmask_dst);
		__k63_rshift2(krev_dst);
		if (ci < 4) {
			knum_src |= ci;
			krev_src |= ((k31key_t)(ci ^ 3) << lmc_src);
			knum_dst.bin[0] |= ci;
			/* please make sure 128 < lmc <= 64 */
			krev_dst.bin[1] |= (uint64_t)(ci ^ 3) << (lmc_dst - 64);
			++last;
			if (last >= ksize_src) {
				if (knum_src < krev_src)
					k = k31hash_get(src, knum_src);
				else
					k = k31hash_get(src, krev_src);
				if (k != KMHASH_END(src))
					++cnt_small;
				else
					cnt_small = 0;
			} else {
				cnt_small = 0;
			}
		} else {
			last = 0;
			cnt_small = 0;
		}
		if (cnt_small >= num_kdst_diff) {
			if (__k63_lt(knum_dst, krev_dst))
				k63hash_put_adj(dst, knum_dst, lock_hash);
			else
				k63hash_put_adj(dst, krev_dst, lock_hash);
		}
		if (cnt_small >= num_kedge_diff) {
			ck = nt4_table[(int)seq[i - ksize_dst]] ^ 3;
			if (__k63_lt(pknum_dst, pkrev_dst))
				k63hash_add_edge(dst, pknum_dst, ci, lock_hash);
			else
				k63hash_add_edge(dst, pkrev_dst, ci + 4, lock_hash);
			if (__k63_lt(knum_dst, krev_dst))
				k63hash_add_edge(dst, knum_dst, ck + 4, lock_hash);
			else
				k63hash_add_edge(dst, krev_dst, ck, lock_hash);
		}
		pknum_dst = knum_dst;
		pkrev_dst = krev_dst;
	}
}

void k63_count_from_k63(struct read_t *r, struct kmer_count_bundle_t *bundle)
{
	struct k63hash_t *src;
	struct k63hash_t *dst;
	int ksize_src, ksize_dst, ksize_edge;
	pthread_mutex_t *lock_hash = bundle->lock_hash;
	dst = (struct k63hash_t *)bundle->dst;
	src = (struct k63hash_t *)bundle->src;
	ksize_src = bundle->ksize_src;
	ksize_dst = bundle->ksize_dst;
	ksize_edge = ksize_dst + 1;

	kmint_t k;
	int i, last, cnt_small, ci, ck, len;
	int lmc_src, lmc_dst, num_kdst_diff, num_kedge_diff;
	char *seq;
	k63key_t knum_src, krev_src, kmask_src;
	k63key_t knum_dst, krev_dst, kmask_dst, pknum_dst, pkrev_dst;
	len = r->len;
	seq = r->seq;

	kmask_src.bin[0] = (uint64_t)-1;
	kmask_src.bin[1] = (1ull << ((ksize_src << 1) - 64)) - 1;
	knum_src = krev_src = (k63key_t){{0ull, 0ull}};
	kmask_dst.bin[0] = (uint64_t)-1;
	kmask_dst.bin[1] = (1ull << ((ksize_dst - 1) - 64)) - 1;
	knum_dst = krev_dst = pknum_dst = pkrev_dst = (k63key_t){{0ull, 0ull}};
	last = cnt_small = 0;
	num_kdst_diff = ksize_dst - ksize_src + 1;
	num_kedge_diff = ksize_edge - ksize_src + 1;
	lmc_src = (ksize_src - 1) << 1;
	lmc_dst = (ksize_dst - 1) << 1;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		__k63_lshift2(knum_src); __k63_and(knum_src, kmask_src);
		__k63_rshift2(krev_src);
		__k63_lshift2(knum_dst); __k63_and(knum_dst, kmask_dst);
		__k63_rshift2(krev_dst);
		if (ci < 4) {
			knum_src.bin[0] |= ci;
			krev_src.bin[1] |= (uint64_t)(ci ^ 3) << (lmc_src - 64);
			knum_dst.bin[0] |= ci;
			krev_dst.bin[1] |= (uint64_t)(ci ^ 3) << (lmc_dst - 64);
			++last;
			if (last >= ksize_src) {
				if (__k63_lt(knum_src, krev_src))
					k = k63hash_get(src, knum_src);
				else
					k = k63hash_get(src, krev_src);
				if (k != KMHASH_END(src))
					++cnt_small;
				else
					cnt_small = 0;
			} else {
				cnt_small = 0;
			}
		} else {
			last = 0;
			cnt_small = 0;
		}
		if (cnt_small >= num_kdst_diff) {
			if (__k63_lt(knum_dst, krev_dst))
				k63hash_put_adj(dst, knum_dst, lock_hash);
			else
				k63hash_put_adj(dst, krev_dst, lock_hash);
		}
		if (cnt_small >= num_kedge_diff) {
			ck = nt4_table[(int)seq[i - ksize_dst]] ^ 3;
			if (__k63_lt(pknum_dst, pkrev_dst))
				k63hash_add_edge(dst, pknum_dst, ci, lock_hash);
			else
				k63hash_add_edge(dst, pkrev_dst, ci + 4, lock_hash);
			if (__k63_lt(knum_dst, krev_dst))
				k63hash_add_edge(dst, knum_dst, ck + 4, lock_hash);
			else
				k63hash_add_edge(dst, krev_dst, ck, lock_hash);
		}
		pknum_dst = knum_dst;
		pkrev_dst = krev_dst;
	}
}

static void *buffer_process(void *data)
{
	struct kmer_count_bundle_t *bundle = (struct kmer_count_bundle_t *)data;
	struct dqueue_t *q = bundle->q;

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	int64_t n_reads;
	int64_t *gcnt_reads;
	gcnt_reads = bundle->n_reads;
	void (*read_process)(struct read_t *, struct kmer_count_bundle_t *) = bundle->read_process_func;
	uint64_t (*barcode_calculator)(struct read_t *, struct read_t *) = bundle->barcode_calculator;

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
			barcode_calculator(&read1, &read2);
			read_process(&read1, bundle);
			read_process(&read2, bundle);

			if (rc1 == READ_END)
				break;
		}
		n_reads = __sync_add_and_fetch(gcnt_reads, n_reads);
		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

static void kmer_start_count(struct opt_proc_t *opt, struct kmer_count_bundle_t *dummy)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	int64_t n_reads = 0;
	int i;

	struct producer_bundle_t *producer_bundles;
	producer_bundles = init_fastq_PE(opt->n_threads, opt->n_files,
						opt->files_1, opt->files_2);

	struct kmer_count_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct kmer_count_bundle_t));

	// struct k31hash_t *dst = dummy->dst;
	// k31hash_init(dst, (kmint_t)opt->hash_size - 1, opt->n_threads, 1);
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].src = dummy->src;
		worker_bundles[i].dst = dummy->dst;
		worker_bundles[i].ksize_dst = dummy->ksize_dst;
		worker_bundles[i].ksize_src = dummy->ksize_src;
		worker_bundles[i].n_reads = &n_reads;
		worker_bundles[i].lock_hash = dummy->lock_hash + i;
		worker_bundles[i].read_process_func = dummy->read_process_func;
		worker_bundles[i].barcode_calculator = dummy->barcode_calculator;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, buffer_process,
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

static void count_k31_from_scratch(struct opt_proc_t *opt,
				struct k31hash_t *h, int ksize)
{
	struct kmer_count_bundle_t skeleton;
	k31hash_init(h, opt->hash_size - 1, opt->n_threads, 1);
	skeleton.dst = (void *)h;
	skeleton.lock_hash = h->locks;
	skeleton.ksize_dst = ksize;
	skeleton.read_process_func = k31_count_from_scratch;
	skeleton.barcode_calculator = barcode_calculators[opt->lib_type];
	kmer_start_count(opt, &skeleton);
}

static void count_k31_from_k31(struct opt_proc_t *opt, struct k31hash_t *dst,
			struct k31hash_t *src, int ksize_dst, int ksize_src)
{
	struct kmer_count_bundle_t skeleton;
	k31hash_init(dst, opt->hash_size - 1, opt->n_threads, 1);
	skeleton.dst = (void *)dst;
	skeleton.lock_hash = dst->locks;
	skeleton.src = (void *)src;
	skeleton.ksize_dst = ksize_dst;
	skeleton.ksize_src = ksize_src;
	skeleton.read_process_func = k31_count_from_k31;
	skeleton.barcode_calculator = barcode_calculators[opt->lib_type];
	kmer_start_count(opt, &skeleton);
}

static void count_k63_from_scratch(struct opt_proc_t *opt,
				struct k63hash_t *h, int ksize)
{
	struct kmer_count_bundle_t skeleton;
	k63hash_init(h, opt->hash_size - 1, opt->n_threads, 1);
	skeleton.dst = (void *)h;
	skeleton.lock_hash = h->locks;
	skeleton.ksize_dst = ksize;
	skeleton.read_process_func = k63_count_from_scratch;
	skeleton.barcode_calculator = barcode_calculators[opt->lib_type];
	kmer_start_count(opt, &skeleton);
}

static void count_k63_from_k31(struct opt_proc_t *opt, struct k63hash_t *dst,
			struct k31hash_t *src, int ksize_dst, int ksize_src)
{
	struct kmer_count_bundle_t skeleton;
	k63hash_init(dst, opt->hash_size - 1, opt->n_threads, 1);
	skeleton.dst = (void *)dst;
	skeleton.lock_hash = dst->locks;
	skeleton.src = (void *)src;
	skeleton.ksize_dst = ksize_dst;
	skeleton.ksize_src = ksize_src;
	skeleton.read_process_func = k63_count_from_k31;
	skeleton.barcode_calculator = barcode_calculators[opt->lib_type];
	kmer_start_count(opt, &skeleton);
}

static void count_k63_from_k63(struct opt_proc_t *opt, struct k63hash_t *dst,
			struct k63hash_t *src, int ksize_dst, int ksize_src)
{
	struct kmer_count_bundle_t skeleton;
	k63hash_init(dst, opt->hash_size - 1, opt->n_threads, 1);
	skeleton.dst = (void *)dst;
	skeleton.lock_hash = dst->locks;
	skeleton.src = (void *)src;
	skeleton.ksize_dst = ksize_dst;
	skeleton.ksize_src = ksize_src;
	skeleton.read_process_func = k63_count_from_k63;
	skeleton.barcode_calculator = barcode_calculators[opt->lib_type];
	kmer_start_count(opt, &skeleton);
}

struct ce_bundle_t {
	void *h;
	int ksize;
	int thread_no;
	int n_threads;
};

void *k31_correct_edge(void *data)
{
	struct ce_bundle_t *bundles = (struct ce_bundle_t *)data;
	struct k31hash_t *h = bundles->h;
	int ksize = bundles->ksize;
	int thread_no = bundles->thread_no;
	int n_threads = bundles->n_threads;
	kmint_t cap, l, r;
	cap = h->size / n_threads + 1;
	l = cap * thread_no;
	r = __min(cap * (thread_no + 1), h->size);

	k31key_t knum, krev, nknum, nkrev, kmask;
	kmint_t i, k;
	int c, lmc;
	kmask = ((k31key_t)1 << (ksize << 1)) - 1;
	lmc = (ksize - 1) << 1;

	for (i = l; i < r; ++i) {
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
				if (k == KMHASH_END(h))
					h->adjs[i] &= ~((uint8_t)1 << c);
			}
			if ((h->adjs[i] >> (c + 4)) & 1) {
				nknum = ((krev << 2) & kmask) | c;
				nkrev = (knum >> 2) | ((k31key_t)(c ^ 3) << lmc);
				if (nknum < nkrev)
					k = k31hash_get(h, nknum);
				else
					k = k31hash_get(h, nkrev);
				if (k == KMHASH_END(h))
					h->adjs[i] &= ~((uint8_t)1 << (c + 4));
			}
		}
	}
	return NULL;
}

void *k63_correct_edge(void *data)
{
	struct ce_bundle_t *bundles = (struct ce_bundle_t *)data;
	struct k63hash_t *h = bundles->h;
	int ksize = bundles->ksize;
	int thread_no = bundles->thread_no;
	int n_threads = bundles->n_threads;
	kmint_t cap, l, r;
	cap = h->size / n_threads + 1;
	l = cap * thread_no;
	r = __min(cap * (thread_no + 1), h->size);

	k63key_t knum, krev, nknum, nkrev, kmask;
	kmint_t i, k;
	int c, lmc;
	kmask.bin[0] = (uint64_t)-1;
	kmask.bin[1] = (1ull << ((ksize << 1) - 64)) - 1;
	knum = krev = (k63key_t){{0ull, 0ull}};
	lmc = (ksize - 1) << 1;

	for (i = l; i < r; ++i) {
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
				if (k == KMHASH_END(h))
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
				if (k == KMHASH_END(h))
					h->adjs[i] &= ~((uint8_t)1 << (c + 4));
			}
		}
	}
	return NULL;
}

static void kmer_correct_edge(void *h, int ksize, int n_threads,
						void *(*correct_func)(void *))
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	struct ce_bundle_t *worker_bundles;
	worker_bundles = malloc(n_threads * sizeof(struct ce_bundle_t));
	int i;
	for (i = 0; i < n_threads; ++i) {
		worker_bundles[i].h = h;
		worker_bundles[i].ksize = ksize;
		worker_bundles[i].n_threads = n_threads;
		worker_bundles[i].thread_no = i;
	}
	pthread_t *threads;
	threads = calloc(n_threads, sizeof(pthread_t));

	for (i = 0; i < n_threads; ++i)
		pthread_create(threads + i, &attr, correct_func,
							worker_bundles + i);
	for (i = 0; i < n_threads; ++i)
		pthread_join(threads[i], NULL);

	free(worker_bundles);
	free(threads);
}

void build_k31_table_from_scratch(struct opt_proc_t *opt, struct k31hash_t *h, int ksize)
{
	__VERBOSE("\nCounting %d-mer\n", ksize);
	count_k31_from_scratch(opt, h, ksize);
	__VERBOSE("\n");
	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer: %llu\n", ksize,
					(long long unsigned)h->n_item);

	/* Filter singleton kmer */
	__VERBOSE("Filtering singleton %d-mer\n", ksize);
	k31hash_filter(h, 1);
	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %llu\n",
					ksize, (long long unsigned)h->n_item);

	/* Correct edges */
	kmer_correct_edge(h, ksize, opt->n_threads, k31_correct_edge);
}

void build_k31_table_from_k31_table(struct opt_proc_t *opt,
	struct k31hash_t *dst, struct k31hash_t *src, int ksize_dst, int ksize_src)
{
	__VERBOSE("\nCounting %d-mer from pre-counted %d-mer\n", ksize_dst, ksize_src);
	count_k31_from_k31(opt, dst, src, ksize_dst, ksize_src);
	__VERBOSE("\n");
	__VERBOSE("Number of pre-filtered %d-mer: %lu\n", ksize_dst, dst->n_item);

	/* Filter singleton kmer */
	__VERBOSE("Filtering singleton\n");
	k31hash_filter(dst, 1);
	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %lu\n",
		ksize_dst, dst->n_item);

	/* Correct edges */
	kmer_correct_edge(dst, ksize_dst, opt->n_threads, k31_correct_edge);
}

// void build_k31_2step(struct opt_proc_t *opt, struct k31hash_t *h,
// 						int ksize_dst, int ksize_src)
// {
// 	if (ksize_src >= ksize_dst) {
// 		build_k31_table_from_scratch(opt, h, ksize_dst);
// 		return;
// 	}
// 	k31hash_t tmp;
// 	__VERBOSE("\n2-step counting %d-mer (from %d-mer)\n", ksize_dst, ksize_src);
// 	count_k31_from_scratch(opt, &tmp, ksize);
// 	__VERBOSE("\n");
// 	__VERBOSE_LOG("KMER COUNT" "Number of %d-mer: %llu\n", ksize_src,
// 						(long long unsigned)h->n_item);
// 	__VERBOSE("Filtering singleton %d-mer\n", ksize_src);
// 	k31hash_filter(&tmp, 1);
// 	__VERBOSE("Number of non-singleton %d-mer: %lu\n", 

// }

void build_k63_table_from_scratch(struct opt_proc_t *opt, struct k63hash_t *h, int ksize)
{
	__VERBOSE("\nCounting %d-mer\n", ksize);
	count_k63_from_scratch(opt, h, ksize);
	__VERBOSE("\n");
	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer: %lu\n", ksize,
							h->n_item);

	/* Filter singleton kmer */
	__VERBOSE("Filtering singleton %d-mer\n", ksize);
	k63hash_filter(h, 1);
	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %lu\n",
							ksize, h->n_item);

	/* Correct edges */
	kmer_correct_edge(h, ksize, opt->n_threads, k63_correct_edge);
}

void build_k63_table_from_k31_table(struct opt_proc_t *opt,
	struct k63hash_t *dst, struct k31hash_t *src, int ksize_dst, int ksize_src)
{
	__VERBOSE("\nCounting %d-mer from pre-counted %d-mer\n", ksize_dst, ksize_src);
	count_k63_from_k31(opt, dst, src, ksize_dst, ksize_src);
	__VERBOSE("\n");
	__VERBOSE("Number of pre-filtered %d-mer: %lu\n", ksize_dst, dst->n_item);

	/* Filter singleton kmer */
	__VERBOSE("Filtering singleton %d-mer\n", ksize_dst);
	k63hash_filter(dst, 1);
	k63hash_test(dst);
	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %lu\n",
							ksize_dst, dst->n_item);

	/* Correct edges */
	kmer_correct_edge(dst, ksize_dst, opt->n_threads, k63_correct_edge);
}

void build_k63_table_from_k63_table(struct opt_proc_t *opt,
	struct k63hash_t *dst, struct k63hash_t *src, int ksize_dst, int ksize_src)
{
	__VERBOSE("\nCounting %d-mer from pre-counted %d-mer\n", ksize_dst, ksize_src);
	count_k63_from_k63(opt, dst, src, ksize_dst, ksize_src);
	__VERBOSE("\n");
	__VERBOSE("Number of pre-filtered %d-mer: %lu\n", ksize_dst, dst->n_item);

	/* Filter singleton kmer */
	__VERBOSE("Filtering singleton %d-mer\n", ksize_dst);
	k63hash_filter(dst, 1);
	__VERBOSE_LOG("KMER COUNT", "Number of non-singleton %d-mer: %lu\n",
							ksize_dst, dst->n_item);

	/* Correct edges */
	kmer_correct_edge(dst, ksize_dst, opt->n_threads, k63_correct_edge);
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

// static void k31_check_sum(struct k31hash_t *h)
// {
// 	kmint_t i;
// 	uint64_t sum = 0;
// 	for (i = 0; i < h->size; ++i) {
// 		if (h->keys[i] != K31_NULL) {
// 			sum += h->keys[i];
// 		}
// 	}
// 	fprintf(stderr, "check sum = %lu\n", sum);
// }

// static void k31_check_edge(struct k31hash_t *h, int ksize)
// {
// 	kmint_t i, k;
// 	k31key_t knum, krev, nknum, nkrev, kmask;
// 	int c, lmc, rc;

// 	kmask = ((k31key_t)1 << (ksize << 1)) - 1;
// 	lmc = (ksize - 1) << 1;
// 	kmint_t sum_deg = 0;
// 	for (i = 0; i < h->size; ++i) {
// 		if (h->keys[i] == K31_NULL)
// 			continue;
// 		knum = h->keys[i];
// 		__k31_revc_num(knum, krev, ksize, kmask);
// 		for (c = 0; c < 4; ++c) {
// 			if ((h->adjs[i] >> c) & 1) {
// 				nknum = ((knum << 2) & kmask) | c;
// 				nkrev = (krev >> 2) | ((k31key_t)(c ^ 3) << lmc);
// 				if (nknum < nkrev)
// 					k = k31hash_get(h, nknum);
// 				else
// 					k = k31hash_get(h, nkrev);
// 				assert(k != KMHASH_END(h));
// 				if (nknum == krev) {
// 					assert(nkrev == knum);
// 				}
// 				rc = knum >> lmc;
// 				assert(rc >= 0 && rc < 4);
// 				if (nknum < nkrev)
// 					assert((h->adjs[k] >> ((rc ^ 3) + 4)) & 1);
// 				else
// 					assert((h->adjs[k] >> (rc ^ 3)) & 1);
// 				++sum_deg;
// 			}

// 			if ((h->adjs[i] >> (c + 4)) & 1) {
// 				nknum = ((krev << 2) & kmask) | c;
// 				nkrev = (knum >> 2) | ((k31key_t)(c ^ 3) << lmc);
// 				if (nknum < nkrev)
// 					k = k31hash_get(h, nknum);
// 				else
// 					k = k31hash_get(h, nkrev);
// 				assert(k != KMHASH_END(h));
// 				if (nknum == knum) {
// 					assert(nkrev == krev);
// 				}
// 				rc = krev >> lmc;
// 				assert(rc >= 0 && rc < 4);
// 				if (nknum < nkrev)
// 					assert((h->adjs[k] >> ((rc ^ 3) + 4)) & 1);
// 				else
// 					assert((h->adjs[k] >> (rc ^ 3)) & 1);
// 				++sum_deg;
// 			}
// 		}
// 	}

// 	__VERBOSE("Check edge kmhash done. Sum degree = %lu\n", sum_deg);
// }

