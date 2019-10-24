#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "atomic.h"
#include "barcode_hash.h"
#include "fastq_producer.h"
#include "io_utils.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"
#include "include/bwa.h"
#include "include/bwamem.h"
#include "basic_resolve.h"
#include "scaffolding/global_params.h"
#include "barcode_builder.h"

struct bccount_bundle_t {
    struct asm_graph_t *g;
    struct dqueue_t *q;
    uint32_t aux_build;
    bwaidx_t *bwa_idx;
    mem_opt_t *bwa_opt;
    pthread_mutex_t *lock;
    int mapper_algo;
};

struct pathcount_bundle_t {
    struct dqueue_t *q;
    bwaidx_t *bwa_idx;
    mem_opt_t *bwa_opt;
    khash_t(contig_count) *count_cand;
    khash_t(contig_count) *count_err;
};

mem_opt_t *asm_memopt_init() {
	mem_opt_t *o;
	o = calloc(1, sizeof(mem_opt_t));
	o->flag = 0;
	o->a = 1;
	o->b = 2;
	o->o_del = o->o_ins = 3;
	o->e_del = o->e_ins = 1;
	o->w = 100;
	o->T = 30;
	o->zdrop = 100;
	o->pen_unpaired = 17;
	o->pen_clip5 = o->pen_clip3 = 5;

	o->max_mem_intv = 20;

	o->min_seed_len = 19;
	o->split_width = 10;
	o->max_occ = 500;
	o->max_chain_gap = 10000;
	o->max_ins = 10000;
	o->mask_level = 0.50;
	o->drop_ratio = 0.50;
	o->XA_drop_ratio = 0.80;
	o->split_factor = 1.5;
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->max_XA_hits = 5;
	o->max_XA_hits_alt = 200;
	o->max_matesw = 50;
	o->mask_level_redun = 0.95;
	o->min_chain_weight = 0;
	o->max_chain_extend = 1 << 30;
	o->mapQ_coef_len = 50;
	o->mapQ_coef_fac = log(o->mapQ_coef_len);
	bwa_fill_scmat(o->a, o->b, o->mat);
	return o;
}

void *barcode_buffer_iterator(void *data);

void *pathcount_buffer_iterator(void *data);

uint64_t get_barcode(struct read_t *r) {
	char *s, *p;
	s = r->info;
	if (s == NULL)
		return (uint64_t) -1;
	int i;
	uint64_t ret = 0;
	p = strstr(s, "BX:Z:");
	if (p == NULL)
		return (uint64_t) -1;
	p += 5;
	for (i = 0; p[i] && !__is_sep(p[i]); ++i)
		ret = ret * 5 + nt4_table[(int) p[i]];
	return ret;
}

static uint32_t count_shared_bc(struct barcode_hash_t *t1, struct barcode_hash_t *t2) {
	uint32_t i, k, ret = 0;
	for (i = 0; i < t1->size; ++i) {
		if (t1->keys[i] == (uint64_t) -1)
			continue;
		k = barcode_hash_get(t2, t1->keys[i]);
		ret += (int) (k != BARCODE_HASH_END(t2));
	}
	return ret;
}

uint32_t count_shared_bc_unique(struct barcode_hash_t *t1, struct barcode_hash_t *t2) {
	uint32_t i, k, ret = 0;
	for (i = 0; i < t1->size; ++i) {
		if (t1->keys[i] == (uint64_t) -1 || t1->cnts[i] == 0)
			continue;
		k = barcode_hash_get(t2, t1->keys[i]);
		if (k != BARCODE_HASH_END(t2))
			ret += (int) (t2->cnts[k] == 1);
	}
	return ret;
}

void print_test_barcode_superior(struct asm_graph_t *g, gint_t e1,
                                 gint_t e2, gint_t e2a) {
	printf("--------------- TEST %ld <-> (%ld, %ld) -----------------\n",
	       e1, e2, e2a);
	struct barcode_hash_t *h1, *h2, *h2a;
	h1 = g->edges[e1].barcodes;
	h2 = g->edges[e2].barcodes;
	h2a = g->edges[e2a].barcodes;
	printf("Number of barcode of %ld: %u\n", e1, h1->n_item);
	printf("Number of barcode of %ld: %u\n", e2, h2->n_item);
	printf("Number of barcode of %ld: %u\n", e2a, h2a->n_item);

	uint32_t share_1_2, share_1_2a, share_1_2_2a, i, k2, k2a;
	share_1_2 = share_1_2a = share_1_2_2a = 0;
	for (i = 0; i < h1->size; ++i) {
		if (h1->keys[i] == (uint64_t) -1)
			continue;
		k2 = barcode_hash_get(h2, h1->keys[i]);
		k2a = barcode_hash_get(h2a, h1->keys[i]);
		share_1_2 += (k2 != BARCODE_HASH_END(h2));
		share_1_2a += (k2a != BARCODE_HASH_END(h2a));
		share_1_2_2a += (k2 != BARCODE_HASH_END(h2) &&
		                 k2a != BARCODE_HASH_END(h2a));
	}

	printf("**** Number of share barcode of %ld and %ld: %u\n", e1, e2, share_1_2);
	printf("**** Number of share barcode of %ld and %ld: %u\n", e1, e2a, share_1_2a);
	printf("**** Number of share barcode of the three: %u\n", share_1_2_2a);
	printf("-----------------------------------------------------------\n");
}

void print_test_barcode_edge(struct asm_graph_t *g, gint_t e1, gint_t e2) {
	printf("---------------- TEST %ld <-> %ld-------------------\n", e1, e2);
	struct barcode_hash_t *h1, *h2;
	h1 = g->edges[e1].barcodes + 1;
	h2 = g->edges[e2].barcodes + 1;
	printf("Number of barcode of %ld: %u\n", e1, h1->n_item);
	printf("Number of barcode of %ld: %u\n", e2, h2->n_item);

	uint32_t cnt = count_shared_bc(h1, h2);
	printf("Number of shared barcode: %u\n", cnt);
	printf("-----------------------------------------------------------\n");
	// printf("Ratio = %.3f\n", get_barcode_ratio(g, e1, e2));

	// printf("Number of unique mapped barcode of %ld: %u\n", e1, h1->n_unique);
	// printf("Number of unique mapped barcode of %ld: %u\n", e2, h2->n_unique);
	// cnt = count_shared_bc_unique(h1, h2);
	// printf("Number of unique mapped shared barcode: %u\n", cnt);
	// printf("Ratio = %.3f\n", get_barcode_ratio_unique(g, e1, e2));
}

void init_barcode_graph(struct asm_graph_t *g, int mapper_algo) {
	g->aux_flag |= ASM_HAVE_BARCODE;
	if (mapper_algo == FOR_SCAFFOLD)
		g->aux_flag |= ASM_HAVE_BARCODE_SCAF;
	gint_t e;
	for (e = 0; e < g->n_e; ++e) {
		pthread_mutex_init(&(g->edges[e].lock), NULL);
		g->edges[e].barcodes = calloc(3, sizeof(struct barcode_hash_t));
		barcode_hash_init(g->edges[e].barcodes, 4);
		barcode_hash_init(g->edges[e].barcodes + 1, 4);
		barcode_hash_init(g->edges[e].barcodes + 2, 4);
		if (mapper_algo == FOR_SCAFFOLD) {
			barcode_hash_init(&g->edges[e].barcodes_scaf, 4);
			barcode_hash_init(&g->edges[e].barcodes_scaf2, 4);
		}
	}
}

void count_readpair_path(int n_threads, struct read_path_t *rpath,
                         const char *fasta_path, khash_t(contig_count) *count_cand) {
	bwa_idx_build(fasta_path, fasta_path, BWTALGO_AUTO, 500000000);
	bwaidx_t *bwa_idx = bwa_idx_load(fasta_path, BWA_IDX_ALL);
	mem_opt_t *bwa_opt = asm_memopt_init();
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int i;
	struct producer_bundle_t *producer_bundles;
	__VERBOSE_LOG("DEBUG", "NTHREAD = %d\n", n_threads);
	producer_bundles = init_fastq_pair(n_threads, 1,
	                                   &(rpath->R1_path), &(rpath->R2_path));

	struct pathcount_bundle_t *worker_bundles;
	worker_bundles = malloc(n_threads * sizeof(struct pathcount_bundle_t));

	for (i = 0; i < n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].bwa_idx = bwa_idx;
		worker_bundles[i].bwa_opt = bwa_opt;
		worker_bundles[i].count_cand = count_cand;
		worker_bundles[i].count_err = NULL;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(1, sizeof(pthread_t));
	worker_threads = calloc(n_threads, sizeof(pthread_t));

	pthread_create(producer_threads, &attr, fastq_producer, producer_bundles);

	for (i = 0; i < n_threads; ++i)
		pthread_create(worker_threads + i, &attr, pathcount_buffer_iterator,
		               worker_bundles + i);

	pthread_join(producer_threads[0], NULL);

	for (i = 0; i < n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	free_fastq_pair(producer_bundles, 1);
	free(worker_bundles);

	free(producer_threads);
	free(worker_threads);

	bwa_idx_destroy(bwa_idx);
	free(bwa_opt);
}

void count_readpair_err_path(int n_threads, struct read_path_t *rpath,
                             const char *fasta_path, khash_t(contig_count) *count_cand,
                             khash_t(contig_count) *count_err) {
	bwa_idx_build(fasta_path, fasta_path, BWTALGO_AUTO, 500000000);
	bwaidx_t *bwa_idx = bwa_idx_load(fasta_path, BWA_IDX_ALL);
	mem_opt_t *bwa_opt = asm_memopt_init();
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int i;
	struct producer_bundle_t *producer_bundles;
	__VERBOSE_LOG("DEBUG", "NTHREAD = %d\n", n_threads);
	producer_bundles = init_fastq_pair(n_threads, 1,
	                                   &(rpath->R1_path), &(rpath->R2_path));

	struct pathcount_bundle_t *worker_bundles;
	worker_bundles = malloc(n_threads * sizeof(struct pathcount_bundle_t));

	for (i = 0; i < n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].bwa_idx = bwa_idx;
		worker_bundles[i].bwa_opt = bwa_opt;
		worker_bundles[i].count_cand = count_cand;
		worker_bundles[i].count_err = count_err;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(1, sizeof(pthread_t));
	worker_threads = calloc(n_threads, sizeof(pthread_t));

	pthread_create(producer_threads, &attr, fastq_producer, producer_bundles);

	for (i = 0; i < n_threads; ++i)
		pthread_create(worker_threads + i, &attr, pathcount_buffer_iterator,
		               worker_bundles + i);

	pthread_join(producer_threads[0], NULL);

	for (i = 0; i < n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	free_fastq_pair(producer_bundles, 1);
	free(worker_bundles);

	free(producer_threads);
	free(worker_threads);

	bwa_idx_destroy(bwa_idx);
	free(bwa_opt);
}

void construct_aux_info(struct opt_proc_t *opt, struct asm_graph_t *g,
                        struct read_path_t *rpath, const char *fasta_path, uint32_t aux_build, int mapper_algo) {
	log_info("Construct aux info with aux_build: %d", aux_build);
	if (aux_build | ASM_BUILD_BARCODE)
		init_barcode_graph(g, mapper_algo);
	bwa_idx_build(fasta_path, fasta_path, BWTALGO_AUTO, 500000000);
	bwaidx_t *bwa_idx = bwa_idx_load(fasta_path, BWA_IDX_ALL);
	mem_opt_t *bwa_opt = asm_memopt_init();
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int i;
	struct producer_bundle_t *producer_bundles;
	producer_bundles = init_fastq_pair(opt->n_threads, 1,
	                                   &(rpath->R1_path), &(rpath->R2_path));

	struct bccount_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct bccount_bundle_t));

	pthread_mutex_t lock;
	pthread_mutex_init(&lock, NULL);
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].g = g;
		worker_bundles[i].bwa_idx = bwa_idx;
		worker_bundles[i].bwa_opt = bwa_opt;
		worker_bundles[i].lock = &lock;
		worker_bundles[i].aux_build = aux_build;
		worker_bundles[i].mapper_algo = mapper_algo;
	}

	pthread_t *producer_threads, *worker_threads;
	/* FIXME: not actually opt->n_files, noob */

	opt->n_files = 1;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_producer,
		               producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, barcode_buffer_iterator,
		               worker_bundles + i);

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	free_fastq_pair(producer_bundles, 1);
	free(worker_bundles);

	free(producer_threads);
	free(worker_threads);

	bwa_idx_destroy(bwa_idx);
	free(bwa_opt);
}

static inline int infer_bw(int l1, int l2, int score, int a, int q, int r) {
	int w;
	if (l1 == l2 && l1 * a - score < (q + r - a) << 1)
		return 0; // to get equal alignment length, we need at least two gaps
	w = ((double) ((l1 < l2 ? l1 : l2) * a - score - q) / r + 2.);
	if (w < abs(l1 - l2)) w = abs(l1 - l2);
	return w;
}

static int asm_get_score(const int8_t mat[25], int o_del, int e_del, int o_ins, int e_ins,
                         int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb,
                         int64_t re) {
	uint8_t *rseq;
	int64_t rlen;
	int i, score;
	if (l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac)) return 0;
	rseq = bns_get_seq(l_pac, pac, rb, re, &rlen);
	score = 0;
	if (re - rb != rlen) goto ret_get_score;
	if (l_query == re - rb && w_ == 0) { /* linear check, haha */
		for (i = 0, score = 0; i < l_query; ++i)
			score += mat[rseq[i] * 5 + query[i]];
	} else {
		int w, max_gap, max_ins, max_del, min_w;
		max_ins = (int) ((double) (((l_query + 1) >> 1) * mat[0] - o_ins) / e_ins + 1.);
		max_del = (int) ((double) (((l_query + 1) >> 1) * mat[0] - o_del) / e_del + 1.);
		max_gap = __max(max_ins, max_del);
		max_gap = __max(max_gap, 1);
		w = (max_gap + abs((int) rlen - l_query) + 1) >> 1;
		w = __min(w, w_);
		min_w = abs((int) rlen - l_query) + 3;
		w = __min(min_w, w);
		score = ksw_global2(l_query, query, rlen, rseq, 5, mat, o_del, e_del, o_ins, e_ins, w, NULL, NULL);
	}
	ret_get_score:
	free(rseq);
	return score;
}

struct asm_align_t asm_reg2aln(const mem_opt_t *opt, const bntseq_t *bns,
                               const uint8_t *pac, int l_query, const uint8_t *query, const mem_alnreg_t *ar) {
	struct asm_align_t a;
	memset(&a, 0, sizeof(struct asm_align_t));
	if (ar == 0 || ar->rb < 0 || ar->re < 0) { /* unmapped record */
		a.rid = -1;
		return a;
	}
	// fprintf(stderr, "%ld %ld\n", ar->rb, ar->re);
	// fprintf(stderr, "%d %d\n", ar->qb, ar->qe);
	// fprintf(stderr, "%d %d\n", ar->score, ar->truesc);
	int64_t rb, re, pos;
	int is_rev, tmp, w2, score, last_sc, i, qb, qe, l_ref,
		clip5, clip3, dist5, dist3, ext5, ext3;
	qb = ar->qb;
	qe = ar->qe;
	rb = ar->rb;
	re = ar->re;
	l_ref = re - rb;
	pos = bns_depos(bns, rb < bns->l_pac ? rb : re - 1, &is_rev);
	// fprintf(stderr, "%ld\n", pos);
	a.rid = bns_pos2rid(bns, pos);
	// fprintf(stderr, "%d\n", a.rid);
	clip5 = is_rev ? l_query - qe : qb;
	clip3 = is_rev ? qb : l_query - qe;
	dist5 = pos - bns->anns[a.rid].offset;
	dist3 = bns->anns[a.rid].offset + bns->anns[a.rid].len - pos - l_ref;
	ext5 = __min(clip5, dist5);
	ext3 = __min(clip3, dist3);
	// fprintf(stderr, "lquery = %d; ext5 = %d, ext3 = %d\n", l_query, ext5, ext3);
	if (ext5 >= 10 || ext3 >= 10) {
		a.rid = -1;
		return a;
	}
	clip5 -= ext5;
	clip3 -= ext3;
	if (is_rev) {
		qe += ext5;
		re += ext5;
		qb -= ext3;
		rb -= ext3;
	} else {
		qe += ext3;
		re += ext3;
		qb -= ext5;
		rb -= ext5;
	}
	/* query should be in nt4 encoding */
	tmp = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_del, opt->e_del);
	w2 = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_ins, opt->e_ins);
	w2 = __max(w2, tmp);
	if (w2 > opt->w) w2 = __min(w2, ar->w);
	i = 0;
	last_sc = -(1 << 30);
	do {
		w2 = __min(w2, opt->w << 2);
		// fprintf(stderr, "w2 = %d\n", w2);
		score = asm_get_score(opt->mat, opt->o_del, opt->e_del,
		                      opt->o_ins, opt->e_ins, w2, bns->l_pac, pac, qe - qb,
		                      (uint8_t *) &query[qb], rb, re);
		// fprintf(stderr, "score = %d\n", score);
		if (score == last_sc || w2 == (opt->w << 2))
			break;
		last_sc = score;
		w2 <<= 1;
	} while (++i < 3 && score < ar->truesc - opt->a);
	if (score > 0) {
		a.score = score;
		a.aligned = l_query - clip5 - clip3;
		pos = bns_depos(bns, rb < bns->l_pac ? rb : re - 1, &is_rev);
		a.pos = pos - bns->anns[a.rid].offset;
	} else {
		a.rid = -1;
	}
	a.strand = is_rev;
	return a;
}

#define FASTA_REF_SEQ        0
#define FASTA_REF_QRY        1
#define FASTA_REF_UNKNOWN        2

struct fasta_ref_t {
    int type;
    gint_t e1;
    gint_t e2;
    int pos1;
    int pos2;
};

static inline void parse_int_array(const char *s, gint_t *a) {
	int i, k;
	gint_t num;
	for (i = k = 0, num = 0; s[i]; ++i) {
		if (s[i] == '_') {
			a[k++] = num;
			num = 0;
		} else {
			if (s[i] < '0' || s[i] > '9')
				__ERROR("Parse QRY fail: %s", s);
			num = num * 10 + s[i] - '0';
		}
	}
}

static inline struct fasta_ref_t parse_fasta_ref(const char *name) {
	struct fasta_ref_t ret;
	if (!strncmp(name, "SEQ", 3)) {
		ret.type = FASTA_REF_SEQ;
		const char *p = name + 4;
		ret.e1 = atol(p);
	} else if (!strncmp(name, "QRY", 3)) {
		ret.type = FASTA_REF_QRY;
		gint_t tmp[4];
		parse_int_array(name + 4, tmp);
		ret.e1 = tmp[0];
		ret.e2 = tmp[1];
		ret.pos1 = tmp[2];
		ret.pos2 = tmp[3];
	} else {
		ret.type = FASTA_REF_UNKNOWN;
	}
	return ret;
}

static inline void add_barcode_edge(struct asm_graph_t *g, gint_t e,
                                    int lvl, uint64_t bc) {
	pthread_mutex_lock(&g->edges[e].lock);
	barcode_hash_add(g->edges[e].barcodes + lvl, bc);
	pthread_mutex_unlock(&g->edges[e].lock);
}

static inline void add_barcode_scaffold(struct asm_graph_t *g, gint_t e, uint64_t bc) {
	pthread_mutex_lock(&g->edges[e].lock);
	barcode_hash_add(&g->edges[e].barcodes_scaf, bc);
	pthread_mutex_unlock(&g->edges[e].lock);
}

static inline void add_barcode_scaffold2(struct asm_graph_t *g, gint_t e, uint64_t bc) {
	pthread_mutex_lock(&g->edges[e].lock);
	barcode_hash_add(&g->edges[e].barcodes_scaf2, bc);
	pthread_mutex_unlock(&g->edges[e].lock);
}

static inline void add_read_count_candidate(struct asm_graph_t *g, gint_t e1, gint_t e2) {
	struct pair_contig_t key = (struct pair_contig_t) {e1, e2};
	khint_t k = kh_get(pair_contig_count, g->candidates, key);
	//if (k == kh_end(g->candidates))
	//	__ERROR("candidate is not initilized");
	if (k == kh_end(g->candidates))
		return;
	atomic_add_and_fetch32(&(kh_value(g->candidates, k).n_read), 1);
}

static inline void add_readpair_count_candidate(struct asm_graph_t *g, gint_t e1, gint_t e2) {
	struct pair_contig_t key = (struct pair_contig_t) {e1, e2};
	khint_t k = kh_get(pair_contig_count, g->candidates, key);
	// if (k == kh_end(g->candidates))
	// 	__ERROR("candidate is not initilized");
	if (k == kh_end(g->candidates))
		return;
	atomic_add_and_fetch32(&(kh_value(g->candidates, k).n_pair), 1);
}

void read_mapper(struct read_t *r1, struct read_t *r2, uint64_t bc,
                 struct bccount_bundle_t *bundle) {
	struct asm_graph_t *g = bundle->g;
	bwaidx_t *idx = bundle->bwa_idx;
	mem_opt_t *opt = bundle->bwa_opt;
	mem_alnreg_v ar1, ar2;
	uint8_t *r1_seq, *r2_seq;
	int i, k, n1, n2, count, best_score_1, best_score_2;
	ar1 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r1->len, r1->seq);
	ar2 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r2->len, r2->seq);
	// fprintf(stderr, "found alignments: n1 = %lu; n2 = %lu\n", ar1.n, ar2.n);
	// 			STAGE 1 			
	r1_seq = malloc(r1->len);
	r2_seq = malloc(r2->len);
	for (i = 0; i < r1->len; ++i)
		r1_seq[i] = nst_nt4_table[(int) r1->seq[i]];
	for (i = 0; i < r2->len; ++i)
		r2_seq[i] = nst_nt4_table[(int) r2->seq[i]];
	struct asm_align_t *p1, *p2;
	p1 = alloca(ar1.n * sizeof(struct asm_align_t));
	p2 = alloca(ar2.n * sizeof(struct asm_align_t));
	n1 = n2 = 0;
	best_score_1 = best_score_2 = -(1 << 30);
	for (i = 0; i < (int) ar1.n; ++i) {
		struct asm_align_t a;
		a = asm_reg2aln(opt, idx->bns, idx->pac, r1->len, r1_seq, ar1.a + i);
		if (a.rid == -1)
			continue;
		if (a.score > best_score_1) {
			best_score_1 = a.score;
			p1[0] = a;
			n1 = 1;
		} else if (a.score == best_score_1) {
			p1[n1++] = a;
		}
	}
	for (i = 0; i < (int) ar2.n; ++i) {
		struct asm_align_t a;
		a = asm_reg2aln(opt, idx->bns, idx->pac, r2->len, r2_seq, ar2.a + i);
		if (a.rid == -1)
			continue;
		struct fasta_ref_t r = parse_fasta_ref(idx->bns->anns[a.rid].name);
		if (a.score > best_score_2) {
			best_score_2 = a.score;
			p2[0] = a;
			n2 = 1;
		} else if (a.score == best_score_2) {
			p2[n2++] = a;
		}
	}
	if (ar1.n != 2 || ar2.n != 2|| ar1.a->score < 50 || ar2.a->score < 50)
		goto free_stage_1;
	if (bundle->aux_build & ASM_BUILD_CANDIDATE) { /* read information */
		for (i = 0; i < n1; ++i) {
			struct fasta_ref_t ref;
			ref = parse_fasta_ref(idx->bns->anns[p1[i].rid].name);
			if (ref.type != FASTA_REF_QRY)
				continue;
			if (p1[i].pos < ref.pos1 - 10 && p1[i].pos + p1[i].aligned > ref.pos2 + 10)
				add_read_count_candidate(g, ref.e1, ref.e2);
		}
		/* read pair information on candidate bridges */
		for (i = 0; i < n1; ++i) {
			struct fasta_ref_t ref;
			ref = parse_fasta_ref(idx->bns->anns[p1[i].rid].name);
			if (ref.type != FASTA_REF_QRY)
				continue;
			for (k = 0; k < n2; ++k) {
				if (p2[k].rid != p1[i].rid) /* map to the same candidate */
					continue;
				int l, r;
				l = __min(p1[i].pos, p2[i].pos);
				r = __max(p1[i].pos + p1[i].aligned, p2[i].pos + p2[i].aligned);
				if (l < ref.pos1 - 10 && r > ref.pos2 + 10)
					add_readpair_count_candidate(g, ref.e1, ref.e2);
			}
		}
		for (i = 0; i < n1; ++i) {
			struct fasta_ref_t ref1, ref2;
			ref1 = parse_fasta_ref(idx->bns->anns[p1[i].rid].name);
			if (ref1.type != FASTA_REF_SEQ)
				continue;
			for (k = 0; k < n2; ++k) {
				ref2 = parse_fasta_ref(idx->bns->anns[p2[k].rid].name);
				if (ref2.type != FASTA_REF_SEQ)
					continue;
				if (p1[i].pos + p2[k].pos < MAX_READ_FRAG_LEN &&
				    ref1.e1 != ref2.e1 &&
				    ref1.e1 != g->edges[ref2.e1].rc_id) {
					add_readpair_count_candidate(g, ref1.e1, ref2.e1);
					add_readpair_count_candidate(g, ref2.e1, ref1.e1);
				}
			}
		}
	}
	if ((bundle->aux_build & ASM_BUILD_BARCODE) && bc != (uint64_t) -1) {
		for (i = 0; i < n1; ++i) {
			struct fasta_ref_t ref;
			ref = parse_fasta_ref(idx->bns->anns[p1[i].rid].name);
			if (ref.type != FASTA_REF_SEQ)
				continue;
			if (p1[i].pos <= CONTIG_LEVEL_0) {
				add_barcode_edge(g, ref.e1, 0, bc);
				add_barcode_edge(g, ref.e1, 1, bc);
				add_barcode_edge(g, ref.e1, 2, bc);
			} else if (p1[i].pos <= CONTIG_LEVEL_1) {
				add_barcode_edge(g, ref.e1, 1, bc);
				add_barcode_edge(g, ref.e1, 2, bc);
			} else if (p1[i].pos <= CONTIG_LEVEL_2) {
				add_barcode_edge(g, ref.e1, 2, bc);
			}
		}
		for (i = 0; i < n2; ++i) {
			struct fasta_ref_t ref;
			ref = parse_fasta_ref(idx->bns->anns[p2[i].rid].name);
			if (ref.type != FASTA_REF_SEQ)
				continue;
			if (p2[i].pos <= CONTIG_LEVEL_0) {
				add_barcode_edge(g, ref.e1, 0, bc);
				add_barcode_edge(g, ref.e1, 1, bc);
				add_barcode_edge(g, ref.e1, 2, bc);
			} else if (p2[i].pos <= CONTIG_LEVEL_1) {
				add_barcode_edge(g, ref.e1, 1, bc);
				add_barcode_edge(g, ref.e1, 2, bc);
			} else if (p2[i].pos <= CONTIG_LEVEL_2) {
				add_barcode_edge(g, ref.e1, 2, bc);
			}
		}
	}
	if (bundle->aux_build & ASM_BUILD_COVERAGE) {
		for (i = 0; i < n1; ++i) {
			struct fasta_ref_t ref;
			ref = parse_fasta_ref(idx->bns->anns[p1[i].rid].name);
			if (ref.type != FASTA_REF_SEQ)
				continue;
			int added_count = __max(p1[i].aligned - g->ksize, 1);
			atomic_add_and_fetch32(&g->edges[ref.e1].count, added_count);
			atomic_add_and_fetch32(&g->edges[g->edges[ref.e1].rc_id].count, added_count);
		}
		for (i = 0; i < n2; ++i) {
			struct fasta_ref_t ref;
			ref = parse_fasta_ref(idx->bns->anns[p2[i].rid].name);
			if (ref.type != FASTA_REF_SEQ)
				continue;
			int added_count = __max(p1[i].aligned - g->ksize, 1);
			atomic_add_and_fetch32(&g->edges[ref.e1].count, added_count);
			atomic_add_and_fetch32(&g->edges[g->edges[ref.e1].rc_id].count, added_count);
		}
	}
free_stage_1:
	free(ar1.a);
	free(ar2.a);
	free(r1_seq);
	free(r2_seq);
}

void read_mapper_scaffold(struct read_t *r1, struct read_t *r2, uint64_t bc,
                          struct bccount_bundle_t *bundle) {
	struct asm_graph_t *g = bundle->g;
	bwaidx_t *idx = bundle->bwa_idx;
	mem_opt_t *opt = bundle->bwa_opt;
	mem_alnreg_v ar1, ar2;
	uint8_t *r1_seq, *r2_seq;
	int i, k, n1, n2, count, best_score_1, best_score_2;
	ar1 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r1->len, r1->seq);
	ar2 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r2->len, r2->seq);
	// fprintf(stderr, "found alignments: n1 = %lu; n2 = %lu\n", ar1.n, ar2.n);
	r1_seq = malloc(r1->len);
	r2_seq = malloc(r2->len);
	for (i = 0; i < r1->len; ++i)
		r1_seq[i] = nst_nt4_table[(int) r1->seq[i]];
	for (i = 0; i < r2->len; ++i)
		r2_seq[i] = nst_nt4_table[(int) r2->seq[i]];
	struct asm_align_t *p1, *p2;
	p1 = alloca(ar1.n * sizeof(struct asm_align_t));
	p2 = alloca(ar2.n * sizeof(struct asm_align_t));
	n1 = n2 = 0;
	best_score_1 = best_score_2 = -(1 << 30);
	for (i = 0; i < (int) ar1.n; ++i) {
		struct asm_align_t a;
		a = asm_reg2aln(opt, idx->bns, idx->pac, r1->len, r1_seq, ar1.a + i);
		if (a.rid == -1)
			continue;
		if (a.score > best_score_1) {
			best_score_1 = a.score;
			p1[0] = a;
			n1 = 1;
		} else if (a.score == best_score_1) {
			p1[n1++] = a;
		}
	}
	for (i = 0; i < (int) ar2.n; ++i) {
		struct asm_align_t a;
		a = asm_reg2aln(opt, idx->bns, idx->pac, r2->len, r2_seq, ar2.a + i);
		if (a.rid == -1)
			continue;
		struct fasta_ref_t r = parse_fasta_ref(idx->bns->anns[a.rid].name);
		if (a.score > best_score_2) {
			best_score_2 = a.score;
			p2[0] = a;
			n2 = 1;
		} else if (a.score == best_score_2) {
			p2[n2++] = a;
		}
	}

	//Something Thang need but I don't know what is it.
	for (i = 0; i < n1; ++i) {
		struct fasta_ref_t ref;
		ref = parse_fasta_ref(idx->bns->anns[p1[i].rid].name);
		if (ref.type != FASTA_REF_SEQ)
			continue;
		if (p1[i].pos <= CONTIG_LEVEL_2) {
			add_barcode_edge(g, ref.e1, 2, bc);
		}
	}
	for (i = 0; i < n2; ++i) {
		struct fasta_ref_t ref;
		ref = parse_fasta_ref(idx->bns->anns[p2[i].rid].name);
		if (ref.type != FASTA_REF_SEQ)
			continue;
		if (p2[i].pos <= CONTIG_LEVEL_2) {
			add_barcode_edge(g, ref.e1, 2, bc);
		}
	}

	//-----------------build barcode scaffold -----------------------
	// todo verify if n1<2 is best
	if (ar1.n <= 2) {
		for (int i = 0; i < n1; i++) {
			struct fasta_ref_t ref;
			ref = parse_fasta_ref(idx->bns->anns[p1[i].rid].name);
			if (ref.type != FASTA_REF_SEQ)
				continue;
			if ((uint32_t) p1[i].pos < MIN(MIN_CONTIG_BARCODE, g->edges[ref.e1].seq_len / 2)) {
				add_barcode_scaffold(g, ref.e1, bc);
			}
//            if (p1[i].pos < MIN_CONTIG_BARCODE2) {
//                add_barcode_scaffold2(g, ref.e1, bc);
//            }
		}
	}
	if (ar2.n <= 2) {
		for (int i = 0; i < n2; i++) {
			struct fasta_ref_t ref;
			ref = parse_fasta_ref(idx->bns->anns[p2[i].rid].name);
			if (ref.type != FASTA_REF_SEQ)
				continue;
			if ((uint32_t) p2[i].pos < MIN(MIN_CONTIG_BARCODE, g->edges[ref.e1].seq_len / 2)) {
				add_barcode_scaffold(g, ref.e1, bc);
			}
//            if (p2[i].pos < MIN_CONTIG_BARCODE2) {
//                add_barcode_scaffold2(g, ref.e1, bc);
//            }
		}
	}

	if (bundle->aux_build & ASM_BUILD_COVERAGE) {
		for (i = 0; i < n1; ++i) {
			struct fasta_ref_t ref;
			ref = parse_fasta_ref(idx->bns->anns[p1[i].rid].name);
			if (ref.type != FASTA_REF_SEQ)
				continue;
			int added_count = __max(p1[i].aligned - g->ksize, 1);
			atomic_add_and_fetch32(&g->edges[ref.e1].count, added_count);
			atomic_add_and_fetch32(&g->edges[g->edges[ref.e1].rc_id].count, added_count);
		}
		for (i = 0; i < n2; ++i) {
			struct fasta_ref_t ref;
			ref = parse_fasta_ref(idx->bns->anns[p2[i].rid].name);
			if (ref.type != FASTA_REF_SEQ)
				continue;
			int added_count = __max(p1[i].aligned - g->ksize, 1);
			atomic_add_and_fetch32(&g->edges[ref.e1].count, added_count);
			atomic_add_and_fetch32(&g->edges[g->edges[ref.e1].rc_id].count, added_count);
		}
	}
	free(ar1.a);
	free(ar2.a);
	free(r1_seq);
	free(r2_seq);
}

void *barcode_buffer_iterator(void *data) {
	struct bccount_bundle_t *bundle = (struct bccount_bundle_t *) data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format, mapper_algo;

	int64_t n_reads;
	int64_t *gcnt_reads;
	uint64_t barcode;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf) {
			break;
		}
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		buf1 = ext_buf->R1_buf;
		buf2 = ext_buf->R2_buf;
		input_format = ext_buf->input_format;
		mapper_algo = bundle->mapper_algo;

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
			      get_read_from_fq(&read1, buf1, &pos1) :
			      get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
			      get_read_from_fq(&read2, buf2, &pos2) :
			      get_read_from_fa(&read2, buf2, &pos2);


			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file when build barcode scaffold\n");

			barcode = get_barcode(&read1);
			if (barcode == (uint64_t) -1) {
				continue;
			}
			if (mapper_algo == FOR_SCAFFOLD)
				read_mapper_scaffold(&read1, &read2, barcode, bundle);
			else
				read_mapper(&read1, &read2, barcode, bundle);

			if (rc1 == READ_END)
				break;
		}
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

void path_mapper(struct read_t *r1, struct read_t *r2, struct pathcount_bundle_t *bundle) {
	khash_t(contig_count) *count_cand = bundle->count_cand;
	khash_t(contig_count) *count_err = bundle->count_err;
	bwaidx_t *idx = bundle->bwa_idx;
	mem_opt_t *opt = bundle->bwa_opt;
	mem_alnreg_v ar1, ar2;
	uint8_t *r1_seq, *r2_seq;
	int i, k, n1, n2, count, best_score_1, best_score_2;
	ar1 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r1->len, r1->seq);
	ar2 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r2->len, r2->seq);
	// fprintf(stderr, "found alignments: n1 = %lu; n2 = %lu\n", ar1.n, ar2.n);
	r1_seq = malloc(r1->len);
	r2_seq = malloc(r2->len);
	for (i = 0; i < r1->len; ++i)
		r1_seq[i] = nst_nt4_table[(int) r1->seq[i]];
	for (i = 0; i < r2->len; ++i)
		r2_seq[i] = nst_nt4_table[(int) r2->seq[i]];
	struct asm_align_t *p1, *p2;
	p1 = alloca(ar1.n * sizeof(struct asm_align_t));
	p2 = alloca(ar2.n * sizeof(struct asm_align_t));
	n1 = n2 = 0;
	best_score_1 = best_score_2 = -(1 << 30);
	for (i = 0; i < (int) ar1.n; ++i) {
		struct asm_align_t a;
		a = asm_reg2aln(opt, idx->bns, idx->pac, r1->len, r1_seq, ar1.a + i);
		if (a.rid == -1)
			continue;
		if (a.score > best_score_1) {
			best_score_1 = a.score;
			p1[0] = a;
			n1 = 1;
		} else if (a.score == best_score_1) {
			p1[n1++] = a;
		}
	}
	for (i = 0; i < (int) ar2.n; ++i) {
		struct asm_align_t a;
		a = asm_reg2aln(opt, idx->bns, idx->pac, r2->len, r2_seq, ar2.a + i);
		if (a.rid == -1)
			continue;
		if (a.score > best_score_2) {
			best_score_2 = a.score;
			p2[0] = a;
			n2 = 1;
		} else if (a.score == best_score_2) {
			p2[n2++] = a;
		}
	}
	//FILE *f = fopen("err.cnt", "a");
	for (i = 0; i < n1; ++i) {
		if (p1[i].aligned < r1->len)
			continue;
		int c1, c2;
		c1 = atoi(idx->bns->anns[p1[i].rid].name);
		int paired_map = 0;
		for (k = 0; k < n2; ++k) {
			if (p2[k].aligned < r2->len)
				continue;
			c2 = atoi(idx->bns->anns[p2[k].rid].name);
			if (c1 == c2 && __abs(p2[k].pos - p1[i].pos) < MAX_READ_FRAG_LEN
				&& p1[i].strand != p2[k].strand) {
				khiter_t it = kh_get(contig_count, count_cand, c1);
				if (it != kh_end(count_cand)) {
					atomic_add_and_fetch32(&kh_value(count_cand, it), 1);
					paired_map = 1;
				}
			}
			if (count_err != NULL && c1 == c2
				&& __abs(p2[k].pos - p1[i].pos) < MAX_READ_FRAG_LEN
				&& p1[i].strand == p2[k].strand) {
				khiter_t it = kh_get(contig_count, count_err, c1);
				if (it != kh_end(count_err))
					atomic_add_and_fetch32(&kh_value(count_err, it), 1);
			}
		}
		if (count_err != NULL && !paired_map) {
			khiter_t it = kh_get(contig_count, count_err, c1);
			if (it != kh_end(count_err))
				atomic_add_and_fetch32(&kh_value(count_err, it), 1);
		}
	}
	for (int i = 0; i < n2; ++i) {
		if (p2[i].aligned < r2->len)
			continue;
		int c1, c2;
		c2 = atoi(idx->bns->anns[p2[i].rid].name);
		int paired_map = 0;
		for (int j = 0; j < n1; ++j) {
			if (p1[j].aligned < r1->len)
				continue;
			c1 = atoi(idx->bns->anns[p1[j].rid].name);
			if (c1 == c2 && __abs(p2[i].pos - p1[j].pos) < MAX_READ_FRAG_LEN
			    && p1[j].strand != p2[i].strand) {
				paired_map = 1;
				break;
			}
		}
		if (count_err != NULL && !paired_map) {
			khiter_t it = kh_get(contig_count, count_err, c2);
			if (it != kh_end(count_err))
				atomic_add_and_fetch32(&kh_value(count_err, it), 1);
		}
	}
	free(ar1.a);
	free(ar2.a);
	free(r1_seq);
	free(r2_seq);
}

void *pathcount_buffer_iterator(void *data) {
	struct pathcount_bundle_t *bundle = (struct pathcount_bundle_t *) data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format, mapper_algo;

	int64_t n_reads;
	int64_t *gcnt_reads;
	uint64_t barcode;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf) {
			break;
		}
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		buf1 = ext_buf->R1_buf;
		buf2 = ext_buf->R2_buf;
		input_format = ext_buf->input_format;

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
			      get_read_from_fq(&read1, buf1, &pos1) :
			      get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
			      get_read_from_fq(&read2, buf2, &pos2) :
			      get_read_from_fa(&read2, buf2, &pos2);

			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file in pathcount\n");

			path_mapper(&read1, &read2, bundle);

			if (rc1 == READ_END)
				break;
		}
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

struct read_path_t parse_read_path_from_opt(struct opt_proc_t *opt)
{
	struct read_path_t res;
	res.R1_path = opt->files_1[0];
	res.R2_path = opt->files_2[0];
	res.idx_path = opt->files_I[0];
	return res;
}
