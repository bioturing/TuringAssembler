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
#include "../include/bwa.h"
#include "../include/bwamem.h"

struct bccount_bundle_t {
	struct asm_graph_t *g;
	struct dqueue_t *q;
	uint32_t aux_build;
	bwaidx_t *bwa_idx;
	mem_opt_t *bwa_opt;
	pthread_mutex_t *lock;
};

mem_opt_t *asm_memopt_init()
{
	mem_opt_t *o;
	o = calloc(1, sizeof(mem_opt_t));
	o->flag = 0;
	o->a = 1; o->b = 2;
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
	o->max_chain_extend = 1<<30;
	o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
	bwa_fill_scmat(o->a, o->b, o->mat);
	return o;
}

void *barcode_buffer_iterator(void *data);

uint64_t get_barcode(struct read_t *r)
{
	char *s, *p;
	s = r->info;
	if (s == NULL)
		return (uint64_t)-1;
	int i;
	uint64_t ret = 0;
	p = strstr(s, "BX:Z:");
	if (p == NULL)
		return (uint64_t)-1;
	p += 5;
	for (i = 0; p[i] && !__is_sep(p[i]); ++i)
		ret = ret * 5 + nt4_table[(int)p[i]];
	return ret;
}

static uint32_t count_shared_bc(struct barcode_hash_t *t1, struct barcode_hash_t *t2)
{
	uint32_t i, k, ret = 0;
	for (i = 0; i < t1->size; ++i) {
		if (t1->keys[i] == (uint64_t)-1)
			continue;
		k = barcode_hash_get(t2, t1->keys[i]);
		ret += (int)(k != BARCODE_HASH_END(t2));
	}
	return ret;
}

uint32_t count_shared_bc_unique(struct barcode_hash_t *t1, struct barcode_hash_t *t2)
{
	uint32_t i, k, ret = 0;
	for (i = 0; i < t1->size; ++i) {
		if (t1->keys[i] == (uint64_t)-1 || t1->cnts[i] == 0)
			continue;
		k = barcode_hash_get(t2, t1->keys[i]);
		if (k != BARCODE_HASH_END(t2))
			ret += (int)(t2->cnts[k] == 1);
	}
	return ret;
}

void print_test_barcode_superior(struct asm_graph_t *g, gint_t e1,
						gint_t e2, gint_t e2a)
{
	printf("--------------- TEST %ld <-> (%ld, %ld) -----------------\n",
								e1, e2, e2a);
	struct barcode_hash_t *h1, *h2, *h2a;
	h1 = &g->edges[e1].barcodes;
	h2 = &g->edges[e2].barcodes;
	h2a = &g->edges[e2a].barcodes;
	printf("Number of barcode of %ld: %u\n", e1, h1->n_item);
	printf("Number of barcode of %ld: %u\n", e2, h2->n_item);
	printf("Number of barcode of %ld: %u\n", e2a, h2a->n_item);

	uint32_t share_1_2, share_1_2a, share_1_2_2a, i, k2, k2a;
	share_1_2 = share_1_2a = share_1_2_2a = 0;
	for (i = 0; i < h1->size; ++i) {
		if (h1->keys[i] == (uint64_t)-1)
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

void print_test_barcode_edge(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	printf("---------------- TEST %ld <-> %ld-------------------\n", e1, e2);
	struct barcode_hash_t *h1, *h2;
	h1 = &g->edges[e1].barcodes;
	h2 = &g->edges[e2].barcodes;
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

void print_test_pair_end(struct asm_graph_t *g, gint_t e)
{
	printf("---------------- TEST %ld-------------------\n", e);
	struct barcode_hash_t *h;
	printf("number of mate = %u\n", g->edges[e].n_mate_contigs);
	gint_t k;
	for (k = 0; k < g->edges[e].n_mate_contigs; ++k) {
		printf("mate = %ld; number of reads = %ld\n",
			g->edges[e].mate_contigs[k], g->edges[e].mate_counts[k]);
	}
	printf("-----------------------------------------------------------\n");
}

static inline uint32_t dump_edge_seq(char **seq, uint32_t *m_seq,
					struct asm_edge_t *e, uint32_t max_len)
{
	uint32_t i, j, k, len;
	// if (max_len != 0)
	// 	len = __min(max_len, e->seq_len);
	// else
	// 	len = e->seq_len;
	len = e->seq_len;
	if (*m_seq < len + 1) {
		*m_seq = len + 1;
		*seq = realloc(*seq, *m_seq);
	}
	j = k = 0;
	for (i = 0; i < len; ++i)
		(*seq)[i] = nt4_char[__binseq_get(e->seq, i)];
	(*seq)[len] = '\0';
	return len;
}

void init_contig_map_info(struct asm_graph_t *g, const char *path,
				uint32_t min_len, uint32_t max_len)
{
	FILE *fp = xfopen(path, "wb");
	gint_t e;
	char *seq = NULL;
	uint32_t seq_len = 0;
	char *buf = alloca(81);
	for (e = 0; e < g->n_e; ++e) {
		pthread_mutex_init(&(g->edges[e].lock), NULL);
		uint32_t len, k, l;
		len = dump_edge_seq(&seq, &seq_len, g->edges + e, max_len);
		if (g->edges[e].seq_len >= min_len) {
			fprintf(fp, ">%ld\n", e);
			k = 0;
			while (k < len) {
				l = __min(80, len - k);
				memcpy(buf, seq + k, l);
				buf[l] = '\0';
				fprintf(fp, "%s\n", buf);
				k += l;
			}
		}
		if (g->aux_flag & ASM_HAVE_BARCODE)
			barcode_hash_init(&g->edges[e].barcodes, 4);
		g->edges[e].n_mate_contigs = 0;
		g->edges[e].mate_contigs = NULL;
		g->edges[e].mate_counts = NULL;
		// g->edges[e].mate_barcodes = NULL;
	}
	fclose(fp);
	bwa_idx_build(path, path, BWTALGO_AUTO, 25000000);
	__VERBOSE("Done indexing contigs\n");
}

void construct_aux_info(struct opt_proc_t *opt, struct asm_graph_t *g,
	struct read_path_t *rpath, const char *fasta_path, uint32_t aux_build)
{
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
	}

	pthread_t *producer_threads, *worker_threads;
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

struct asm_align_t {
	int rid;
	int pos;
	int score;
	int aligned;
};

static inline int infer_bw(int l1, int l2, int score, int a, int q, int r)
{
	int w;
	if (l1 == l2 && l1 * a - score < (q + r - a)<<1) return 0; // to get equal alignment length, we need at least two gaps
	w = ((double)((l1 < l2? l1 : l2) * a - score - q) / r + 2.);
	if (w < abs(l1 - l2)) w = abs(l1 - l2);
	return w;
}

static int asm_get_score(const int8_t mat[25], int o_del, int e_del, int o_ins, int e_ins,
	int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re)
{
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
		max_ins = (int)((double)(((l_query + 1) >> 1) * mat[0] - o_ins) / e_ins + 1.);
		max_del = (int)((double)(((l_query + 1) >> 1) * mat[0] - o_del) / e_del + 1.);
		max_gap = __max(max_ins, max_del);
		max_gap = __max(max_gap, 1);
		w = (max_gap + abs((int)rlen - l_query) + 1) >> 1;
		w = __min(w, w_);
		min_w = abs((int)rlen - l_query) + 3;
		w = __min(min_w, w);
		score = ksw_global2(l_query, query, rlen, rseq, 5, mat, o_del, e_del, o_ins, e_ins, w, NULL, NULL);
	}
ret_get_score:
	free(rseq);
	return score;
}

struct asm_align_t asm_reg2aln(const mem_opt_t *opt, const bntseq_t *bns,
	const uint8_t *pac, int l_query, const uint8_t *query, const mem_alnreg_t *ar)
{
	struct asm_align_t a;
	memset(&a, 0, sizeof(mem_aln_t));
	if (ar == 0 || ar->rb < 0 || ar->re < 0) { /* unmapped record */
		a.rid = -1;
		return a;
	}
	int64_t qb, qe, pos;
	int is_rev, tmp, w2, score, last_sc, i, rb, re, l_ref,
				clip5, clip3, dist5, dist3, ext5, ext3;
	qb = ar->qb; qe = ar->qe;
	rb = ar->rb; re = ar->re;
	l_ref = re - rb;
	pos = bns_depos(bns, rb < bns->l_pac ? rb : re - 1, &is_rev);
	a.rid = bns_pos2rid(bns, pos);
	clip5 = is_rev ? l_query - qe : qb;
	clip3 = is_rev ? qb : l_query - qe;
	dist5 = pos - bns->anns[a.rid].offset;
	dist3 = bns->anns[a.rid + 1].offset - pos - l_ref;
	ext5 = __min(clip5, dist5);
	ext3 = __min(clip3, dist3);
	if (ext5 >= 10 || ext3 >= 10) {
		a.rid = -1;
		return a;
	}
	clip5 -= ext5; clip3 -= ext3;
	if (is_rev) {
		qe += ext5; re += ext5;
		qb -= ext3; rb -= ext3;
	} else {
		qe += ext3; re += ext3;
		qb -= ext5; rb -= ext5;
	}
	/* query should be in nt4 encoding */
	tmp = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_del, opt->e_del);
	w2 = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_ins, opt->e_ins);
	w2 = __max(w2, tmp);
	if (w2 > opt->w) w2 = __min(w2, ar->w);
	i = 0;
	last_sc = - (1 << 30);
	do {
		w2 = __min(w2, opt->w << 2);
		score = asm_get_score(opt->mat, opt->o_del, opt->e_del,
			opt->o_ins, opt->e_ins, w2, bns->l_pac, pac, qe - qb,
						(uint8_t*)&query[qb], rb, re);
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
	return a;
}

#define FASTA_REF_SEQ	0
#define FASTA_REF_QRY	1
#define FASTA_REF_UNKNOWN	2

struct fasta_ref_t {
	int type;
	gint_t e;
	gint_t e1;
	gint_t e2;
	int pos1;
	int pos2;
};

static inline struct fasta_ref_t parse_fasta_ref(const char *name)
{
	struct fasta_ref_t ret;
	if (!strncmp(name, "SEQ", 3)) {
		ret.type = FASTA_REF_SEQ;
		const char *p = name + 4;
		ret.e = atol(p);
	} else if (!strncmp(name, "QRY", 3)) {
		ret.type = FASTA_REF_QRY;
	} else {
		ret.type = FASTA_REF_UNKNOWN;
	}
	return ret;
}

static inline void add_barcode_edge(struct asm_graph_t *g, gint_t e,
							int lvl, uint64_t bc)
{
	pthread_mutex_lock(&g->edges[e].lock);
	barcode_hash_add(&g->edges[e].barcodes + lvl, bc);
	pthread_mutex_unlock(&g->edges[e].lock);
}

void read_mapper(struct read_t *r1, struct read_t *r2, uint64_t bc,
						struct bccount_bundle_t *bundle)
{
	struct asm_graph_t *g = bundle->g;
	bwaidx_t *idx = bundle->bwa_idx;
	mem_opt_t *opt = bundle->bwa_opt;
	mem_alnreg_v ar1, ar2;
	uint8_t *r1_seq, *r2_seq;
	int i, n1, n2, count, best_score_1, best_score_2;
	ar1 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r1->len, r1->seq);
	ar2 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r2->len, r2->seq);
	r1_seq = malloc(r1->len);
	r2_seq = malloc(r2->len);
	for (i = 0; i < r1->len; ++i)
		r1_seq[i] = nst_nt4_table[(int)r1->seq[i]];
	for (i = 0; i < r2->len; ++i)
		r2_seq[i] = nst_nt4_table[(int)r2->seq[i]];
	struct asm_align_t *p1, *p2;
	p1 = alloca(ar1.n * sizeof(struct asm_align_t));
	p2 = alloca(ar2.n * sizeof(struct asm_align_t));
	n1 = n2 = 0;
	best_score_1 = best_score_2 = - (1 << 30);
	for (i = 0; i < (int)ar1.n; ++i) {
		struct asm_align_t a;
		a = asm_reg2aln(opt, idx->bns, idx->pac, r1->len, r1_seq, ar1.a + i);
		if (a.rid == -1)
			continue;
		struct fasta_ref_t r = parse_fasta_ref(idx->bns->anns[a.rid].name);
		if (r.type == FASTA_REF_SEQ) {
			if (bundle->aux_build & ASM_BUILD_COVERAGE) {
				count = __max(a.aligned, 1);
				atomic_add_and_fetch64(&g->edges[r.e].count, count);
				atomic_add_and_fetch64(&g->edges[g->edges[r.e].rc_id].count, count);
			}
			if (a.score > best_score_1) {
				best_score_1 = a.score;
				p1[0] = a;
				n1 = 1;
			} else if (a.score == best_score_1) {
				p1[n1++] = a;
			}
		}
	}
	for (i = 0; i < (int)ar2.n; ++i) {
		struct asm_align_t a;
		a = asm_reg2aln(opt, idx->bns, idx->pac, r2->len, r2_seq, ar2.a + i);
		if (a.rid == -1)
			continue;
		struct fasta_ref_t r = parse_fasta_ref(idx->bns->anns[a.rid].name);
		if (r.type == FASTA_REF_SEQ) {
			if (bundle->aux_build & ASM_BUILD_COVERAGE) {
				count = __max(a.aligned, 1);
				atomic_add_and_fetch64(&g->edges[r.e].count, count);
				atomic_add_and_fetch64(&g->edges[g->edges[r.e].rc_id].count, count);
			}
			if (a.score > best_score_2) {
				best_score_2 = a.score;
				p2[0] = a;
				n2 = 1;
			} else if (a.score == best_score_2) {
				p2[n2++] = a;
			}
		}
	}
	if ((bundle->aux_build & ASM_BUILD_BARCODE) && bc != (uint64_t)-1) {
		for (i = 0; i < n1; ++i) {
			if (p1[i].pos <= CONTIG_LEVEL_1) {
				add_barcode_edge(g, p1[i].rid, 0, bc);
				add_barcode_edge(g, p1[i].rid, 1, bc);
			} else if (p1[i].pos <= CONTIG_LEVEL_2) {
				add_barcode_edge(g, p1[i].rid, 1, bc);
			}
		}
		for (i = 0; i < n2; ++i) {
			if (p2[i].pos <= CONTIG_LEVEL_1) {
				add_barcode_edge(g, p2[i].rid, 0, bc);
				add_barcode_edge(g, p2[i].rid, 1, bc);
			} else if (p2[i].pos <= CONTIG_LEVEL_2) {
				add_barcode_edge(g, p2[i].rid, 1, bc);
			}
		}
	}
	free(ar1.a);
	free(ar2.a);
	free(r1_seq);
	free(r2_seq);
}

void *barcode_buffer_iterator(void *data)
{
	struct bccount_bundle_t *bundle = (struct bccount_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	int64_t n_reads;
	int64_t *gcnt_reads;
	uint64_t barcode;

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

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, buf1, &pos1) :
				get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, buf2, &pos2) :
				get_read_from_fa(&read2, buf2, &pos2);


			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file\n");

			barcode = get_barcode(&read1);
			read_mapper(&read1, &read2, barcode, bundle);

			if (rc1 == READ_END)
				break;
		}
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

