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
	int64_t *n_reads;
	bwaidx_t *bwa_idx;
	mem_opt_t *bwa_opt;
	uint64_t (*barcode_calculator)(struct read_t *, struct read_t *);

	pthread_mutex_t *lock;
	uint64_t *hash_sum;
};

void barcode_start_count(struct opt_proc_t *opt, struct bccount_bundle_t *ske);

uint64_t ust_get_barcode(struct read_t *r1, struct read_t *r2)
{
	char *s = r1->info;
	if (s == NULL)
		return (uint64_t)-1;
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

double get_barcode_ratio(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	struct barcode_hash_t *h1, *h2;
	h1 = &g->edges[e1].barcodes;
	h2 = &g->edges[e2].barcodes;
	uint32_t cnt = count_shared_bc(h1, h2);
	if (h1->n_item + h2->n_item - cnt == 0)
		return -1;
	return cnt * 1.0 / (h1->n_item + h2->n_item - cnt);
}

double get_barcode_ratio_unique(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	struct barcode_hash_t *h1, *h2;
	h1 = &g->edges[e1].barcodes;
	h2 = &g->edges[e2].barcodes;
	uint32_t cnt = count_shared_bc_unique(h1, h2);
	if (h1->n_unique + h2->n_unique - cnt == 0)
		return -1;
	return cnt * 1.0 / (h1->n_unique + h2->n_unique - cnt);
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

// static inline gint_t find_best_mate(struct barcode_hash_t *h)
// {
// 	gint_t best, second_best;
// 	uint32_t cnt_best, cnt_2nd_best, i;
// 	best = second_best = -1;
// 	cnt_best = cnt_2nd_best = 0;
// 	for (i = 0; i < h->size; ++i) {
// 		if (h->keys[i] == (uint64_t)-1)
// 			continue;
// 		if (h->cnts[i] > cnt_best) {
// 			cnt_2nd_best = cnt_best;
// 			second_best = best;
// 			cnt_best = h->cnts[i];
// 			best = (gint_t)h->keys[i];
// 		} else if (h->cnts[i] > cnt_2nd_best) {
// 			cnt_2nd_best = h->cnts[i];
// 			second_best = (gint_t)h->keys[i];
// 		}
// 	}
// 	if (best != -1) {
// 		if (second_best == -1 || cnt_best > cnt_2nd_best * 2)
// 			return best;
// 		return -1;
// 	} else {
// 		return -1;
// 	}
// }

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

// void graph_aux_refine(struct asm_graph_t *g)
// {
// 	gint_t e;
// 	for (e = 0; e < g->n_e; ++e) {
// 		// if (g->aux_flag & ASM_HAVE_BARCODE)
// 		// 	barcode_hash_filter(&g->edges[e].barcodes, 1);
// 		if (g->aux_flag & ASM_HAVE_READPAIR) {
// 			g->edges[e].best_mate_contigs = find_best_mate(&g->edges[e].mate_contigs);
// 			// barcode_hash_filter(&g->edges[e].mate_contigs, 0);
// 		}
// 	}
// }

void construct_aux_information(struct opt_proc_t *opt, struct asm_graph_t *g, uint32_t aux_build)
{
	extern uint64_t (*barcode_calculators[])(struct read_t *, struct read_t *);
	char fasta_prefix[1024];
	strcpy(fasta_prefix, opt->out_dir);
	strcat(fasta_prefix, "/contigs.fasta");
	g->aux_flag = 0;
	if (aux_build & ASM_BUILD_BARCODE)
		g->aux_flag |= ASM_HAVE_BARCODE;
	if (aux_build & ASM_BUILD_READPAIR)
		g->aux_flag |= ASM_HAVE_READPAIR;
	init_contig_map_info(g, fasta_prefix, 0, 0);
	// if (!(g->aux_flag & ASM_BUILD_COVERAGE))
	// 	init_contig_map_info(g, fasta_prefix, MIN_CONTIG_READPAIR, MIN_CONTIG_BARCODE);
	// else
	// 	init_contig_map_info(g, fasta_prefix, 0, 0);
	struct bccount_bundle_t ske;
	ske.g = g;
	ske.aux_build = aux_build;
	ske.barcode_calculator = barcode_calculators[opt->lib_type];
	ske.bwa_idx = bwa_idx_load(fasta_prefix, BWA_IDX_ALL);
	ske.bwa_opt = mem_opt_init();
	ske.bwa_opt->max_XA_hits = 100;
	barcode_start_count(opt, &ske);
	bwa_idx_destroy(ske.bwa_idx);
	free(ske.bwa_opt);
	// graph_aux_refine(g);
}

static inline void add_barcode_edge(struct asm_graph_t *g, gint_t e, uint64_t bc)
{
	pthread_mutex_lock(&g->edges[e].lock);
	barcode_hash_add(&g->edges[e].barcodes, bc);
	pthread_mutex_unlock(&g->edges[e].lock);
}

// static inline void add_barcode_edge_unique(struct asm_graph_t *g, gint_t e, uint64_t bc)
// {
// 	pthread_mutex_lock(&g->edges[e].lock);
// 	barcode_hash_add_unique(&g->edges[e].barcodes, bc);
// 	pthread_mutex_unlock(&g->edges[e].lock);
// }

static inline void add_read_pair_edge(struct asm_graph_t *g, gint_t e, gint_t next_e, uint64_t bc)
{
	if (next_e == e || next_e == g->edges[e].rc_id)
		return;
	pthread_mutex_lock(&g->edges[e].lock);
	gint_t i;
	for (i = 0; i < g->edges[e].n_mate_contigs; ++i)
		if (g->edges[e].mate_contigs[i] == next_e)
			break;
	if (i == g->edges[e].n_mate_contigs) {
		++g->edges[e].n_mate_contigs;
		g->edges[e].mate_contigs = realloc(g->edges[e].mate_contigs,
			g->edges[e].n_mate_contigs * sizeof(gint_t));
		g->edges[e].mate_counts = realloc(g->edges[e].mate_counts,
			g->edges[e].n_mate_contigs * sizeof(gint_t));
		// g->edges[e].mate_barcodes = realloc(g->edges[e].mate_barcodes,
		// 	g->edges[e].n_mate_contigs * sizeof(struct barcode_hash_t));
		g->edges[e].mate_contigs[i] = next_e;
		g->edges[e].mate_counts[i] = 0;
		// barcode_hash_init(g->edges[e].mate_barcodes + i, 4);
	}
	++g->edges[e].mate_counts[i];
	// barcode_hash_add(g->edges[e].mate_barcodes + i, bc);
	pthread_mutex_unlock(&g->edges[e].lock);
}

struct ref_contig_t {
	gint_t e;
	int pos;
	int strand;
};

static inline int count_M_cigar(int n_cigar, uint32_t *cigar)
{
	int i, ret = 0;
	for (i = 0; i < n_cigar; ++i)
		ret += (int)((cigar[i] & 0xf) == 0) * (cigar[i] >> 4);
	return ret;
}

static inline int check_clip_both_end(int n_cigar, uint32_t *cigar)
{
	return (n_cigar > 2 && (cigar[0] & 0xf) == 3 && (cigar[n_cigar - 1] & 0xf) == 3);
}

int is_contain_N(struct read_t *r)
{
	int i;
	for (i = 0; i < r->len; ++i)
		if (nt4_table[(int)r->seq[i]] >= 4)
			return 1;
	return 0;
}

void barcode_read_mapper(struct read_t *r1, struct read_t *r2, uint64_t bc,
						struct bccount_bundle_t *bundle)
{
	int i, k, n1, n2, best_score1, best_score2;
	struct asm_graph_t *g = bundle->g;
	bwaidx_t *idx = bundle->bwa_idx;
	mem_opt_t *opt = bundle->bwa_opt;
	mem_alnreg_v ar1, ar2;
	ar1 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r1->len, r1->seq);
	ar2 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r2->len, r2->seq);
	struct ref_contig_t *p1, *p2;
	p1 = alloca(ar1.n * sizeof(struct ref_contig_t));
	p2 = alloca(ar2.n * sizeof(struct ref_contig_t));
	n1 = n2 = 0;
	best_score1 = best_score2 = -1000 * 1000 * 1000;
	for (i = 0; i < (int)ar1.n; ++i) {
		mem_aln_t a;
		a = mem_reg2aln(opt, idx->bns, idx->pac, r1->len, r1->seq, &ar1.a[i]);
		int aligned = count_M_cigar(a.n_cigar, a.cigar);
		free(a.cigar);
		if (check_clip_both_end(a.n_cigar, a.cigar) ||
			aligned + 20 < r1->len)
			continue;
		gint_t e = atol(idx->bns->anns[a.rid].name);
		if (bundle->aux_build & ASM_BUILD_COVERAGE) {
			if (aligned > g->ksize)
				atomic_add_and_fetch64(&g->edges[e].count, aligned - g->ksize);
		}
		if (ar1.a[i].score > best_score1 && ar1.a[i].score + 5 >= aligned) {
			best_score1 = ar1.a[i].score;
			p1[0] = (struct ref_contig_t){e, (int)a.pos, (int)a.is_rev};
			n1 = 1;
		} else if (ar1.a[i].score == best_score1) {
			for (k = 0; k < n1; ++k) {
				if (p1[k].e == e) {
					p1[k].pos = __min(p1[k].pos, (int)a.pos);
					break;
				}
			}
			if (k == n1)
				p1[n1++] = (struct ref_contig_t){e, (int)a.pos, (int)a.is_rev};
		}
	}
	for (i = 0; i < (int)ar2.n; ++i) {
		mem_aln_t a;
		a = mem_reg2aln(opt, idx->bns, idx->pac, r2->len, r2->seq, &ar2.a[i]);
		int aligned = count_M_cigar(a.n_cigar, a.cigar);
		free(a.cigar);
		if (check_clip_both_end(a.n_cigar, a.cigar) ||
			aligned + 20 < r2->len)
			continue;
		gint_t e = atol(idx->bns->anns[a.rid].name);
		if (bundle->aux_build & ASM_BUILD_COVERAGE) {
			if (aligned > g->ksize)
				atomic_add_and_fetch64(&g->edges[e].count, aligned - g->ksize);
		}
		if (ar2.a[i].score > best_score2 && ar2.a[i].score + 5 >= aligned) {
			best_score2 = ar2.a[i].score;
			p2[0] = (struct ref_contig_t){e, (int)a.pos, (int)a.is_rev};
			n2 = 1;
		} else if (ar2.a[i].score == best_score2) {
			for (k = 0; k < n2; ++k) {
				if (p2[k].e == e) {
					p2[k].pos = __min(p2[k].pos, (int)a.pos);
					break;
				}
			}
			if (k == n2)
				p2[n2++] = (struct ref_contig_t){e, (int)a.pos, (int)a.is_rev};
		}
	}
	if ((bundle->aux_build & ASM_BUILD_BARCODE) && bc != (uint64_t)-1) {
		for (i = 0; i < n1; ++i)
			if (p1[i].pos <= MIN_CONTIG_BARCODE)
				add_barcode_edge(g, p1[i].e, bc);
		for (i = 0; i < n2; ++i)
			if (p2[i].pos <= MIN_CONTIG_BARCODE)
				add_barcode_edge(g, p2[i].e, bc);
		// for (i = 0; i < n1; ++i) {
		// 	if (p1[i].pos > MIN_CONTIG_BARCODE)
		// 		continue;
		// 	for (k = 0; k < n2; ++k) {
		// 		if (p2[k].pos <= MIN_CONTIG_BARCODE &&
		// 			p1[i].e == p2[k].e) {
		// 			add_barcode_edge(g, p1[i].e, bc);
		// 			break;
		// 		}
		// 	}
		// }
	}
	if ((bundle->aux_build & ASM_BUILD_READPAIR) && bc != (uint64_t)-1) {
		for (i = 0; i < n1; ++i) {
			for (k = 0; k < n2; ++k) {
				if (p1[i].e != p2[k].e && p1[i].strand == p2[k].strand
					&& p1[i].pos + p2[k].pos < MAX_PAIR_LEN) {
					add_read_pair_edge(g, p1[i].e, p2[k].e, bc);
					add_read_pair_edge(g, p2[k].e, p1[i].e, bc);
				}
			}
		}
	}
	free(ar1.a);
	free(ar2.a);
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
	gcnt_reads = bundle->n_reads;
	// void (*read_process)(struct read_t *, uint64_t, struct bccount_bundle_t *) = bundle->read_process_func;
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
			barcode = barcode_calculator(&read1, &read2);
			barcode_read_mapper(&read1, &read2, barcode, bundle);

			if (rc1 == READ_END)
				break;
		}
		n_reads = atomic_add_and_fetch64(gcnt_reads, n_reads);
		// __VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

void barcode_start_count(struct opt_proc_t *opt, struct bccount_bundle_t *ske)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int i;
	struct producer_bundle_t *producer_bundles;
	producer_bundles = init_fastq_PE(opt->n_threads, opt->n_files,
						opt->files_1, opt->files_2);

	struct bccount_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct bccount_bundle_t));
	int64_t n_reads = 0;

	pthread_mutex_t lock;
	pthread_mutex_init(&lock, NULL);
	uint64_t hash_sum;
	hash_sum = 0;
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].n_reads = &n_reads;
		worker_bundles[i].g = ske->g;
		worker_bundles[i].aux_build = ske->aux_build;
		worker_bundles[i].bwa_idx = ske->bwa_idx;
		worker_bundles[i].bwa_opt = ske->bwa_opt;
		worker_bundles[i].barcode_calculator = ske->barcode_calculator;
		worker_bundles[i].lock = &lock;
		worker_bundles[i].hash_sum = &hash_sum;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, barcode_buffer_iterator,
				worker_bundles + i);

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	// __VERBOSE("hash sum = %lu\n", hash_sum);

	free_fastq_PE(producer_bundles, opt->n_files);
	free(worker_bundles);

	free(producer_threads);
	free(worker_threads);
}

