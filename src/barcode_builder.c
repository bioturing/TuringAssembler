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
	int64_t *n_reads;
	bwaidx_t *bwa_idx;
	mem_opt_t *bwa_opt;
	uint64_t (*barcode_calculator)(struct read_t *, struct read_t *);
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
	if (len != 18)
		return (uint64_t)-1;
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

double get_barcode_ratio_small(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
}

double get_barcode_ratio(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
}

int test_edge_barcode(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
}

void print_test_barcode_edge(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
}

static inline uint32_t dump_edge_seq(char **seq, uint32_t *m_seq,
					struct asm_edge_t *e, uint32_t max_len)
{
	uint32_t i, j, k, len;
	len = __min(max_len, e->seq_len);
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

void init_barcode_map(struct asm_graph_t *g, const char *path, uint32_t max_len)
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
		fprintf(fp, ">%ld\n", e);
		k = 0;
		while (k < len) {
			l = __min(80, len - k);
			memcpy(buf, seq + k, l);
			buf[l] = '\0';
			fprintf(fp, "%s\n", buf);
			k += l;
		}
		barcode_hash_init(&g->edges[e].barcodes, 4);
	}
	fclose(fp);
	bwa_idx_build(path, path, BWTALGO_AUTO, 25000000);
	__VERBOSE("Done indexing contigs\n");
}

void construct_barcode_map(struct opt_proc_t *opt, struct asm_graph_t *g, int need_count)
{
	char fasta_prefix[1024];
	strcpy(fasta_prefix, opt->out_dir);
	strcat(fasta_prefix, "/contigs.fasta");
	init_barcode_map(g, fasta_prefix, opt->split_len);
	struct bccount_bundle_t ske;
	ske.g = g;
	extern uint64_t (*barcode_calculators[])(struct read_t *, struct read_t *);
	ske.barcode_calculator = barcode_calculators[opt->lib_type];
	ske.bwa_idx = bwa_idx_load(fasta_prefix, BWA_IDX_ALL);
	ske.bwa_opt = mem_opt_init();
	barcode_start_count(opt, &ske);
}

static inline void add_barcode_to_edge(struct asm_graph_t *g, gint_t e, uint64_t bc)
{
	pthread_mutex_lock(&g->edges[e].lock);
	barcode_hash_inc_count(&g->edges[e].barcodes, bc);
	pthread_mutex_unlock(&g->edges[e].lock);
}

void barcode_read_mapper(struct read_t *r1, struct read_t *r2, uint64_t bc,
						struct bccount_bundle_t *bundle)
{
	struct asm_graph_t *g = bundle->g;
	bwaidx_t *idx = bundle->bwa_idx;
	mem_opt_t *opt = bundle->bwa_opt;
	mem_alnreg_v ar1, ar2;
	ar1 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r1->len, r1->seq);
	ar2 = mem_align1(opt, idx->bwt, idx->bns, idx->pac, r2->len, r2->seq);
	int i;
	gint_t prev_e = -1;
	for (i = 0; i < ar1.n; ++i) {
		mem_aln_t a;
		if (ar1.a[i].secondary >= 0)
			continue;
		a = mem_reg2aln(opt, idx->bns, idx->pac, r1->len, r1->seq, &ar1.a[i]);
		gint_t e = atol(idx->bns->anns[a.rid].name);
		if (e != prev_e)
			add_barcode_to_edge(g, e, bc);
		prev_e = e;
		free(a.cigar);
	}
	prev_e = -1;
	for (i = 0; i < ar2.n; ++i) {
		mem_aln_t a;
		if (ar2.a[i].secondary >= 0)
			continue;
		a = mem_reg2aln(opt, idx->bns, idx->pac, r2->len, r2->seq, &ar2.a[i]);
		gint_t e = atol(idx->bns->anns[a.rid].name);
		if (e != prev_e)
			add_barcode_to_edge(g, e, bc);
		prev_e = e;
		free(a.cigar);
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
			// if (barcode != (uint64_t)-1) {
			// 	bcread_iterator(&read1, barcode, bundle);
			// 	bcread_iterator(&read2, barcode, bundle);
			// }

			if (rc1 == READ_END)
				break;
		}
		n_reads = atomic_add_and_fetch64(gcnt_reads, n_reads);
		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

void barcode_start_count(struct opt_proc_t *opt, struct bccount_bundle_t *ske)
{
	__VERBOSE("Start count barcode\n");
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
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].n_reads = &n_reads;
		worker_bundles[i].g = ske->g;
		worker_bundles[i].bwa_idx = ske->bwa_idx;
		worker_bundles[i].bwa_opt = ske->bwa_opt;
		worker_bundles[i].barcode_calculator = ske->barcode_calculator;
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

	free_fastq_PE(producer_bundles, opt->n_files);
	free(worker_bundles);

	free(producer_threads);
	free(worker_threads);
}

