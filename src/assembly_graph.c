#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "fastq_producer.h"
#include "io_utils.h"
#include "k31hash.h"
#include "k63hash.h"
#include "k31_count.h"
#include "k63_count.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

static void dump_bin_seq(char *seq, uint32_t *bin, gint_t len)
{
	gint_t i;
	for (i = 0; i < len; ++i)
		seq[i] = nt4_char[(bin[i >> 4] >> ((i & 15) << 1)) & 3];
	seq[len] = '\0';
}

#define __bin_seq_get_char(seq, l) (((seq)[(l) >> 4] >> (((l) & 15) << 1)) & (uint32_t)0x3)


#define TIPS_THRESHOLD			5.0
#define TIPS_RATIO_THRESHOLD		0.0625

#define NON_TIPS_LEN			250

void remove_tips(struct asm_graph_t *g0, struct asm_graph_t *g);
void remove_tips_topology(struct asm_graph_t *g0, struct asm_graph_t *g);
void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void write_gfa(struct asm_graph_t *g, const char *path);
void dump_fasta(struct asm_graph_t *g, const char *path);
void remove_bubble_simple(struct asm_graph_t *g0, struct asm_graph_t *g);
static int dfs_dead_end(struct asm_graph_t *g, gint_t u,
			gint_t len, gint_t max_len, int ksize);

void k63_build0(struct opt_count_t *opt, int ksize, struct asm_graph_t *g0)
{
	struct k63hash_t *kmer_hash;
	kmer_hash = calloc(1, sizeof(struct k63hash_t));
	__VERBOSE("Estimating kmer\n");
	build_k63_table_lazy(opt, kmer_hash, ksize);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Estimating kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nBuilding assembly graph\n");
	build_asm_graph_from_k63(opt, ksize, kmer_hash, g0);
	__VERBOSE_LOG("kmer_%d_graph_#0", "Number of nodes: %lld\n", ksize,
							(long long)g0->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#0", "Number of edges: %lld\n", ksize,
							(long long)g0->n_e);
}

void k31_build0(struct opt_count_t *opt, int ksize, struct asm_graph_t *g0)
{
	struct k31hash_t *kmer_hash;
	kmer_hash = calloc(1, sizeof(struct k31hash_t));
	__VERBOSE("Estimating kmer\n");
	build_k31_table_lazy(opt, kmer_hash, ksize);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Estimating kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nBuilding assembly graph\n");
	build_asm_graph_from_k31(opt, ksize, kmer_hash, g0);
	__VERBOSE_LOG("kmer_%d_graph_#0", "Number of nodes: %lld\n", ksize,
							(long long)g0->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#0", "Number of edges: %lld\n", ksize,
							(long long)g0->n_e);
}

void graph_convert_process(struct opt_build_t *opt)
{
	char path[1024];
	init_clock();

	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Dump graph from bin archive\n");
	struct asm_graph_t *g;
	g = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g, opt->in_path);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g->ksize);
	test_asm_graph(g);
	snprintf(path, 1024, "%s/graph_k_%d_loaded.gfa", opt->out_dir, g->ksize);
	write_gfa(g, path);
	snprintf(path, 1024, "%s/graph_k_%d_loaded.fasta", opt->out_dir, g->ksize);
	dump_fasta(g, path);
}

void build0_process(struct opt_count_t *opt)
{
	char path[1024];
	init_clock();

	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Building assembly graph using kmer size %d\n", opt->k0);
	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	if (opt->k0 < 32)
		k31_build0(opt, opt->k0, g0);
	else if (opt->k0 > 32 && opt->k0 < 64)
		k63_build0(opt, opt->k0, g0);
	test_asm_graph(g0);
	snprintf(path, 1024, "%s/graph_k_%d_level_0.gfa", opt->out_dir, opt->k0);
	write_gfa(g0, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_0.fasta", opt->out_dir, opt->k0);
	dump_fasta(g0, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_0.bin", opt->out_dir, opt->k0);
	save_asm_graph(g0, path);
}

void build0_1_process(struct opt_build_t *opt)
{
	char path[1024];
	init_clock();

	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_path);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips\n");
	struct asm_graph_t *g1;
	g1 = calloc(1, sizeof(struct asm_graph_t));
	remove_tips(g0, g1);
	__VERBOSE_LOG("kmer_%d_graph_#1", "Number of nodes: %lld\n", g0->ksize,
							(long long)g1->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#1", "Number of edges: %lld\n", g0->ksize,
							(long long)g1->n_e);
	test_asm_graph(g1);
	snprintf(path, 1024, "%s/graph_k_%d_level_1.gfa", opt->out_dir, g0->ksize);
	write_gfa(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_1.fasta", opt->out_dir, g0->ksize);
	dump_fasta(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_1.bin", opt->out_dir, g0->ksize);
	save_asm_graph(g1, path);
}

void build1_2_process(struct opt_build_t *opt)
{
	char path[1024];
	init_clock();

	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_path);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips using graph topology\n");
	struct asm_graph_t *g1;
	g1 = calloc(1, sizeof(struct asm_graph_t));
	remove_tips_topology(g0, g1);
	__VERBOSE_LOG("kmer_%d_graph_#2", "Number of nodes: %lld\n", g0->ksize,
							(long long)g1->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#2", "Number of edges: %lld\n", g0->ksize,
							(long long)g1->n_e);
	test_asm_graph(g1);
	snprintf(path, 1024, "%s/graph_k_%d_level_2.gfa", opt->out_dir, g0->ksize);
	write_gfa(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_2.fasta", opt->out_dir, g0->ksize);
	dump_fasta(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_2.bin", opt->out_dir, g0->ksize);
	save_asm_graph(g1, path);
}

void build2_3_process(struct opt_build_t *opt)
{
	char path[1024];
	init_clock();

	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_path);
	test_asm_graph(g0);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing simple branching\n");
	struct asm_graph_t *g1;
	g1 = calloc(1, sizeof(struct asm_graph_t));
	remove_bubble_simple(g0, g1);
	__VERBOSE_LOG("kmer_%d_graph_#3", "Number of nodes: %lld\n", g0->ksize,
							(long long)g1->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#3", "Number of edges: %lld\n", g0->ksize,
							(long long)g1->n_e);
	test_asm_graph(g1);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.gfa", opt->out_dir, g0->ksize);
	write_gfa(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.fasta", opt->out_dir, g0->ksize);
	dump_fasta(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.bin", opt->out_dir, g0->ksize);
	save_asm_graph(g1, path);
}

void build2_3a_process(struct opt_build_t *opt)
{
	char path[1024];
	init_clock();

	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_path);
	test_asm_graph(g0);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Building barcode information\n");
	construct_barcode_map(g0, opt);
}

void assembly_process(struct opt_count_t *opt)
{
	char path[1024];
	init_clock();

	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Building assembly graph using kmer size %d\n", opt->k0);
	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	if (opt->k0 < 32)
		k31_build0(opt, opt->k0, g0);
	else if (opt->k0 > 32 && opt->k0 < 64)
		k63_build0(opt, opt->k0, g0);
	test_asm_graph(g0);
	snprintf(path, 1024, "%s/graph_k_%d_level_0.gfa", opt->out_dir, opt->k0);
	write_gfa(g0, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_0.fasta", opt->out_dir, opt->k0);
	dump_fasta(g0, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_0.bin", opt->out_dir, opt->k0);
	save_asm_graph(g0, path);

	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips\n");
	struct asm_graph_t *g1;
	g1 = calloc(1, sizeof(struct asm_graph_t));
	remove_tips(g0, g1);
	__VERBOSE_LOG("kmer_%d_graph_#1", "Number of nodes: %lld\n", opt->k0,
							(long long)g1->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#1", "Number of edges: %lld\n", opt->k0,
							(long long)g1->n_e);
	test_asm_graph(g1);
	snprintf(path, 1024, "%s/graph_k_%d_level_1.gfa", opt->out_dir, opt->k0);
	write_gfa(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_1.fasta", opt->out_dir, opt->k0);
	dump_fasta(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_1.bin", opt->out_dir, opt->k0);
	save_asm_graph(g1, path);

	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips #2\n");
	struct asm_graph_t *g2;
	g2 = calloc(1, sizeof(struct asm_graph_t));
	remove_tips_topology(g1, g2);
	__VERBOSE_LOG("kmer_%d_graph_#2", "Number of nodes: %lld\n", opt->k0,
							(long long)g2->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#2", "Number of edges: %lld\n", opt->k0,
							(long long)g2->n_e);
	test_asm_graph(g2);
	snprintf(path, 1024, "%s/graph_k_%d_level_2.gfa", opt->out_dir, opt->k0);
	write_gfa(g2, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_2.fasta", opt->out_dir, opt->k0);
	dump_fasta(g2, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_2.bin", opt->out_dir, opt->k0);
	save_asm_graph(g2, path);

	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing bubbles\n");
	struct asm_graph_t *g3;
	g3 = calloc(1, sizeof(struct asm_graph_t));
	remove_bubble_simple(g2, g3);
	__VERBOSE_LOG("kmer_%d_graph_#3", "Number of nodes: %lld\n", opt->k0,
							(long long)g3->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#3", "Number of edges: %lld\n", opt->k0,
							(long long)g3->n_e);
	test_asm_graph(g3);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.gfa", opt->out_dir, opt->k0);
	write_gfa(g3, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.fasta", opt->out_dir, opt->k0);
	dump_fasta(g3, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.bin", opt->out_dir, opt->k0);
	save_asm_graph(g3, path);
}

void k63_process(struct opt_count_t *opt)
{
	char path[1024];
	init_clock();

	struct k63hash_t *kmer_hash;
	kmer_hash = calloc(1, sizeof(struct k63hash_t));
	__VERBOSE("Estimating kmer\n");
	build_k63_table_lazy(opt, kmer_hash, opt->k0);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Estimating kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nBuilding assembly graph\n");
	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	build_asm_graph_from_k63(opt, opt->k0, kmer_hash, g0);
	__VERBOSE_LOG("Graph #1", "Number of nodes: %lld\n", (long long)g0->n_v);
	__VERBOSE_LOG("Graph #1", "Number of edges: %lld\n", (long long)g0->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_0.gfa");
	write_gfa(g0, path);
	test_asm_graph(g0);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_0.fasta");
	dump_fasta(g0, path);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_0.bin");
	save_asm_graph(g0, path);

	__VERBOSE("\nRemoving tips\n");
	struct asm_graph_t *g1;
	g1 = calloc(1, sizeof(struct asm_graph_t));
	remove_tips(g0, g1);
	__VERBOSE_LOG("Graph #2", "Number of nodes: %lld\n", (long long)g1->n_v);
	__VERBOSE_LOG("Graph #2", "Number of edges: %lld\n", (long long)g1->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_1.gfa");
	write_gfa(g1, path);
	test_asm_graph(g1);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_1.fasta");
	dump_fasta(g1, path);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_1.bin");
	save_asm_graph(g1, path);


	__VERBOSE("\nRemoving tips #2\n");
	struct asm_graph_t *g2;
	g2 = calloc(1, sizeof(struct asm_graph_t));
	remove_tips_topology(g1, g2);
	__VERBOSE_LOG("Graph #3", "Number of nodes: %lld\n", (long long)g2->n_v);
	__VERBOSE_LOG("Graph #3", "Number of edges: %lld\n", (long long)g2->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_2.gfa");
	write_gfa(g2, path);
	test_asm_graph(g2);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_2.fasta");
	dump_fasta(g2, path);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_2.bin");
	save_asm_graph(g2, path);

	__VERBOSE("\nRemoving bubble\n");
	struct asm_graph_t *g3;
	g3 = calloc(1, sizeof(struct asm_graph_t));
	remove_bubble_simple(g2, g3);
	__VERBOSE_LOG("Graph #4", "Number of nodes: %lld\n", (long long)g3->n_v);
	__VERBOSE_LOG("Graph #4", "Number of edges: %lld\n", (long long)g3->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_3.gfa");
	write_gfa(g3, path);
	test_asm_graph(g3);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_3.fasta");
	dump_fasta(g3, path);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_3.bin");
	save_asm_graph(g3, path);
}

void k31_process(struct opt_count_t *opt)
{
	char path[1024];
	init_clock();

	struct k31hash_t *kmer_hash;
	kmer_hash = calloc(1, sizeof(struct k31hash_t));
	__VERBOSE("Estimating kmer\n");
	build_k31_table_lazy(opt, kmer_hash, opt->k0);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Estimating kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nBuilding assembly graph\n");
	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	build_asm_graph_from_k31(opt, opt->k0, kmer_hash, g0);
	__VERBOSE_LOG("Graph #1", "Number of nodes: %lld\n", (long long)g0->n_v);
	__VERBOSE_LOG("Graph #1", "Number of edges: %lld\n", (long long)g0->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_0.gfa");
	write_gfa(g0, path);
	test_asm_graph(g0);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_0.fasta");
	dump_fasta(g0, path);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_0.bin");
	save_asm_graph(g0, path);

	__VERBOSE("\nRemoving tips\n");
	struct asm_graph_t *g1;
	g1 = calloc(1, sizeof(struct asm_graph_t));
	remove_tips(g0, g1);
	__VERBOSE_LOG("Graph #2", "Number of nodes: %lld\n", (long long)g1->n_v);
	__VERBOSE_LOG("Graph #2", "Number of edges: %lld\n", (long long)g1->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_1.gfa");
	write_gfa(g1, path);
	test_asm_graph(g1);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_1.fasta");
	dump_fasta(g1, path);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_1.bin");
	save_asm_graph(g1, path);

	__VERBOSE("\nRemoving tips #2\n");
	struct asm_graph_t *g2;
	g2 = calloc(1, sizeof(struct asm_graph_t));
	remove_tips_topology(g1, g2);
	__VERBOSE_LOG("Graph #3", "Number of nodes: %lld\n", (long long)g2->n_v);
	__VERBOSE_LOG("Graph #3", "Number of edges: %lld\n", (long long)g2->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_2.gfa");
	write_gfa(g2, path);
	test_asm_graph(g2);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_2.fasta");
	dump_fasta(g2, path);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_2.bin");
	save_asm_graph(g2, path);

	__VERBOSE("\nRemoving bubble\n");
	struct asm_graph_t *g3;
	g3 = calloc(1, sizeof(struct asm_graph_t));
	remove_bubble_simple(g2, g3);
	__VERBOSE_LOG("Graph #4", "Number of nodes: %lld\n", (long long)g3->n_v);
	__VERBOSE_LOG("Graph #4", "Number of edges: %lld\n", (long long)g3->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_3.gfa");
	write_gfa(g3, path);
	test_asm_graph(g3);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_3.fasta");
	dump_fasta(g3, path);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_3.bin");
	save_asm_graph(g3, path);
}

double get_uni_genome_coverage(struct asm_graph_t *g)
{
	gint_t max_len, e;
	double ret_cov = 0.0;
	max_len = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].seq_len > max_len) {
			max_len = g->edges[e].seq_len;
			ret_cov = g->edges[e].count * 1.0 /
				(g->edges[e].seq_len - g->ksize);
		}
	}
	return ret_cov;
}

static uint32_t *append_bin_seq(uint32_t *dst, gint_t dlen, uint32_t *src,
					gint_t slen, int skip)
{
	{
		gint_t i;
		for (i = 0; i < skip; ++i) {
			uint32_t c1, c2;
			c1 = __bin_seq_get_char(src, i);
			c2 = __bin_seq_get_char(dst, dlen - (skip - i));
			if (c1 != c2) {
				char *seq1, *seq2;
				seq1 = malloc(dlen + 1);
				dump_bin_seq(seq1, dst, dlen);
				fprintf(stderr, "dst = %s\n", seq1);
				seq2 = malloc(slen + skip + 1);
				dump_bin_seq(seq2, src, slen + skip);
				fprintf(stderr, "src = %s\n", seq2);
				assert(0 && "Fail when append");
			}
		}
	}
	uint32_t *new_ptr;
	gint_t cur_m, m, i, k;
	cur_m = (dlen + 15) >> 4;
	m = (dlen + slen + 15) >> 4;
	if (m > cur_m) {
		new_ptr = realloc(dst, m * sizeof(uint32_t));
		memset(new_ptr + cur_m, 0, (m - cur_m) * sizeof(uint32_t));
	} else {
		new_ptr = dst;
	}
	for (i = skip; i < skip + slen; ++i) {
		k = dlen + (i - skip);
		new_ptr[k >> 4] |= ((src[i >> 4] >> ((i & 15) << 1)) & (uint32_t)3) << ((k & 15) << 1);
		assert(__bin_seq_get_char(new_ptr, k) == ((src[i >> 4] >> ((i & 15) << 1)) & (uint32_t)3));
	}
	return new_ptr;
}

static uint32_t *asm_clone_binseq(uint32_t *src, gint_t len)
{
	uint32_t *ret;
	ret = calloc((len + 15) >> 4, sizeof(uint32_t));
	memcpy(ret, src, ((len + 15) >> 4) * sizeof(uint32_t));
	return ret;
}

static gint_t asm_find_edge_index(gint_t *adj, gint_t deg, gint_t id)
{
	gint_t i, ret;
	ret = -1;
	for (i = 0; i < deg; ++i) {
		if (adj[i] == id)
			ret = i;
	}
	return ret;
}

static gint_t asm_only_positive_edge(gint_t *adj, int deg)
{
	int i, ret;
	ret = -1;
	for (i = 0; i < deg; ++i) {
		if (adj[i] != -1) {
			if (ret == -1)
				ret = adj[i];
			else
				ret = -2;
		}
	}
	return ret;
}

int asm_is_edge_rc(uint32_t *seq1, gint_t l1, uint32_t *seq2, gint_t l2)
{
	if (l1 != l2)
		return 0;
	gint_t i, k;
	uint32_t c1, c2;
	for (i = 0; i < l1; ++i) {
		k = l1 - i - 1;
		c1 = __bin_seq_get_char(seq1, i);
		c2 = __bin_seq_get_char(seq2, k);
		if (c1 != (c2 ^ 3)) {
			return 0;
		}
	}
	return 1;
}

void asm_edge_cc(struct asm_graph_t *g, gint_t *id_edge, gint_t **ret_size)
{
	memset(id_edge, 255, g->n_e * sizeof(gint_t));
	gint_t m_cc = 0x10000;
	gint_t *size = malloc(m_cc * sizeof(gint_t));
	gint_t n_cc = 0;
	gint_t k, l, r;
	gint_t *q = malloc(g->n_e * sizeof(gint_t));

	for (k = 0; k < g->n_e; ++k) {
		if (id_edge[k] != -1)
			continue;
		id_edge[k] = id_edge[g->edges[k].rc_id] = n_cc;
		gint_t cur_size = 0;
		l = r = 0;
		q[0] = k;
		while (l <= r) {
			gint_t e = q[l++];
			gint_t e_rc = g->edges[e].rc_id;
			cur_size += 2 * (g->edges[e].seq_len - g->ksize);

			gint_t u, c;
			u = g->edges[e].target;
			if (g->nodes[u].deg == 0)
				cur_size += g->ksize;
			for (c = 0; c < g->nodes[u].deg; ++c) {
				gint_t next_e = g->nodes[u].adj[c];
				gint_t next_e_rc = g->edges[next_e].rc_id;
				if (id_edge[next_e] == -1) {
					id_edge[next_e] = id_edge[next_e_rc] = n_cc;
					q[++r] = next_e;
				}
			}

			u = g->edges[e_rc].target;
			if (g->nodes[u].deg == 0)
				cur_size += g->ksize;
			for (c = 0; c < g->nodes[u].deg; ++c) {
				gint_t next_e = g->nodes[u].adj[c];
				gint_t next_e_rc = g->edges[next_e].rc_id;
				if (id_edge[next_e] == -1) {
					id_edge[next_e] = id_edge[next_e_rc] = n_cc;
					q[++r] = next_e;
				}
			}
		}
		if (n_cc + 1 == m_cc) {
			m_cc <<= 1;
			size = realloc(size, m_cc * sizeof(gint_t));
		}
		size[n_cc++] = cur_size / 2;
	}
	size = realloc(size, n_cc * sizeof(gint_t));
	*ret_size = size;
	free(q);
}

void dump_fasta(struct asm_graph_t *g, const char *path)
{
	gint_t *id_edge, *cc_size;
	id_edge = malloc(g->n_e * sizeof(gint_t));
	cc_size = NULL;
	asm_edge_cc(g, id_edge, &cc_size);

	FILE *fp = xfopen(path, "w");
	char *seq = NULL;
	int seq_len = 0;
	char *buf = alloca(81);
	gint_t e, e_rc;
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		gint_t cc_id = id_edge[e];
		if (cc_size[cc_id] < 250 || g->edges[e].seq_len < 100)
			continue;
		if (seq_len < g->edges[e].seq_len + 1) {
			seq_len = g->edges[e].seq_len + 1;
			seq = realloc(seq, seq_len);
		}
		dump_bin_seq(seq, g->edges[e].seq, g->edges[e].seq_len);
		fprintf(fp, ">SEQ_%lld_length_%lld_count_%llu\n", (long long)e,
			(long long)g->edges[e].seq_len, (long long unsigned)g->edges[e].count);
		gint_t k = 0;
		while (k < g->edges[e].seq_len) {
			gint_t l = __min(80, g->edges[e].seq_len - k);
			memcpy(buf, seq + k, l);
			buf[l] = '\0';
			fprintf(fp, "%s\n", buf);
			k += l;
		}
	}
	fclose(fp);

	free(id_edge);
	free(cc_size);
}


void write_gfa(struct asm_graph_t *g, const char *path)
{
	gint_t *id_edge, *cc_size;
	id_edge = malloc(g->n_e * sizeof(gint_t));
	cc_size = NULL;
	asm_edge_cc(g, id_edge, &cc_size);

	FILE *fp = xfopen(path, "w");
	char *seq = NULL;
	int seq_len = 0;
	gint_t e, e_rc;
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		gint_t cc_id = id_edge[e];
		if (cc_size[cc_id] < 250)
			continue;
		if (seq_len < g->edges[e].seq_len + 1) {
			seq_len = g->edges[e].seq_len + 1;
			seq = realloc(seq, seq_len);
		}
		dump_bin_seq(seq, g->edges[e].seq, g->edges[e].seq_len);
		double cov = g->edges[e].count * 1.0 / (g->edges[e].seq_len - g->ksize);
		uint64_t fake_count = (uint64_t)(cov * g->edges[e].seq_len);
		/* print fake count for correct coverage display on Bandage */
		fprintf(fp, "S\t%lld\t%s\tKC:i:%llu\n", (long long)e, seq,
					(long long unsigned)fake_count);
	}
	for (e = 0; e < g->n_e; ++e) {
		gint_t cc_id = id_edge[e];
		if (cc_size[cc_id] < 250)
			continue;
		e_rc = g->edges[e].rc_id;
		gint_t pe, next_pe;
		char ce, next_ce;
		if (e > e_rc) {
			pe = e_rc;
			ce = '-';
		} else {
			pe = e;
			ce = '+';
		}
		gint_t n = g->edges[e].target;
		gint_t k;
		for (k = 0; k < g->nodes[n].deg; ++k) {
			gint_t next_e, next_e_rc;
			next_e = g->nodes[n].adj[k];
			next_e_rc = g->edges[next_e].rc_id;
			if (next_e > next_e_rc) {
				next_pe = next_e_rc;
				next_ce = '-';
			} else {
				next_pe = next_e;
				next_ce = '+';
			}
			fprintf(fp, "L\t%lld\t%c\t%lld\t%c\t%dM\n",
				(long long)pe, ce, (long long)next_pe, next_ce,
				g->ksize);
		}
	}
	fclose(fp);
	free(id_edge);
	free(cc_size);
}

void remove_bubble_simple(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	/* "Don't you want a balloon?"
	 * 1g = 1 genome walk
	 *                   0.8g
	 *                +--------+
	 *               /          \
	 *       1g     /            v    1g
	 * o---------->U             V--------->o
	 *              \            ^
	 *               \   0.3g   /
	 *                +--------+
	 */
	double uni_cov = get_uni_genome_coverage(g0);
	__VERBOSE("1 genome walk coverage: %lf\n", uni_cov);
	gint_t cnt_rm = 0;
	gint_t u;
	for (u = 0; u < g0->n_v; ++u) {
		if (u % 1000 == 0)
			fprintf(stderr, "\rRemoving bubble node %ld", u);
		gint_t ctg;
		double cov;
		/* check topology of node u */
		gint_t u_rc = g0->nodes[u].rc_id;
		if (g0->nodes[u_rc].deg != 1 || g0->nodes[u].deg != 2)
			continue;
		/* check if previous contig is uni genome walk */
		ctg = g0->nodes[u_rc].adj[0];
		cov = g0->edges[ctg].count * 1.0 /
			(g0->edges[ctg].seq_len - g0->ksize);
		if (cov / uni_cov < 0.75 || cov / uni_cov > 1.25)
			continue;

		/* check topology of 2 small contigs (end at same node) */
		gint_t e0, e1;
		e0 = g0->nodes[u].adj[0];
		e1 = g0->nodes[u].adj[1];
		if (e0 == g0->edges[e1].rc_id)
			continue;
		if (g0->edges[e0].target != g0->edges[e1].target)
			continue;
		/* check topology of node v */
		gint_t v, v_rc;
		v = g0->edges[e0].target;
		v_rc = g0->nodes[v].rc_id;
		if (g0->nodes[v_rc].deg != 2 || g0->nodes[v].deg != 1)
			continue;
		/* check if next contig is uni genome walk */
		ctg = g0->nodes[v].adj[0];
		cov = g0->edges[ctg].count * 1.0 /
			(g0->edges[ctg].seq_len - g0->ksize);
		if (cov / uni_cov < 0.75 || cov / uni_cov > 1.25)
			continue;

		gint_t e_rc0, e_rc1;
		e_rc0 = g0->nodes[v_rc].adj[0];
		e_rc1 = g0->nodes[v_rc].adj[1];
		assert(g0->edges[e_rc0].target == g0->edges[e_rc1].target);
		/* keep the longer edge */
		gint_t e;
		if (g0->edges[e0].seq_len > g0->edges[e1].seq_len)
			e = e0;
		else
			e = e1;
		g0->nodes[u].adj[0] = e;
		g0->nodes[u].deg = 1;
		g0->nodes[u].adj = realloc(g0->nodes[u].adj, sizeof(gint_t));
		g0->edges[e].count = g0->edges[e0].count + g0->edges[e1].count;

		assert(g0->edges[e].rc_id == e_rc0 || g0->edges[e].rc_id == e_rc1);
		g0->nodes[v_rc].adj[0] = g0->edges[e].rc_id;
		g0->nodes[v_rc].deg = 1;
		g0->nodes[v_rc].adj = realloc(g0->nodes[v_rc].adj, sizeof(gint_t));
		g0->edges[g0->edges[e].rc_id].count = g0->edges[e_rc0].count + g0->edges[e_rc1].count;
		++cnt_rm;
	}
	__VERBOSE("\nNumber of removed edges: %ld\n", cnt_rm);
	test2_asm_graph(g0);
	asm_condense(g0, g);
}

struct edge_cov_t {
	double cov;
	gint_t e;
};

static inline void ec_sort(struct edge_cov_t *b, struct edge_cov_t *e)
{
	struct edge_cov_t *i, *j, t;
	for (i = b + 1; i < e; ++i) {
		if (i->cov < (i - 1)->cov) {
			t = *i;
			for (j = i; j > b && t.cov < (j - 1)->cov; --j)
				*j = *(j - 1);
			*j = t;
		}
	}
}

void remove_tips_topology(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	/*             o
	 *             ^      o
	 *       o<--+ |      ^
	 *            \|     /         TIP
	 *        o<---o--->o--->o
	 *             ^
	 *            /
	 * o-------->o-------->o------>o------>o
	 */
	gint_t u, e;
	uint32_t *flag = calloc((g0->n_e + 31) >> 5, sizeof(uint32_t));
	struct edge_cov_t *tmp = NULL;
	gint_t tmp_len = 0;
	gint_t count_rm = 0;
	for (u = 0; u < g0->n_v; ++u) {
		if (u % 1000 == 0)
			fprintf(stderr, "\rRemoving dead-end node %ld", u);
		gint_t c, deg, cnt;
		deg = g0->nodes[u].deg;
		if (deg > tmp_len) {
			tmp_len = deg;
			tmp = realloc(tmp, tmp_len * sizeof(struct edge_cov_t));
		}
		cnt = 0;
		for (c = 0; c < g0->nodes[u].deg; ++c) {
			gint_t e_id;
			e_id = g0->nodes[u].adj[c];
			int ret = dfs_dead_end(g0, g0->edges[e_id].target,
					g0->ksize, NON_TIPS_LEN, g0->ksize);
			if (ret) {
				tmp[cnt].e = e_id;
				tmp[cnt].cov = g0->edges[e_id].count * 1.0 /
					(g0->edges[e_id].seq_len - g0->ksize);
				++cnt;
			}
		}
		if (!cnt)
			continue;
		ec_sort(tmp, tmp + cnt);
		for (c = 0; c < cnt && deg > 1; ++c) {
			gint_t e_id, e_rc;
			e_id = tmp[c].e;
			e_rc = g0->edges[e_id].rc_id;
			flag[e_id >> 5] |= (uint32_t)1 << (e_id & 31);
			flag[e_rc >> 5] |= (uint32_t)1 << (e_rc & 31);
			--deg;
			++count_rm;
		}
	}
	fprintf(stderr, "\nNumber of removed edges: %ld\n", count_rm);
	free(tmp);
	for (e = 0; e < g0->n_e; ++e) {
		if (!((flag[e >> 5] >> (e & 31)) & 1))
			continue;
		u = g0->edges[e].source;
		gint_t c = asm_find_edge_index(g0->nodes[u].adj, g0->nodes[u].deg, e);
		assert(c != -1);
		g0->nodes[u].adj[c] = -1;
	}
	for (u = 0; u < g0->n_v; ++u) {
		gint_t c, deg;
		deg = g0->nodes[u].deg;
		g0->nodes[u].deg = 0;
		for (c = 0; c < deg; ++c) {
			gint_t e_id = g0->nodes[u].adj[c];
			if (e_id != -1)
				g0->nodes[u].adj[g0->nodes[u].deg++] = e_id;
		}
		g0->nodes[u].adj = realloc(g0->nodes[u].adj,
					g0->nodes[u].deg * sizeof(gint_t));
	}
	free(flag);

	asm_condense(g0, g);
}

void remove_tips(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	/*
	 *                      o
	 *                      ^
	 *                     /
	 *              cov~1 /
	 *                   /
	 *       cov~10     /   cov~10
	 * o-------------->o------------->o---------->o
	 */
	gint_t u;
	for (u = 0; u < g0->n_v; ++u) {
		double max_cov = 0.0, cov;
		gint_t deg = g0->nodes[u].deg, c;
		gint_t e_id, e_rc;
		for (c = 0; c < deg; ++c) {
			e_id = g0->nodes[u].adj[c];
			if (e_id == -1)
				continue;
			cov = g0->edges[e_id].count * 1.0 /
					(g0->edges[e_id].seq_len - g0->ksize);
			max_cov = __max(max_cov, cov);
		}
		for (c = 0; c < deg; ++c) {
			e_id = g0->nodes[u].adj[c];
			if (e_id == -1)
				continue;
			cov = g0->edges[e_id].count * 1.0 /
					(g0->edges[e_id].seq_len - g0->ksize);
			if (cov / max_cov < TIPS_RATIO_THRESHOLD) {
				gint_t v, v_rc;
				v = g0->edges[e_id].target;
				v_rc = g0->nodes[v].rc_id;
				e_rc = g0->edges[e_id].rc_id;
				gint_t j = asm_find_edge_index(g0->nodes[v_rc].adj,
						g0->nodes[v_rc].deg, e_rc);
				assert(j != -1);
				/* disconnect edge */
				g0->nodes[u].adj[c] = -1;
				g0->nodes[v_rc].adj[j] = -1;
			}
		}
	}
	for (u = 0; u < g0->n_v; ++u) {
		gint_t c, deg;
		deg = g0->nodes[u].deg;
		g0->nodes[u].deg = 0;
		for (c = 0; c < deg; ++c) {
			gint_t e_id = g0->nodes[u].adj[c];
			if (e_id != -1)
				g0->nodes[u].adj[g0->nodes[u].deg++] = e_id;
		}
		g0->nodes[u].adj = realloc(g0->nodes[u].adj,
					g0->nodes[u].deg * sizeof(gint_t));
	}

	asm_condense(g0, g);
}

void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	gint_t *node_id, *edge_id;
	gint_t n_v, n_e;
	gint_t u;
	node_id = malloc(g0->n_v * sizeof(gint_t));
	edge_id = malloc(g0->n_e * sizeof(gint_t));
	memset(node_id, 255, g0->n_v * sizeof(gint_t));
	memset(edge_id, 255, g0->n_e * sizeof(gint_t));

	n_v = n_e = 0;
	for (u = 0; u < g0->n_v; ++u) {
		int deg_fw = g0->nodes[u].deg;
		int deg_rv = g0->nodes[g0->nodes[u].rc_id].deg;
		if ((deg_fw == 1 && deg_rv == 1) || deg_fw + deg_rv == 0)
			continue;
		node_id[u] = n_v++;
		n_e += deg_fw;
	}

	struct asm_node_t *nodes = calloc(n_v, sizeof(struct asm_node_t));
	struct asm_edge_t *edges = calloc(n_e, sizeof(struct asm_edge_t));
	n_e = 0;
	for (u = 0; u < g0->n_v; ++u) {
		gint_t new_id = node_id[u];
		if (new_id == -1)
			continue;
		assert(node_id[g0->nodes[u].rc_id] != -1);
		gint_t new_id_rc = node_id[g0->nodes[u].rc_id];
		nodes[new_id].rc_id = new_id_rc;
		nodes[new_id_rc].rc_id = new_id;
		nodes[new_id].adj = malloc(g0->nodes[u].deg * sizeof(gint_t));
		nodes[new_id].deg = 0;
		gint_t c;
		for (c = 0; c < g0->nodes[u].deg; ++c) {
			gint_t e_id = g0->nodes[u].adj[c], e_rc, n_id;
			assert(e_id != -1);
			edges[n_e].seq = asm_clone_binseq(g0->edges[e_id].seq,
						g0->edges[e_id].seq_len);
			edges[n_e].seq_len = g0->edges[e_id].seq_len;
			edges[n_e].count = g0->edges[e_id].count;
			edges[n_e].source = new_id;
			edges[n_e].rc_id = -1;

			do {
				assert(edge_id[e_id] == -1);
				edge_id[e_id] = n_e;
				n_id = g0->edges[e_id].target;
				if (node_id[n_id] == -1) { /* middle node */
					assert(g0->nodes[n_id].deg == 1);
					e_id = g0->nodes[n_id].adj[0];
					assert(e_id != -1);
					edges[n_e].seq =
						append_bin_seq(edges[n_e].seq,
							edges[n_e].seq_len,
							g0->edges[e_id].seq,
							g0->edges[e_id].seq_len
								- g0->ksize,
							g0->ksize);
					edges[n_e].seq_len +=
						g0->edges[e_id].seq_len - g0->ksize;
					edges[n_e].count += g0->edges[e_id].count;
				} else {
					e_id = -1;
				}
			} while (e_id != -1);
			assert(node_id[n_id] != -1);
			edges[n_e].target = node_id[n_id];
			nodes[new_id].adj[nodes[new_id].deg++] = n_e;
			e_id = g0->nodes[u].adj[c];
			e_rc = g0->edges[e_id].rc_id;
			if (edge_id[e_rc] != -1) {
				edges[n_e].rc_id = edge_id[e_rc];
				edges[edge_id[e_rc]].rc_id = n_e;
			}
			++n_e;
		}
		nodes[new_id].adj = realloc(nodes[new_id].adj, nodes[new_id].deg * sizeof(gint_t));
	}

	for (gint_t e = 0; e < g->n_e; ++e) {
		if (edge_id[e] != -1) {
			gint_t e_rc = g->edges[e].rc_id;
			gint_t new_e = edge_id[e];
			if (edges[new_e].rc_id == -1) {
				gint_t new_e_rc = edge_id[e_rc];
				fprintf(stderr, "%ld <-> %ld\n", e, e_rc);
				fprintf(stderr, "%ld <-> %ld\n", new_e, new_e_rc);
				assert(0 && "Edge not link reverse complement");
			}
		}
	}

	free(node_id);
	free(edge_id);
	g->ksize = g0->ksize;
	g->n_v = n_v;
	g->n_e = n_e;
	g->nodes = nodes;
	g->edges = edges;
}

static int dfs_dead_end(struct asm_graph_t *g, gint_t u,
			gint_t len, gint_t max_len, int ksize)
{
	if (len >= max_len)
		return 0;
	gint_t j, e;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		int k = dfs_dead_end(g, g->edges[e].target,
			len + g->edges[e].seq_len - ksize, max_len, ksize);
		if (k == 0)
			return 0;
	}
	return 1;
}

void test2_asm_graph(struct asm_graph_t *g)
{
	gint_t u, e;
	for (u = 0; u < g->n_v; ++u) {
		gint_t j;
		for (j = 0; j < g->nodes[u].deg; ++j) {
			e = g->nodes[u].adj[j];
			gint_t e_rc = g->edges[e].rc_id;
			if (!asm_is_edge_rc(g->edges[e].seq, g->edges[e].seq_len,
					g->edges[e_rc].seq, g->edges[e_rc].seq_len)) {
				fprintf(stderr, "seq_len = %ld; rc_seq_len = %ld\n",
					g->edges[e].seq_len, g->edges[e_rc].seq_len);
				assert(g->edges[e].seq_len == g->edges[e_rc].seq_len);
				char *seq = malloc(g->edges[e].seq_len + 1);
				dump_bin_seq(seq, g->edges[e].seq, g->edges[e].seq_len);
				fprintf(stderr, "seq    = %s\n", seq);
				dump_bin_seq(seq, g->edges[e_rc].seq, g->edges[e_rc].seq_len);
				fprintf(stderr, "seq_rc = %s\n", seq);
				assert(0 && "Smart error");
			}
		}
	}
}

void test_asm_graph(struct asm_graph_t *g)
{
	gint_t u, e;
	for (u = 0; u < g->n_v; ++u) {
		gint_t j;
		for (j = 0; j < g->nodes[u].deg; ++j) {
			gint_t e_id = g->nodes[u].adj[j];
			if (g->edges[e_id].source != u) {
				fprintf(stderr, "node = %ld; edges = (%ld -> %ld)\n",
					u, g->edges[e_id].source, g->edges[e_id].target);
				assert(0 && "Fail node reverse complement id");
			}
		}
	}
	for (e = 0; e < g->n_e; ++e) {
		gint_t e_rc = g->edges[e].rc_id;
		// if (e == 12584 || e_rc == 12584)
		// 	fprintf(stderr, "e = %ld; e_rc = %ld; source = %ld; target = %ld\n",
		// 		e, e_rc, g->edges[e].source, g->edges[e].target);
		// if (e == 12586 || e_rc == 12586)
		// 	fprintf(stderr, "e = %ld; e_rc = %ld; source = %ld; target = %ld\n",
		// 		e, e_rc, g->edges[e].source, g->edges[e].target);
		// if (e == 86376 || e_rc == 86376)
		// 	fprintf(stderr, "e = %ld; e_rc = %ld; source = %ld; target = %ld\n",
		// 		e, e_rc, g->edges[e].source, g->edges[e].target);

		if (e != g->edges[e_rc].rc_id) {
			fprintf(stderr, "edge = (%ld -> %ld)\n",
				g->edges[e].source, g->edges[e].target);
			fprintf(stderr, "edge_rc = (%ld -> %ld)\n",
				g->edges[e_rc].source, g->edges[e_rc].target);
			assert(0 && "Fail edge reverse complement id");
		}
	}
	for (u = 0; u < g->n_v; ++u) {
		int c;
		for (c = 0; c < g->ksize; ++c) {
			int prev_c = -1;
			gint_t k;
			for (k = 0; k < g->nodes[u].deg; ++k) {
				gint_t e = g->nodes[u].adj[k];
				int ck = __bin_seq_get_char(g->edges[e].seq, c);
				if (ck != prev_c && prev_c != -1) {
					assert(0 && "Mistake in sequence extract");
				}
				prev_c = ck;
			}
		}
	}
	for (u = 0; u < g->n_v; ++u) {
		gint_t j, k;
		for (j = 0; j < g->nodes[u].deg; ++j) {
			gint_t e_id = g->nodes[u].adj[j];
			gint_t v = g->edges[e_id].target;
			gint_t v_rc = g->nodes[v].rc_id;
			int found = 0;
			for (k = 0; k < g->nodes[v_rc].deg; ++k) {
				gint_t ek = g->nodes[v_rc].adj[k];
				if (g->edges[ek].target == g->nodes[u].rc_id) {
					found = 1;
					break;
				}
			}
			if (!found) {
				fprintf(stderr, "Corrupt node: %ld; u_rc = %ld\n", u, g->nodes[u].rc_id);
				fprintf(stderr, "Corrupt at edge number: %ld; Info (%ld -> %ld) %ld\n",
					e_id, g->edges[e_id].source, g->edges[e_id].target, g->edges[e_id].rc_id);
				fprintf(stderr, "v_rc = %ld\n", v_rc);
				for (k = 0; k < g->nodes[v_rc].deg; ++k) {
					gint_t ek = g->nodes[v_rc].adj[k];
					fprintf(stderr, "\tedges: %ld; Info (%ld -> %ld) %ld\n",
						ek, g->edges[ek].source, g->edges[ek].target, g->edges[ek].rc_id);
				}
				assert(0);
			}
		}
	}
	for (e = 0; e < g->n_e; ++e) {
		gint_t e_rc = g->edges[e].rc_id;
		if (!asm_is_edge_rc(g->edges[e].seq, g->edges[e].seq_len,
				g->edges[e_rc].seq, g->edges[e_rc].seq_len)) {
			fprintf(stderr, "seq_len = %ld; rc_seq_len = %ld\n",
				g->edges[e].seq_len, g->edges[e_rc].seq_len);
			assert(g->edges[e].seq_len == g->edges[e_rc].seq_len);
			char *seq = malloc(g->edges[e].seq_len + 1);
			dump_bin_seq(seq, g->edges[e].seq, g->edges[e].seq_len);
			fprintf(stderr, "seq    = %s\n", seq);
			dump_bin_seq(seq, g->edges[e_rc].seq, g->edges[e_rc].seq_len);
			fprintf(stderr, "seq_rc = %s\n", seq);
			assert(0 && "Stupid error");
		}
		assert(e_rc != -1);
		gint_t src, tgt, src_rc, tgt_rc;
		src = g->edges[e].source;
		tgt = g->edges[e].target;
		src_rc = g->edges[e_rc].source;
		tgt_rc = g->edges[e_rc].target;
		if (!(src == g->nodes[tgt_rc].rc_id && g->nodes[src].rc_id == tgt_rc
				&& tgt == g->nodes[src_rc].rc_id && g->nodes[tgt].rc_id == src_rc)) {
			fprintf(stderr, "source = %ld; target = %ld; rc[source] = %ld; rc[target] = %ld\n",
				src, tgt, g->nodes[src].rc_id, g->nodes[tgt].rc_id);
			fprintf(stderr, "source_rc = %ld; target_rc = %ld; rc[source_rc] = %ld; rc[target_rc] = %ld\n",
				src_rc, tgt_rc, g->nodes[src_rc].rc_id, g->nodes[tgt_rc].rc_id);
			assert(0);
		}
		assert(src == g->nodes[tgt_rc].rc_id && g->nodes[src].rc_id == tgt_rc);
		assert(tgt == g->nodes[src_rc].rc_id && g->nodes[tgt].rc_id == src_rc);
	}
}

void save_asm_graph(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "wb");
	char *sig = "asmg";
	xfwrite(sig, 4, 1, fp);
	xfwrite(&g->ksize, sizeof(int), 1, fp);
	xfwrite(&g->n_v, sizeof(gint_t), 1, fp);
	xfwrite(&g->n_e, sizeof(gint_t), 1, fp);
	gint_t u, e;
	for (u = 0; u < g->n_v; ++u) {
		xfwrite(&g->nodes[u].rc_id, sizeof(gint_t), 1, fp);
		xfwrite(&g->nodes[u].deg, sizeof(gint_t), 1, fp);
		if (g->nodes[u].deg)
			xfwrite(g->nodes[u].adj, sizeof(gint_t), g->nodes[u].deg, fp);
	}
	for (e = 0; e < g->n_e; ++e) {
		xfwrite(&g->edges[e].source, sizeof(gint_t), 1, fp);
		xfwrite(&g->edges[e].target, sizeof(gint_t), 1, fp);
		xfwrite(&g->edges[e].rc_id, sizeof(gint_t), 1, fp);
		xfwrite(&g->edges[e].count, sizeof(uint64_t), 1, fp);
		xfwrite(&g->edges[e].seq_len, sizeof(gint_t), 1, fp);
		xfwrite(g->edges[e].seq, sizeof(uint32_t), (g->edges[e].seq_len + 15) >> 4, fp);
	}
	fclose(fp);
}

void load_asm_graph(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "rb");
	char sig[4];
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, "asmg", 4))
		__ERROR("Not assembly graph format file");
	xfread(&g->ksize, sizeof(int), 1, fp);
	xfread(&g->n_v, sizeof(gint_t), 1, fp);
	xfread(&g->n_e, sizeof(gint_t), 1, fp);
	g->nodes = malloc(g->n_v * sizeof(struct asm_node_t));
	g->edges = malloc(g->n_e * sizeof(struct asm_edge_t));
	gint_t u, e;
	for (u = 0; u < g->n_v; ++u) {
		xfread(&g->nodes[u].rc_id, sizeof(gint_t), 1, fp);
		xfread(&g->nodes[u].deg, sizeof(gint_t), 1, fp);
		if (g->nodes[u].deg) {
			g->nodes[u].adj = malloc(g->nodes[u].deg * sizeof(gint_t));
			xfread(g->nodes[u].adj, sizeof(gint_t), g->nodes[u].deg, fp);
		} else {
			g->nodes[u].adj = NULL;
		}
	}
	for (e = 0; e < g->n_e; ++e) {
		xfread(&g->edges[e].source, sizeof(gint_t), 1, fp);
		xfread(&g->edges[e].target, sizeof(gint_t), 1, fp);
		xfread(&g->edges[e].rc_id, sizeof(gint_t), 1, fp);
		xfread(&g->edges[e].count, sizeof(uint64_t), 1, fp);
		xfread(&g->edges[e].seq_len, sizeof(gint_t), 1, fp);
		g->edges[e].seq = calloc((g->edges[e].seq_len + 15) >> 4, sizeof(uint32_t));
		xfread(g->edges[e].seq, sizeof(uint32_t), (g->edges[e].seq_len + 15) >> 4, fp);
	}
	fclose(fp);
}

