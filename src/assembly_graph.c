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

#define __bin_seq_get_char(seq, l) (((seq)[(l) >> 4] >> (((l) & 15) << 1)) & (uint32_t)0x3)

static void dump_edge_seq(char **seq, uint32_t *m_seq, struct asm_edge_t *e)
{
	uint32_t i, j, k, len = e->seq_len;
	for (i = 0; i < e->n_holes; ++i)
		len += e->l_holes[i];
	if (*m_seq < len + 1) {
		*m_seq = len + 1;
		*seq = realloc(*seq, *m_seq);
	}
	j = k = 0;
	for (i = 0; i < e->seq_len; ++i) {
		(*seq)[k++] = nt4_char[__bin_seq_get_char(e->seq, i)];
		/* append holes to the sequences */
		if (j < e->n_holes && e->p_holes[j] == i) {
			/* fill with 'N's */
			memset((*seq) + k, 0x4e, e->l_holes[j]);
			k += e->l_holes[j];
			++j;
		}
	}
	(*seq)[k++] = '\0';
}

static void dump_bin_seq(char *seq, uint32_t *bin, gint_t len)
{
	gint_t i;
	for (i = 0; i < len; ++i)
		seq[i] = nt4_char[(bin[i >> 4] >> ((i & 15) << 1)) & 3];
	seq[len] = '\0';
}

#define TIPS_THRESHOLD			5.0
#define TIPS_RATIO_THRESHOLD		0.1

#define NON_TIPS_LEN			250

void remove_tips(struct asm_graph_t *g0, struct asm_graph_t *g);
void remove_tips_topology(struct asm_graph_t *g0, struct asm_graph_t *g);
void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void write_gfa(struct asm_graph_t *g, const char *path);
void dump_fasta(struct asm_graph_t *g, const char *path);
void remove_bubble_and_loop(struct asm_graph_t *g0, struct asm_graph_t *g);
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
	save_asm_graph_simple(g0, path);
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
	save_asm_graph_simple(g1, path);
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
	save_asm_graph_simple(g1, path);
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
	remove_bubble_and_loop(g0, g1);
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
	save_asm_graph_simple(g1, path);
}

void build_barcode_process(struct opt_build_t *opt)
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
	__VERBOSE("\n");

	snprintf(path, 1024, "%s/graph_k_%d_added_barcode.bin", opt->out_dir, g0->ksize);
	save_asm_graph_barcode(g0, path);
}

void graph_query_process(struct opt_build_t *opt)
{
	init_clock();

	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_path);
	fprintf(stderr, "bin size = %d\n", g0->bin_size);
	test_asm_graph(g0);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Querying pair-edge barcode information\n");
	gint_t u, v;
	FILE *fp;
	fp = xfopen(opt->in_file, "r");
	while (1) {
		int ret = fscanf(fp, "%ld%ld", &u, &v);
		if (ret == EOF || ret == 0)
			break;
		print_test_barcode_edge(g0, u, v);
	}
	fclose(fp);

}

void build2_3a_process(struct opt_build_t *opt)
{
	// char path[1024];
	// init_clock();

	// struct asm_graph_t *g0;
	// g0 = calloc(1, sizeof(struct asm_graph_t));
	// load_asm_graph(g0, opt->in_path);
	// test_asm_graph(g0);
	// __VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	// __VERBOSE("\n+------------------------------------------------------------------------------+\n");
	// __VERBOSE("Building barcode information\n");
	// construct_barcode_map(g0, opt);
	// __VERBOSE("\n");
	// print_test_barcode_edge(g0, 140382, 101945, opt->split_len);
	// print_test_barcode_edge(g0, 109035, 101945, opt->split_len);
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
	save_asm_graph_simple(g0, path);

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
	save_asm_graph_simple(g1, path);

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
	save_asm_graph_simple(g2, path);

	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing bubbles\n");
	struct asm_graph_t *g3;
	g3 = calloc(1, sizeof(struct asm_graph_t));
	remove_bubble_and_loop(g2, g3);
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
	save_asm_graph_simple(g3, path);
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
	save_asm_graph_simple(g0, path);

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
	save_asm_graph_simple(g1, path);


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
	save_asm_graph_simple(g2, path);

	__VERBOSE("\nRemoving bubble\n");
	struct asm_graph_t *g3;
	g3 = calloc(1, sizeof(struct asm_graph_t));
	remove_bubble_and_loop(g2, g3);
	__VERBOSE_LOG("Graph #4", "Number of nodes: %lld\n", (long long)g3->n_v);
	__VERBOSE_LOG("Graph #4", "Number of edges: %lld\n", (long long)g3->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_3.gfa");
	write_gfa(g3, path);
	test_asm_graph(g3);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_3.fasta");
	dump_fasta(g3, path);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_3.bin");
	save_asm_graph_simple(g3, path);
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
	save_asm_graph_simple(g0, path);

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
	save_asm_graph_simple(g1, path);

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
	save_asm_graph_simple(g2, path);

	__VERBOSE("\nRemoving bubble\n");
	struct asm_graph_t *g3;
	g3 = calloc(1, sizeof(struct asm_graph_t));
	remove_bubble_and_loop(g2, g3);
	__VERBOSE_LOG("Graph #4", "Number of nodes: %lld\n", (long long)g3->n_v);
	__VERBOSE_LOG("Graph #4", "Number of edges: %lld\n", (long long)g3->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_3.gfa");
	write_gfa(g3, path);
	test_asm_graph(g3);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_3.fasta");
	dump_fasta(g3, path);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_3.bin");
	save_asm_graph_simple(g3, path);
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

static gint_t asm_find_id(gint_t *adj, gint_t deg, gint_t id)
{
	gint_t i, ret;
	ret = -1;
	for (i = 0; i < deg; ++i) {
		if (adj[i] == id)
			ret = i;
	}
	return ret;
}

int asm_is_edge_rc(uint32_t *seq1, uint32_t l1, uint32_t *seq2, uint32_t l2)
{
	if (l1 != l2)
		return 0;
	uint32_t c1, c2, i, k;
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
	uint32_t seq_len = 0;
	char *buf = alloca(81);
	gint_t e, e_rc;
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		gint_t cc_id = id_edge[e];
		if (cc_size[cc_id] < 250 || g->edges[e].seq_len < 100)
			continue;
		// if (seq_len < g->edges[e].seq_len + 1) {
		// 	seq_len = g->edges[e].seq_len + 1;
		// 	seq = realloc(seq, seq_len);
		// }
		// dump_bin_seq(seq, g->edges[e].seq, g->edges[e].seq_len);
		dump_edge_seq(&seq, &seq_len, g->edges + e);
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

	free(seq);
	free(id_edge);
	free(cc_size);
}

static inline uint64_t get_bandage_count(struct asm_edge_t *e, int ksize)
{
	double cov = e->count * 1.0 / ((e->n_holes + 1) * ksize);
	uint32_t i, len = e->seq_len;
	for (i = 0; i < e->n_holes; ++i)
		len += e->l_holes[i];
	return (uint64_t)(cov * len);
}

void write_gfa(struct asm_graph_t *g, const char *path)
{
	gint_t *id_edge, *cc_size;
	id_edge = malloc(g->n_e * sizeof(gint_t));
	cc_size = NULL;
	asm_edge_cc(g, id_edge, &cc_size);

	FILE *fp = xfopen(path, "w");
	char *seq = NULL;
	uint32_t seq_len = 0;
	gint_t e, e_rc;
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		gint_t cc_id = id_edge[e];
		if (cc_size[cc_id] < 250)
			continue;
		// if (seq_len < g->edges[e].seq_len + 1) {
		// 	seq_len = g->edges[e].seq_len + 1;
		// 	seq = realloc(seq, seq_len);
		// }
		// dump_bin_seq(seq, g->edges[e].seq, g->edges[e].seq_len);
		dump_edge_seq(&seq, &seq_len, g->edges + e);
		uint64_t fake_count = get_bandage_count(g->edges + e, g->ksize);
		// double cov = g->edges[e].count * 1.0 / (g->edges[e].seq_len - g->ksize);
		// uint64_t fake_count = (uint64_t)(cov * g->edges[e].seq_len);
		/* print fake count for correct coverage display on Bandage */
		fprintf(fp, "S\t%lld_%lld\t%s\tKC:i:%llu\n", (long long)e,
			(long long)e_rc, seq, (long long unsigned)fake_count);
	}
	for (e = 0; e < g->n_e; ++e) {
		gint_t cc_id = id_edge[e];
		if (cc_size[cc_id] < 250)
			continue;
		e_rc = g->edges[e].rc_id;
		gint_t pe, pe_rc, next_pe, next_pe_rc;
		char ce, next_ce;
		if (e > e_rc) {
			pe = e_rc;
			pe_rc = e;
			ce = '-';
		} else {
			pe = e;
			pe_rc = e_rc;
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
				next_pe_rc = next_e;
				next_ce = '-';
			} else {
				next_pe = next_e;
				next_pe_rc = next_e_rc;
				next_ce = '+';
			}
			fprintf(fp, "L\t%lld_%lld\t%c\t%lld_%lld\t%c\t%dM\n",
				(long long)pe, (long long)pe_rc, ce,
				(long long)next_pe, (long long)next_pe_rc, next_ce,
				g->ksize);
		}
	}
	fclose(fp);
	free(seq);
	free(id_edge);
	free(cc_size);
}

void remove_self_loop(struct asm_graph_t *g, double uni_cov)
{
	/* Loop is bad
	 *               +---+
	 *              /    |
	 *             +     +  <---------- loop
	 *             |    /
	 * o==========>o<--+
	 *             \\
	 *               ++========>o
	 */
	gint_t e, e_rc, u, u_rc;
	gint_t cnt = 0;
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		if (e != e_rc)
			continue;
		u = g->edges[e].source;
		u_rc = g->nodes[e].rc_id;
		assert(g->edges[e].target == u_rc);
		if (g->nodes[u_rc].deg == 1 && g->nodes[u].deg == 2) {
			gint_t e_true, e_prev;
			if (g->nodes[u].adj[0] == e)
				e_true = g->nodes[u].adj[1];
			else
				e_true = g->nodes[u].adj[0];
			e_prev = g->nodes[u_rc].adj[0];
			double cov1, cov2;
			cov1 = g->edges[e_prev].count / (g->edges[e_prev].seq_len - g->ksize);
			cov2 = g->edges[e_true].count / (g->edges[e_true].seq_len - g->ksize);
			if (cov1 / cov2 < 0.5 || cov1 / cov2 > 1.5)
				continue;
			g->nodes[u].adj[0] = e_true;
			g->nodes[u].deg = 1;
			g->nodes[u].adj = realloc(g->nodes[u].adj, sizeof(gint_t));
			++cnt;
		}
	}
	__VERBOSE("Number of removed loops: %ld\n", cnt);
}

void remove_bubble_simple(struct asm_graph_t *g0, double uni_cov)
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
	gint_t cnt = 0;
	gint_t u;
	for (u = 0; u < g0->n_v; ++u) {
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
		if (cov / uni_cov < 0.51 || cov / uni_cov > 1.51)
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
		++cnt;
	}
	__VERBOSE("\nNumber of collapsed bubble: %ld\n", cnt);
}

void remove_bubble_and_loop(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	double uni_cov = get_uni_genome_coverage(g0);
	__VERBOSE("1 genome walk coverage: %lf\n", uni_cov);
	remove_self_loop(g0, uni_cov);
	remove_bubble_simple(g0, uni_cov);
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
		gint_t c = asm_find_id(g0->nodes[u].adj, g0->nodes[u].deg, e);
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
				gint_t j = asm_find_id(g0->nodes[v_rc].adj,
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

void asm_clone_edge(struct asm_edge_t *dst, struct asm_edge_t *src)
{
	dst->count = src->count;
	dst->seq_len = src->seq_len;
	dst->seq = calloc((dst->seq_len + 15) >> 4, sizeof(uint32_t));
	memcpy(dst->seq, src->seq,
		((dst->seq_len + 15) >> 4) * sizeof(uint32_t));
	dst->n_holes = src->n_holes;
	dst->p_holes = malloc(dst->n_holes * sizeof(uint32_t));
	memcpy(dst->p_holes, src->p_holes, dst->n_holes * sizeof(uint32_t));
	dst->l_holes = malloc(dst->n_holes * sizeof(uint32_t));
	memcpy(dst->l_holes, src->l_holes, dst->n_holes * sizeof(uint32_t));
}

void asm_clone_reverse(struct asm_edge_t *dst, struct asm_edge_t *src)
{
	dst->count = src->count;
	dst->seq_len = src->seq_len;
	dst->seq = calloc((dst->seq_len + 15) >> 4, sizeof(uint32_t));
	uint32_t i, k;
	for (i = 0; i < dst->seq_len; ++i) {
		k = dst->seq_len - i - 1;
		dst->seq[i >> 4] |= (uint32_t)(__bin_seq_get_char(src->seq, k) ^ 3)
							<< ((i & 15) << 1);
	}
	dst->n_holes = src->n_holes;
	dst->l_holes = malloc(dst->n_holes * sizeof(uint32_t));
	dst->p_holes = malloc(dst->n_holes * sizeof(uint32_t));
	for (i = 0; i < dst->n_holes; ++i) {
		dst->l_holes[i] = src->l_holes[dst->n_holes - i - 1];
		dst->p_holes[i] = dst->seq_len - 1
			- (dst->p_holes[dst->n_holes - i - 1] + 1);
	}
}

void asm_append_edge(struct asm_edge_t *dst, struct asm_edge_t *src, uint32_t overlap)
{
	uint32_t seq_len, new_m, m;
	dst->count += src->count;
	/* append the bin seq */
	seq_len = dst->seq_len + src->seq_len - overlap;
	new_m = (seq_len + 15) >> 4;
	m = (dst->seq_len + 15) >> 4;
	if (new_m > m) {
		dst->seq = realloc(dst->seq, new_m * sizeof(uint32_t));
		memset(dst->seq + m, 0, (new_m - m) * sizeof(uint32_t));
	}

	uint32_t i, k;
	for (i = overlap; i < src->seq_len; ++i) {
		k = i - overlap + dst->seq_len;
		dst->seq[k >> 4] |= ((src->seq[i >> 4] >> ((i & 15) << 1)) & 3)
							<< ((k & 15) << 1);
	}
	/* append the gaps */
	if (src->n_holes) {
		uint32_t n_holes = dst->n_holes + src->n_holes;
		dst->p_holes = realloc(dst->p_holes, n_holes * sizeof(uint32_t));
		dst->l_holes = realloc(dst->l_holes, n_holes * sizeof(uint32_t));
		for (i = 0; i < src->n_holes; ++i)
			dst->p_holes[dst->n_holes + i] = src->p_holes[i] + dst->seq_len;
		memcpy(dst->l_holes + dst->n_holes, src->l_holes, src->n_holes * sizeof(uint32_t));
		dst->n_holes = n_holes;
	}
	dst->seq_len = seq_len;
}

void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	gint_t *node_id;
	gint_t n_v, n_e, m_e;
	node_id = malloc(g0->n_v * sizeof(gint_t));
	memset(node_id, 255, g0->n_v * sizeof(gint_t));

	/* nodes on new graph only consist of branching nodes on old graph */
	n_v = 0;
	gint_t u;
	for (u = 0; u < g0->n_v; ++u) {
		gint_t deg_fw = g0->nodes[u].deg;
		gint_t deg_rv = g0->nodes[g0->nodes[u].rc_id].deg;
		/* non-branching node */
		if ((deg_fw == 1 && deg_rv == 1) || deg_fw + deg_rv == 0)
			continue;
		node_id[u] = n_v++;
	}
	struct asm_node_t *nodes = calloc(n_v, sizeof(struct asm_node_t));
	/* set reverse complement link between nodes */
	for (u = 0; u < g0->n_v; ++u) {
		gint_t x, x_rc, u_rc;
		x = node_id[u];
		if (x == -1)
			continue;
		u_rc = g0->nodes[u].rc_id;
		x_rc = node_id[u_rc];
		assert(x_rc != -1);
		nodes[x].rc_id = x_rc;
		nodes[x_rc].rc_id = x;
		nodes[x].adj = NULL;
		nodes[x].deg = 0;
	}
	n_e = 0;
	m_e = 0x10000;
	struct asm_edge_t *edges = calloc(m_e, sizeof(struct asm_edge_t));
	/* construct new edges */
	for (u = 0; u < g0->n_v; ++u) {
		gint_t x, y_rc;
		x = node_id[u];
		if (x == -1)
			continue;
		gint_t c;
		for (c = 0; c < g0->nodes[u].deg; ++c) {
			if (n_e + 2 > m_e) {
				edges = realloc(edges, (m_e << 1) * sizeof(struct asm_edge_t));
				memset(edges + m_e, 0, m_e * sizeof(struct asm_edge_t));
				m_e <<= 1;
			}
			gint_t e = g0->nodes[u].adj[c], e_rc, v, v_rc, p, q;
			if (e == -1)
				continue;
			p = n_e; q = n_e + 1;
			edges[p].rc_id = q;
			edges[q].rc_id = p;
			asm_clone_edge(edges + p, g0->edges + e);

			do {
				v = g0->edges[e].target;
				if (node_id[v] == -1) { /* middle node */
					assert(g0->nodes[v].deg == 1);
					e = g0->nodes[v].adj[0];
					assert(e != -1);
					asm_append_edge(edges + p, g0->edges + e, g0->ksize);
				} else {
					break;
				}
			} while (1);
			edges[p].source = x;
			edges[p].target = node_id[v];
			asm_clone_reverse(edges + q, edges + p);
			v_rc = g0->nodes[v].rc_id;
			e_rc = g0->edges[e].rc_id;
			gint_t j = asm_find_id(g0->nodes[v_rc].adj,
						g0->nodes[v_rc].deg, e_rc);
			assert(j >= 0);
			g0->nodes[v_rc].adj[j] = -1;
			y_rc = node_id[v_rc];
			edges[q].source = y_rc;
			edges[q].target = nodes[x].rc_id;

			nodes[x].adj = realloc(nodes[x].adj, (nodes[x].deg + 1) * sizeof(gint_t));
			nodes[x].adj[nodes[x].deg++] = p;
			nodes[y_rc].adj = realloc(nodes[y_rc].adj,
					(nodes[y_rc].deg + 1) * sizeof(gint_t));
			nodes[y_rc].adj[nodes[y_rc].deg++] = q;
			n_e += 2;
		}
	}
	free(node_id);
	edges = realloc(edges, n_e * sizeof(struct asm_edge_t));
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
				fprintf(stderr, "seq_len = %u; rc_seq_len = %u\n",
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
					for (k = 0; k < g->nodes[u].deg; ++k) {
						gint_t e = g->nodes[u].adj[k];
						char *seq = malloc(g->edges[e].seq_len);
						dump_bin_seq(seq, g->edges[e].seq,
							g->edges[e].seq_len);
						fprintf(stderr, "e = %ld; seq = %s\n",
							e, seq);
					}
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
			fprintf(stderr, "seq_len = %u; rc_seq_len = %u\n",
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

void save_graph(struct asm_graph_t *g, FILE *fp)
{
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
		xfwrite(&g->edges[e].n_holes, sizeof(uint32_t), 1, fp);
		if (g->edges[e].n_holes) {
			xfwrite(g->edges[e].p_holes, sizeof(uint32_t), g->edges[e].n_holes, fp);
			xfwrite(g->edges[e].l_holes, sizeof(uint32_t), g->edges[e].n_holes, fp);
		}
	}
}

void save_asm_graph_simple(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "wb");
	char *sig = "asmg";
	xfwrite(sig, 4, 1, fp);
	save_graph(g, fp);
	fclose(fp);
}

void save_asm_graph_barcode(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "wb");
	char *sig = "asmb";
	xfwrite(sig, 4, 1, fp);
	save_graph(g, fp);
	xfwrite(&g->bin_size, sizeof(int), 1, fp);
	gint_t e;
	for (e = 0; e < g->n_e; ++e) {
		gint_t e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		gint_t n, k;
		n = (g->edges[e].seq_len + g->bin_size - 1) / g->bin_size;
		for (k = 0; k < n; ++k) {
			struct barcode_hash_t *h = g->edges[e].bucks + k;
			xfwrite(&h->size, sizeof(uint32_t), 1, fp);
			xfwrite(&h->n_item, sizeof(uint32_t), 1, fp);
			xfwrite(h->keys, sizeof(uint64_t), h->size, fp);
			xfwrite(h->cnts, sizeof(uint32_t), h->size, fp);
		}
	}
	fclose(fp);
}

void load_graph(struct asm_graph_t *g, FILE *fp)
{
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
		xfread(&g->edges[e].n_holes, sizeof(uint32_t), 1, fp);
		if (g->edges[e].n_holes) {
			g->edges[e].p_holes = malloc(g->edges[e].n_holes * sizeof(uint32_t));
			g->edges[e].l_holes = malloc(g->edges[e].n_holes * sizeof(uint32_t));
			xfread(g->edges[e].p_holes, sizeof(uint32_t), g->edges[e].n_holes, fp);
			xfread(g->edges[e].l_holes, sizeof(uint32_t), g->edges[e].n_holes, fp);
		} else {
			g->edges[e].p_holes = NULL;
			g->edges[e].l_holes = NULL;
		}
	}
}

void load_barcode(struct asm_graph_t *g, FILE *fp)
{
	xfread(&g->bin_size, sizeof(int), 1, fp);
	gint_t e;
	for (e = 0; e < g->n_e; ++e) {
		gint_t e_rc = g->edges[e].rc_id;
		if (e > e_rc) {
			g->edges[e].bucks = g->edges[e_rc].bucks;
			continue;
		}
		gint_t n, k;
		n = (g->edges[e].seq_len + g->bin_size - 1) / g->bin_size;
		g->edges[e].bucks = calloc(n, sizeof(struct barcode_hash_t));
		for (k = 0; k < n; ++k) {
			struct barcode_hash_t *h = g->edges[e].bucks + k;
			xfread(&h->size, sizeof(uint32_t), 1, fp);
			xfread(&h->n_item, sizeof(uint32_t), 1, fp);
			h->keys = malloc(h->size * sizeof(uint64_t));
			h->cnts = malloc(h->size * sizeof(uint32_t));
			xfread(h->keys, sizeof(uint64_t), h->size, fp);
			xfread(h->cnts, sizeof(uint32_t), h->size, fp);
		}
	}
}

void load_asm_graph(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "rb");
	char sig[4];
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, "asmg", 4) && strncmp(sig, "asmb", 4))
		__ERROR("Not assembly graph format file");
	load_graph(g, fp);
	if (strncmp(sig, "asmb", 4) == 0)
		load_barcode(g, fp);
	fclose(fp);
}

void load_asm_graph_simple(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "rb");
	char sig[4];
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, "asmg", 4) && strncmp(sig, "asmb", 4))
		__ERROR("Not assembly graph format file");
	load_graph(g, fp);
	fclose(fp);
}

void load_asm_graph_complex(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "rb");
	char sig[4];
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, "asmb", 4))
		__ERROR("Not assembly graph with barcode format file");
	load_graph(g, fp);
	load_barcode(g, fp);
	fclose(fp);
}
