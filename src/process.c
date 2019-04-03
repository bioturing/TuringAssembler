#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "io_utils.h"
#include "k31hash.h"
#include "k63hash.h"
#include "k31_count.h"
#include "k63_count.h"
#include "process.h"
#include "resolve.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

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
	write_fasta(g, path);
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
	write_fasta(g0, path);
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
	write_fasta(g1, path);
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
	write_fasta(g1, path);
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
	resolve_chain(g0, g1);
	// remove_bubble_and_loop(g0, g1);
	__VERBOSE_LOG("kmer_%d_graph_#3", "Number of nodes: %lld\n", g0->ksize,
							(long long)g1->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#3", "Number of edges: %lld\n", g0->ksize,
							(long long)g1->n_e);
	test_asm_graph(g1);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.gfa", opt->out_dir, g0->ksize);
	write_gfa(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.fasta", opt->out_dir, g0->ksize);
	write_fasta(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.bin", opt->out_dir, g0->ksize);
	save_asm_graph_simple(g1, path);
}

void build2_3a_process(struct opt_build_t *opt)
{
	char path[1024];
	init_clock();

	struct asm_graph_t *g0, *g1;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	g1 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_path);
	test_asm_graph(g0);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Checking complex regions\n");
	// asm_condense(g0, g1);
	// __VERBOSE("Number of nodes: %ld\n", g1->n_v);
	// __VERBOSE("Number of edges: %ld\n", g1->n_e);
	// snprintf(path, 1024, "%s/graph_k_%d_level_3a.gfa", opt->out_dir, g0->ksize);
	// write_gfa(g1, path);
	// snprintf(path, 1024, "%s/graph_k_%d_level_3a.bin", opt->out_dir, g0->ksize);
	// save_asm_graph_simple(g1, path);

	// detect_complex(g0, 2000, 2000);
	// construct_barcode_map(g0, opt);
	// __VERBOSE("\n");
	// print_test_barcode_edge(g0, 140382, 101945, opt->split_len);
	// print_test_barcode_edge(g0, 109035, 101945, opt->split_len);
}

void build3_4_process(struct opt_build_t *opt)
{
	char path[1024];
	init_clock();

	struct asm_graph_t *g0, *g1;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	g1 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_path);
	test_asm_graph(g0);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Collapsing 1-1 jungle\n");
	// resolve_bridge(g0);
	// collapse_1_1_jungle(g0, g1);
	collapse_n_m_jungle(g0, g1);
	__VERBOSE_LOG("kmer_%d_graph_#4", "Number of nodes: %lld\n", g0->ksize,
							(long long)g1->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#4", "Number of edges: %lld\n", g0->ksize,
							(long long)g1->n_e);
	test_asm_graph(g1);
	snprintf(path, 1024, "%s/graph_k_%d_level_4.gfa", opt->out_dir, g0->ksize);
	write_gfa(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_4.fasta", opt->out_dir, g0->ksize);
	write_fasta(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_4.bin", opt->out_dir, g0->ksize);
	save_asm_graph_simple(g1, path);
}

void build4_5_process(struct opt_build_t *opt)
{
	char path[1024];
	init_clock();

	struct asm_graph_t *g0, *g1;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	g1 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_path);
	test_asm_graph(g0);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Resolving 2-2, n-n bridge\n");
	// detect_simple_tandem(g0);
	collapse_n_m_jungle(g0, g1);
	__VERBOSE_LOG("kmer_%d_graph_#5", "Number of nodes: %lld\n", g0->ksize,
							(long long)g1->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#5", "Number of edges: %lld\n", g0->ksize,
							(long long)g1->n_e);
	test_asm_graph(g1);
	snprintf(path, 1024, "%s/graph_k_%d_level_5.gfa", opt->out_dir, g0->ksize);
	write_gfa(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_5.fasta", opt->out_dir, g0->ksize);
	write_fasta(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_5.bin", opt->out_dir, g0->ksize);
	save_asm_graph_simple(g1, path);
}

void build5_6_process(struct opt_build_t *opt)
{
	char path[1024];
	init_clock();

	struct asm_graph_t *g0, *g1;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	g1 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_path);
	test_asm_graph(g0);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Expanding edges\n");
	// graph_expanding(g0);
	resolve_chain(g0, g1);
	// detect_simple_tandem(g0);
	// struct asm_graph_t *g1, *g2;
	// g1 = calloc(1, sizeof(struct asm_graph_t));
	// g2 = calloc(1, sizeof(struct asm_graph_t));
	// asm_condense(g0, g1);
	// resolve_loop_bugle(g1, &g1);
	// asm_condense(g0, g1);
	__VERBOSE_LOG("kmer_%d_graph_#6", "Number of nodes: %lld\n", g0->ksize,
							(long long)g1->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#6", "Number of edges: %lld\n", g0->ksize,
							(long long)g1->n_e);
	test_asm_graph(g1);
	snprintf(path, 1024, "%s/graph_k_%d_level_6.gfa", opt->out_dir, g0->ksize);
	write_gfa(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_6.fasta", opt->out_dir, g0->ksize);
	write_fasta(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_6.bin", opt->out_dir, g0->ksize);
	save_asm_graph_simple(g1, path);
}

void build6_7_process(struct opt_build_t *opt)
{
	char path[1024];
	init_clock();

	struct asm_graph_t *g0, *g1;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	g1 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_path);
	test_asm_graph(g0);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Resolving bridge complex\n");
	// graph_expanding(g0);
	// resolve_chain(g0, g1);
	// resolve_bridge_complex(g0, g1);
	// detect_simple_tandem(g0);
	// struct asm_graph_t *g1, *g2;
	// g1 = calloc(1, sizeof(struct asm_graph_t));
	// g2 = calloc(1, sizeof(struct asm_graph_t));
	// asm_condense(g0, g1);
	// resolve_loop_bugle(g1, &g1);
	// asm_condense(g0, g1);

	// __VERBOSE_LOG("kmer_%d_graph_#6", "Number of nodes: %lld\n", g0->ksize,
	// 						(long long)g1->n_v);
	// __VERBOSE_LOG("kmer_%d_graph_#6", "Number of edges: %lld\n", g0->ksize,
	// 						(long long)g1->n_e);
	// test_asm_graph(g1);
	// snprintf(path, 1024, "%s/graph_k_%d_level_6.gfa", opt->out_dir, g0->ksize);
	// write_gfa(g1, path);
	// snprintf(path, 1024, "%s/graph_k_%d_level_6.fasta", opt->out_dir, g0->ksize);
	// write_fasta(g1, path);
	// snprintf(path, 1024, "%s/graph_k_%d_level_6.bin", opt->out_dir, g0->ksize);
	// save_asm_graph_simple(g1, path);
}

void clean_process(struct opt_build_t *opt)
{
	char path[1024];
	init_clock();

	struct asm_graph_t *g0,  *g1;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	g1 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_path);
	test_asm_graph(g0);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Testing clean\n");
	// resolve_loop_bugle(g0, &g1);
	test_asm_graph(g1);
	snprintf(path, 1024, "%s/graph_k_%d_level_5.gfa", opt->out_dir, g1->ksize);
	write_gfa(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_5.bin", opt->out_dir, g0->ksize);
	save_asm_graph_simple(g1, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_5.fasta", opt->out_dir, g0->ksize);
	write_fasta(g1, path);
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
		int qret = test_edge_barcode(g0, u, v);
		fprintf(stdout, "ret = %d\n", qret);
		print_test_barcode_edge(g0, u, v);
	}
	fclose(fp);
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
	write_fasta(g0, path);
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
	write_fasta(g1, path);
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
	write_fasta(g2, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_2.bin", opt->out_dir, opt->k0);
	save_asm_graph_simple(g2, path);

	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing bubbles\n");
	struct asm_graph_t *g3;
	g3 = calloc(1, sizeof(struct asm_graph_t));
	resolve_chain(g2, g3);
	// remove_bubble_and_loop(g2, g3);
	__VERBOSE_LOG("kmer_%d_graph_#3", "Number of nodes: %lld\n", opt->k0,
							(long long)g3->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#3", "Number of edges: %lld\n", opt->k0,
							(long long)g3->n_e);
	test_asm_graph(g3);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.gfa", opt->out_dir, opt->k0);
	write_gfa(g3, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.fasta", opt->out_dir, opt->k0);
	write_fasta(g3, path);
	snprintf(path, 1024, "%s/graph_k_%d_level_3.bin", opt->out_dir, opt->k0);
	save_asm_graph_simple(g3, path);
}

