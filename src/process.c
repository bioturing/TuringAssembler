#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "io_utils.h"
#include "k31hash.h"
#include "k63hash.h"
#include "kmer_count.h"
#include "process.h"
#include "resolve.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

void graph_convert_process(struct opt_build_t *opt)
{
	char path[1024];
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Dump graph from bin archive\n");
	struct asm_graph_t *g;
	g = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g, opt->in_path);
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g->ksize);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g->ksize);
	test_asm_graph(g);
	snprintf(path, 1024, "%s/graph_k_%d_loaded.gfa", opt->out_dir, g->ksize);
	write_gfa(g, path);
	snprintf(path, 1024, "%s/graph_k_%d_loaded.fasta", opt->out_dir, g->ksize);
	write_fasta(g, path);
}

void build_0_scratch(struct opt_count_t *opt, int ksize, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Building assembly graph from read using kmer size %d\n", ksize);
	if (ksize < 32)
		k31_build_scratch(opt, ksize, g);
	else if (ksize > 32 && ksize < 64)
		k63_build_scratch(opt, ksize, g);
	test_asm_graph(g);
}

void build_0_precount(struct opt_count_t *opt, int ksize_dst, int ksize_src, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Building assembly graph from read using kmer size %d\n", ksize_dst);
	if (ksize_dst < 32)
		k31_build_precount(opt, ksize_dst, ksize_src, g);
	else if (ksize_dst > 32 && ksize_dst < 64)
		k63_build_precount(opt, ksize_dst, ksize_src, g);
	test_asm_graph(g);
}

void build_0_1(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g0->ksize);
	set_time_now();
	remove_tips(g0, g);
	test_asm_graph(g);
	__VERBOSE_LOG("TIMER", "Build graph level 1 time: %.3f\n", sec_from_prev_time());
}

void build_1_2(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips using graph topology\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g0->ksize);
	set_time_now();
	remove_tips_topology(g0, g);
	test_asm_graph(g);
	__VERBOSE_LOG("TIMER", "Build graph level 2 time: %.3f\n", sec_from_prev_time());
}

void build_2_3(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips using graph topology\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g0->ksize);
	set_time_now();
	resolve_chain(g0, g);
	test_asm_graph(g);
	__VERBOSE_LOG("TIMER", "Build graph level 3 time: %.3f\n", sec_from_prev_time());
}

void build_3_4(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips using graph topology\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g0->ksize);
	set_time_now();
	resolve_n_m_simple(g0, g);
	test_asm_graph(g);
	__VERBOSE_LOG("TIMER", "Build graph level 3 time: %.3f\n", sec_from_prev_time());
}

void build_barcode(struct opt_build_t *opt, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Building barcode information\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g->ksize);
	set_time_now();
	construct_barcode_map(opt, g);
	__VERBOSE_LOG("TIMER", "Build barcode information time: %.3f\n", sec_from_prev_time());
}

void get_opt_build(struct opt_build_t *opt_b, struct opt_count_t *opt_c)
{
	opt_b->n_threads = opt_c->n_threads;
	opt_b->hash_size = opt_c->hash_size;
	opt_b->n_files = opt_c->n_files;
	opt_b->files_1 = opt_c->files_1;
	opt_b->files_2 = opt_c->files_2;
	opt_b->out_dir = opt_c->out_dir;
	opt_b->split_len = opt_c->split_len;
}

void build_barcode_opt_count(struct opt_count_t *opt, struct asm_graph_t *g)
{
	struct opt_build_t opt_b;
	get_opt_build(&opt_b, opt);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Building barcode information\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g->ksize);
	set_time_now();
	construct_barcode_map(&opt_b, g);
	__VERBOSE_LOG("TIMER", "Build barcode information time: %.3f\n", sec_from_prev_time());
}

void graph_query_process(struct opt_build_t *opt)
{
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

void save_graph_info(const char *out_dir, struct asm_graph_t *g, const char *suffix, int is_simple)
{
	char path[1024];
	__VERBOSE_LOG("graph_k_%d_%s", "Number of nodes: %ld\n",
						g->ksize, suffix, g->n_v);
	__VERBOSE_LOG("graph_k_%d_%s", "Number of edges: %ld\n",
						g->ksize, suffix, g->n_e);
	snprintf(path, 1024, "%s/graph_k_%d_%s.gfa",
						out_dir, g->ksize, suffix);
	write_gfa(g, path);
	snprintf(path, 1024, "%s/graph_k_%d_%s.fasta",
						out_dir, g->ksize, suffix);
	write_fasta(g, path);
	snprintf(path, 1024, "%s/graph_k_%d_%s.bin",
						out_dir, g->ksize, suffix);
	if (is_simple)
		save_asm_graph_simple(g, path);
	else
		save_asm_graph_barcode(g, path);
}

void assembly_process(struct opt_count_t *opt)
{
	struct asm_graph_t g1, g2;

	/* Count k0-mer and build graph */
	build_0_scratch(opt, opt->k0, &g1);
	save_graph_info(opt->out_dir, &g1, "level_0", 1);

	build_0_1(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_1", 1);
	asm_graph_destroy(&g1);

	build_1_2(&g2, &g1);
	save_graph_info(opt->out_dir, &g1, "level_2", 1);
	asm_graph_destroy(&g2);

	build_2_3(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_3", 1);
	asm_graph_destroy(&g1);

	build_barcode_opt_count(opt, &g2);
	build_3_4(&g2, &g1);
	save_graph_info(opt->out_dir, &g1, "level_4", 1);
	asm_graph_destroy(&g1);
	asm_graph_destroy(&g2);
}

void build_0_process(struct opt_count_t *opt)
{
	struct asm_graph_t g;

	build_0_scratch(opt, opt->k0, &g);
	save_graph_info(opt->out_dir, &g, "level_0", 1);
	asm_graph_destroy(&g);
}

void build_0_1_process(struct opt_build_t *opt)
{
	struct asm_graph_t g1, g2;
	load_asm_graph(&g1, opt->in_file);
	build_0_1(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_1", 1);
	asm_graph_destroy(&g1);
	asm_graph_destroy(&g2);
}

void build_1_2_process(struct opt_build_t *opt)
{
	struct asm_graph_t g1, g2;
	load_asm_graph(&g1, opt->in_file);
	build_1_2(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_2", 1);
	asm_graph_destroy(&g1);
	asm_graph_destroy(&g2);
}

void build_2_3_process(struct opt_build_t *opt)
{
	struct asm_graph_t g1, g2;
	load_asm_graph(&g1, opt->in_file);
	build_2_3(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_3", 1);
	asm_graph_destroy(&g1);
	asm_graph_destroy(&g2);
}

void build_3_4_process(struct opt_build_t *opt)
{
	struct asm_graph_t g1, g2;
	load_asm_graph(&g1, opt->in_file);
	build_barcode(opt, &g1);
	build_2_3(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_4", 1);
	asm_graph_destroy(&g1);
	asm_graph_destroy(&g2);
}

void build_3_4_no_bc_rebuild_process(struct opt_build_t *opt)
{
	struct asm_graph_t g1, g2;
	load_asm_graph(&g1, opt->in_file);
	build_2_3(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_4", 1);
	asm_graph_destroy(&g1);
	asm_graph_destroy(&g2);
}

void build_barcode_process(struct opt_build_t *opt)
{
	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	build_barcode(opt, &g);
	save_graph_info(opt->out_dir, &g, "added_barcode", 0);
	asm_graph_destroy(&g);
}
