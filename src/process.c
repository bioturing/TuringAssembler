#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "io_utils.h"
#include "process.h"
#include "resolve.h"
#include "sort_read.h"
#include "time_utils.h"
#include "utils.h"
#include "verbose.h"

void graph_convert_process(struct opt_proc_t *opt)
{
	char path[1024];
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Dump graph from bin archive\n");
	struct asm_graph_t *g;
	g = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g, opt->in_file);
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g->ksize);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g->ksize);
	test_asm_graph(g);
	snprintf(path, 1024, "%s/graph_k_%d_loaded.gfa", opt->out_dir, g->ksize);
	write_gfa(g, path);
	snprintf(path, 1024, "%s/graph_k_%d_loaded.fasta", opt->out_dir, g->ksize);
	write_fasta(g, path);
}

void build_0_KMC(struct opt_proc_t *opt, int ksize, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Building assembly graph from read using kmer size %d\n", ksize);
	build_initial_graph(opt, ksize, g);
	// graph_build_KMC(opt, ksize, g);
	test_asm_graph(g);
}

void build_0_1(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Resolve graph using small operation\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g0->ksize);
	set_time_now();
	resolve_graph_operation(g0, g);
	// remove_tips(g0, g);
	test_asm_graph(g);
	__VERBOSE_LOG("TIMER", "Build graph level 1 time: %.3f\n", sec_from_prev_time());
}

void build_1_2b(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Resolving graph using barcode (local assembly included)\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g0->ksize);
	set_time_now();
	// resolve_tips_topo(g0, g);
	// remove_tips_topology(g0, g);
	// test_asm_graph(g);
	__VERBOSE_LOG("TIMER", "Build graph level 2 time: %.3f\n", sec_from_prev_time());

}

void build_1_2(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips using graph topology\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g0->ksize);
	set_time_now();
	// resolve_tips_topo(g0, g);
	// remove_tips_topology(g0, g);
	test_asm_graph(g);
	__VERBOSE_LOG("TIMER", "Build graph level 2 time: %.3f\n", sec_from_prev_time());
}

void build_2_3(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Resolving small complex structure\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g0->ksize);
	set_time_now();
	// resolve_chain(g0, g);
	test_asm_graph(g);
	__VERBOSE_LOG("TIMER", "Build graph level 3 time: %.3f\n", sec_from_prev_time());
}

void build_3_4(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Resolving graph using read pair + barcode (simple repetitive)\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g0->ksize);
	set_time_now();
	resolve_n_m_simple(g0, g);
	test_asm_graph(g);
	__VERBOSE_LOG("TIMER", "Build graph level 4 time: %.3f\n", sec_from_prev_time());
}

void build_4_5(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Resolving graph using read pair + barcode (complex jungle)\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g0->ksize);
	set_time_now();
	resolve_complex(g0, g);
	test_asm_graph(g);
	__VERBOSE_LOG("TIMER", "Build graph level 5 time: %.3f\n", sec_from_prev_time());
}

void build_barcode(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Building barcode information\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g->ksize);
	set_time_now();
	// construct_aux_information(opt, g, ASM_BUILD_BARCODE | ASM_BUILD_READPAIR);
	__VERBOSE_LOG("TIMER", "Build barcode information time: %.3f\n", sec_from_prev_time());
}

void build_barcode_read(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Building barcode information\n");
	__VERBOSE_LOG("INFO", "Input graph kmer size: %d\n", g->ksize);
	set_time_now();
	// construct_aux_information(opt, g, ASM_BUILD_BARCODE | ASM_BUILD_READPAIR | ASM_BUILD_COVERAGE);
	__VERBOSE_LOG("TIMER", "Build barcode information time: %.3f\n", sec_from_prev_time());
}

void graph_query_process(struct opt_proc_t *opt)
{
	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_file);
	fprintf(stderr, "bin size = %d\n", g0->bin_size);
	test_asm_graph(g0);
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Querying pair-edge barcode information\n");
	gint_t u, v, v2;
	FILE *fp;
	fp = xfopen(opt->in_fasta, "r");
	while (1) {
		char c;
		int ret = fscanf(fp, "%c %ld", &c, &u);
		if (ret == EOF || ret == 0)
			break;
		if (c == 'L') {
			fscanf(fp, "%ld\n", &v);
			print_test_barcode_edge(g0, u, v);
		} else if (c == 'P') {
			fscanf(fp, "%ld\n", &v);
			test_local_assembly(opt, g0, u, v);
		} else if (c == 'S') {
			fscanf(fp, "%ld %ld\n", &v, &v2);
			print_test_barcode_superior(g0, u, v, v2);
			// if (check_medium_pair_superior(g0, u, v, v2)) {
			// 	printf("success\n");
			// } else {
			// 	printf("failed\n");
			// }
		}
		// int qret = test_edge_barcode(g0, u, v);
		// fprintf(stdout, "ret = %d\n", qret);
	}
	fclose(fp);
}

void save_graph_info(const char *out_dir, struct asm_graph_t *g, const char *suffix)
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
	save_asm_graph(g, path);
}

void build_barcode_info(struct opt_proc_t *opt)
{
	struct asm_graph_t g;
	struct read_path_t read_sorted_path;

	load_asm_graph(&g, opt->in_file);
	char fasta_path[MAX_PATH];
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
		sort_read(opt, &read_sorted_path);
	}
	sprintf(fasta_path, "%s/barcode_build_dir", opt->out_dir);
	mkdir(fasta_path, 0755);
	sprintf(fasta_path, "%s/barcode_build_dir/contigs_tmp.fasta", opt->out_dir);
	write_fasta_seq(&g, fasta_path);
	construct_aux_info(opt, &g, &read_sorted_path, fasta_path, ASM_BUILD_BARCODE);
	save_graph_info(opt->out_dir, &g, "added_barcode");
	asm_graph_destroy(&g);
}

void assembly_process(struct opt_proc_t *opt)
{
	struct asm_graph_t g1, g2;
	struct read_path_t read_sorted_path;

	load_asm_graph(&g1, opt->in_file);
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
		sort_read(opt, &read_sorted_path);
	}
	resolve_n_m_local(opt, &read_sorted_path, &g1, &g2);
	// save_graph_info(opt->out_dir, &g2, "level_2");
}

void assembly3_process(struct opt_proc_t *opt)
{
	struct asm_graph_t g1, g2;
	build_0_KMC(opt, opt->k0, &g1);
	save_graph_info(opt->out_dir, &g1, "level_0");

	build_0_1(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_1");
	asm_graph_destroy(&g1);

	// build_barcode(opt, &g2);
	// build_3_4(&g2, &g1);
	// save_graph_info(opt->out_dir, &g1, "level_4");
	// asm_graph_destroy(&g2);

	// build_barcode(opt, &g1);
	// build_4_5(&g1, &g2);
	// save_graph_info(opt->out_dir, &g2, "level_5");

	// asm_graph_destroy(&g1);
	// asm_graph_destroy(&g2);
}

void build_0_process(struct opt_proc_t *opt)
{
	struct asm_graph_t g;
	build_0_KMC(opt, opt->k0, &g);
	save_graph_info(opt->out_dir, &g, "level_0");
	asm_graph_destroy(&g);
}

void build_0_1_process(struct opt_proc_t *opt)
{
	struct asm_graph_t g1, g2;
	load_asm_graph(&g1, opt->in_file);
	build_0_1(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_1");
	// asm_graph_destroy(&g1);
	// asm_graph_destroy(&g2);
}

void build_1_2_process(struct opt_proc_t *opt)
{
	struct asm_graph_t g1, g2;
	load_asm_graph(&g1, opt->in_file);
	do_something_local(opt, &g1);
	// build_1_2(&g1, &g2);
	// save_graph_info(opt->out_dir, &g2, "level_2");
	// asm_graph_destroy(&g1);
	// asm_graph_destroy(&g2);
}

void build_2_3_process(struct opt_proc_t *opt)
{
	struct asm_graph_t g1, g2;
	load_asm_graph(&g1, opt->in_file);
	build_2_3(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_3");
	asm_graph_destroy(&g1);
	asm_graph_destroy(&g2);
}

void build_3_4_process(struct opt_proc_t *opt)
{
	struct asm_graph_t g1, g2;
	load_asm_graph(&g1, opt->in_file);
	build_barcode(opt, &g1);
	build_3_4(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_4");
	asm_graph_destroy(&g1);
	asm_graph_destroy(&g2);
}

void build_3_4_no_bc_rebuild_process(struct opt_proc_t *opt)
{
	struct asm_graph_t g1, g2;
	load_asm_graph(&g1, opt->in_file);
	build_3_4(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_4");
	asm_graph_destroy(&g1);
	asm_graph_destroy(&g2);
}

void build_4_5_process(struct opt_proc_t *opt)
{
	struct asm_graph_t g1, g2;
	load_asm_graph(&g1, opt->in_file);
	build_barcode(opt, &g1);
	build_4_5(&g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_5");
	asm_graph_destroy(&g1);
	asm_graph_destroy(&g2);
}

void build_barcode_process(struct opt_proc_t *opt)
{
	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	build_barcode(opt, &g);
	save_graph_info(opt->out_dir, &g, "added_barcode");
	asm_graph_destroy(&g);
}

void build_barcode_process_fasta(struct opt_proc_t *opt)
{
	struct asm_graph_t g;
	load_asm_graph_fasta(&g, opt->in_fasta, opt->k0);
	build_barcode_read(opt, &g);
	save_graph_info(opt->out_dir, &g, "added_barcode");
	asm_graph_destroy(&g);
}
