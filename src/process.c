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
#include "build_bridge.h"
#include "read_list.h"
#include "scaffolding/scaffolding.h"
#include "barcode_resolve2.h"
#include "basic_resolve.h"
#include "fastg.h"
#include "fastq_reducer.h"
#include "resolve_big.h"
#include "complex_resolve.h"

void graph_convert_process(struct opt_proc_t *opt)
{
	char path[1024];
	log_info("Dump graph from bin archive");
	struct asm_graph_t *g;
	g = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g, opt->in_file);
	log_info("Input graph kmer size: %d", g->ksize);
	log_info("kmer size: %d", g->ksize);
	test_asm_graph(g);
	snprintf(path, 1024, "%s/graph_k_%d_loaded.gfa", opt->out_dir, g->ksize);
	write_gfa(g, path);
	snprintf(path, 1024, "%s/graph_k_%d_loaded.fasta", opt->out_dir, g->ksize);
	write_fasta(g, path);
}

void build_0_KMC(struct opt_proc_t *opt, int ksize, struct asm_graph_t *g)
{
	log_info("Building assembly graph from read using kmer size %d", ksize);
	build_initial_graph(opt, ksize, g);
	// graph_build_KMC(opt, ksize, g);
	test_asm_graph(g);
}

void build_local_0_1(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	log_info("Resolve graph using small operation");
	log_info("Input graph kmer size: %d", g0->ksize);
	set_time_now();
	resolve_local_graph_operation(g0, g);
	// remove_tips(g0, g);
	test_asm_graph(g);
	log_info("Build graph level 1 time: %.3f", sec_from_prev_time());
}

void build_0_1(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	log_info("Resolve graph using small operation");
	log_info("Input graph kmer size: %d", g0->ksize);
	set_time_now();
	resolve_graph_operation(g0, g);
	// remove_tips(g0, g);
	test_asm_graph(g);
	log_info("Build graph level 1 time: %.3f", sec_from_prev_time());
}

void build_1_2b(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	log_info("Resolving graph using barcode (local assembly included)");
	log_info("Input graph kmer size: %d", g0->ksize);
	set_time_now();
	// resolve_tips_topo(g0, g);
	// remove_tips_topology(g0, g);
	// test_asm_graph(g);
	log_info("Build graph level 2 time: %.3f", sec_from_prev_time());

}

void build_1_2(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	log_info("Removing tips using graph topology");
	log_info("Input graph kmer size: %d", g0->ksize);
	set_time_now();
	// resolve_tips_topo(g0, g);
	// remove_tips_topology(g0, g);
	test_asm_graph(g);
	log_info("Build graph level 2 time: %.3f", sec_from_prev_time());
}

struct asm_graph_t* create_and_load_graph(struct opt_proc_t *opt)
{
	struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_file);
	test_asm_graph(g0);
	return g0;
}

void build_scaffolding_1_2_process(struct opt_proc_t *opt)
{
	char *log_file = str_concate(opt->out_dir, "/build_scaffolding_1_2.log");
	init_logger(opt->log_level, log_file);
	init_clock();
	struct asm_graph_t *g0 = create_and_load_graph(opt);
	log_info("Build scaffolding with kmer size: %d", g0->ksize);

	char *out_name = str_concate(opt->out_dir, "/scaffolds.fasta");
	FILE *out_file = fopen(out_name, "w");
	
	scaffolding(out_file, g0, opt);

	free(out_name);
	fclose(out_file);
	asm_graph_destroy(g0);
	free(g0);
	free(log_file);
	close_logger();
}

void build_scaffolding_test_process(struct opt_proc_t *opt)
{
	init_clock();
	FILE *fp;
	struct asm_graph_t *g0 = create_and_load_graph(opt);

	scaffolding_test(g0, opt);

	asm_graph_destroy(g0);
}

void build_2_3(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	log_info("Resolving small complex structure");
	log_info("Input graph kmer size: %d", g0->ksize);
	set_time_now();
	// resolve_chain(g0, g);
	test_asm_graph(g);
	log_info("Build graph level 3 time: %.3f", sec_from_prev_time());
}

void build_3_4(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	log_info("Resolving graph using read pair + barcode (simple repetitive)");
	log_info("Input graph kmer size: %d", g0->ksize);
	set_time_now();
	resolve_n_m_simple(g0, g);
	test_asm_graph(g);
	log_info("Build graph level 4 time: %.3f", sec_from_prev_time());
}

void build_4_5(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	log_info("Resolving graph using read pair + barcode (complex jungle)");
	log_info("Input graph kmer size: %d", g0->ksize);
	set_time_now();
	resolve_complex(g0, g);
	test_asm_graph(g);
	log_info("Build graph level 5 time: %.3f", sec_from_prev_time());
}

void build_barcode(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	log_info("Building barcode information");
	log_info("Input graph kmer size: %d", g->ksize);
	set_time_now();
	// construct_aux_information(opt, g, ASM_BUILD_BARCODE | ASM_BUILD_READPAIR);
	log_info("Build barcode information time: %.3f", sec_from_prev_time());
}

void build_barcode_read(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	log_info("Building barcode information");
	log_info("Input graph kmer size: %d", g->ksize);
	set_time_now();
	// construct_aux_information(opt, g, ASM_BUILD_BARCODE | ASM_BUILD_READPAIR | ASM_BUILD_COVERAGE);
	log_info("Build barcode information time: %.3f", sec_from_prev_time());
}

void build_bridge_process(struct opt_proc_t *opt)
{
	build_bridge(opt);
}

void resolve_complex_bulges_process(struct opt_proc_t *opt)
{
	char path[1024];
	sprintf(path, "%s/resolve_bulges.log", opt->out_dir);
	init_logger(opt->log_level, path);
	set_log_stage("Resolve complex bulges");
	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	asm_resolve_complex_bulges_ite(opt, &g);

	save_graph_info(opt->out_dir, &g, "level_3");
	asm_graph_destroy(&g);
}

void resolve_bulges_process(struct opt_proc_t *opt)
{
	char path[1024];
	sprintf(path, "%s/resolve_bulges.log", opt->out_dir);
	init_logger(opt->log_level, path);
	set_log_stage("Resolve bulges");
	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	asm_resolve_simple_bulges_ite(opt, &g);

	save_graph_info(opt->out_dir, &g, "level_2");
	asm_graph_destroy(&g);
}

void reduce_read_process(struct opt_proc_t *opt)
{
	struct read_path_t org_rpath;
	char path[1024];
	sprintf(path, "%s", opt->files_1[0]);
	org_rpath.R1_path = strdup(path);
	sprintf(path, "%s", opt->files_2[0]);
	org_rpath.R2_path = strdup(path);


	struct read_path_t reduced_rpath;
	sprintf(path, "%s/R1.added_barcode.reduced.fastq", opt->out_dir);
	reduced_rpath.R1_path = strdup(path);
	sprintf(path, "%s/R2.added_barcode.reduced.fastq", opt->out_dir);
	reduced_rpath.R2_path = strdup(path);
	__VERBOSE("REDUCING READS\n");
	fastq_reducer(opt, &org_rpath, &reduced_rpath);
}

/*
 * @brief: main function for resolving graph before scaffolding
 * @param opt: application options
 */
void resolve_local_process(struct opt_proc_t *opt)
{
	char log_file[1024];
	struct asm_graph_t g0;
	load_asm_graph(&g0, opt->in_file);
	int resolved_loop = 0;
	/*asm_resolve_dump_loop_ite(&g0);
	asm_resolve_dump_jungle_ite(opt, &g0);
	do_some_resolve_bridge(&g0);*/
	asm_resolve_simple_bulges_ite(opt, &g0);
	asm_resolve_complex_bulges_ite(opt, &g0);
	char path[1024];
	sprintf(path, "%s/graph_k_%d_level_2.bin", opt->out_dir, g0.ksize);
	struct asm_graph_t g1;
	//asm_condense(&g0, &g1); // Remove barcode
	asm_condense(&g0, &g1);
	asm_graph_destroy(&g0);

	save_graph_info(opt->out_dir, &g1, "level_2");
	asm_graph_destroy(&g1);
}

void save_graph_info(const char *out_dir, struct asm_graph_t *g, const char *suffix)
{
	char path[1024];
	log_info("graph_k_%d_%s: Number of nodes: %ld",
						g->ksize, suffix, g->n_v);
	log_info("graph_k_%d_%s: Number of edges: %ld",
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
	construct_aux_info(opt, &g, &read_sorted_path, fasta_path, ASM_BUILD_BARCODE, NOT_FOR_SCAFF);
	save_graph_info(opt->out_dir, &g, "added_barcode");
	asm_graph_destroy(&g);
}

void build_barcode_scaffold(struct opt_proc_t *opt)
{
	char *log_file = str_concate(opt->out_dir, "/build_barcode_scaffold.log");
	init_logger(opt->log_level, log_file);
	struct asm_graph_t g;
	struct read_path_t read_path;

	load_asm_graph(&g, opt->in_file);
	char fasta_path[MAX_PATH];
//	if (opt->lib_type == LIB_TYPE_SORTED) {
//		read_sorted_path.R1_path = opt->files_1[0];
//		read_sorted_path.R2_path = opt->files_2[0];
//		read_sorted_path.idx_path = opt->files_I[0];
//	} else {
//		log_info("Read library is not sorted (type %d). Sorting reads by barcodes", opt->lib_type);
//		sort_read(opt, &read_sorted_path);
//	}
	if (opt->lib_type == LIB_TYPE_BIOT || opt->lib_type== LIB_TYPE_UST) {
		read_path.R1_path = opt->files_1[0];
		read_path.R2_path = opt->files_2[0];
	} else{
		__ERROR("MUST BE SORTED OR UST TYPE");
	}
	sprintf(fasta_path, "%s/barcode_build_dir", opt->out_dir);
	mkdir(fasta_path, 0755);
	sprintf(fasta_path, "%s/barcode_build_dir/contigs_tmp.fasta", opt->out_dir);
	write_fasta_seq(&g, fasta_path);
	log_info("Aligning reads on two heads of each contigs using BWA");
	construct_aux_info(opt, &g, &read_path, fasta_path, ASM_BUILD_BARCODE, FOR_SCAFFOLD);
	save_graph_info(opt->out_dir, &g, "added_barcode");
	asm_graph_destroy(&g);
	close_logger();
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

/**
 * @brief Main process for the assembly. Any step further please add into this.
 * Only add after a pull request acceptance
 * @param opt Main option structure of the process
 */
void assembly3_process(struct opt_proc_t *opt)
{
	char *log_file = str_concate(opt->out_dir, "/assembly3.log");
	init_logger(opt->log_level, log_file);
	set_log_stage("General");
	log_info("Version: %s\n", GIT_SHA);
	struct read_path_t read_sorted_path; /* To sort fastq file if needed */

	/**
	 * Build Assembly graph from scratch
	 */
	struct asm_graph_t g_lv0, g_lv1;
	set_log_stage("KmerCounting");
	build_0_KMC(opt, opt->k0, &g_lv0); /* Do kmer counting using KMC */
	save_graph_info(opt->out_dir, &g_lv0, "level_0");

	set_log_stage("GraphConstruction");
	build_0_1(&g_lv0, &g_lv1); /* Simplify graph level 0 to graph level 1 */
	save_graph_info(opt->out_dir, &g_lv1, "level_1");
	if (g_lv1.n_e == 0) {
		log_error("graph after lv1 has 0 edges");
	}
	asm_graph_destroy(&g_lv1);

	/**
	 * Rearrange reads in fastq files. Reads from the same barcodes are grouped together
	 */
	char fasta_path[MAX_PATH];
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
		log_info("Read library is not sorted (type %d). Sorting reads by barcodes", opt->lib_type);
		set_log_stage("SortReads");
		sort_read(opt, &read_sorted_path);
		opt->lib_type = LIB_TYPE_SORTED;
		opt->n_files = 1;
		opt->files_1 = alloca(sizeof(char *));
		opt->files_2 = alloca(sizeof(char *));
		opt->files_I = alloca(sizeof(char *));
		opt->files_1[0] = read_sorted_path.R1_path;
		opt->files_2[0] = read_sorted_path.R2_path;
		opt->files_I[0] = read_sorted_path.idx_path;
	}

	/**
	  * Resolve process
	  */
	set_log_stage("Resolve process");
	log_info("Start resolve process with kmer size %d", opt->lk);
	char lv1_path[1024];
	sprintf(lv1_path, "%s/graph_k_%d_level_1.bin", opt->out_dir, opt->k0);
	opt->in_file = lv1_path;
	resolve_local_process(opt);

	/**
	 * Use bwa to align reads on two ends of each contigs, barcode awared.
	 */
	char lv2_path[1024];
	sprintf(lv2_path, "%s/graph_k_%d_level_2.bin", opt->out_dir, opt->k0);
	struct asm_graph_t g_lv2;
	load_asm_graph(&g_lv2, lv2_path);
	set_log_stage("BWAIndex");
	sprintf(fasta_path, "%s/barcode_build_dir", opt->out_dir); /* Store temporary contigs for indexing two heads */
	mkdir(fasta_path, 0755);
	sprintf(fasta_path, "%s/barcode_build_dir/contigs_tmp.fasta", opt->out_dir);
	log_info("Write down temporary contigs for indexing.");
	write_fasta_seq(&g_lv2, fasta_path);
	log_info("Aligning reads on two heads of each contigs using BWA");
	set_log_stage("MapReads");
	construct_aux_info(opt, &g_lv2, &read_sorted_path, fasta_path, ASM_BUILD_BARCODE, FOR_SCAFFOLD);
	log_info("Done alignment. Serializing assembly graph to disk.");
	save_graph_info(opt->out_dir, &g_lv2, "added_barcode"); /* Barcode-added assembly graph */
	asm_graph_destroy(&g_lv2);

	/**
	* Scaffolding
	*/
	char lv2_scaff_added_path[1024];
	sprintf(lv2_scaff_added_path, "%s/graph_k_%d_added_barcode.bin", opt->out_dir,
			opt->k0);
	struct asm_graph_t g_scaff_lv2_added;
	load_asm_graph(&g_scaff_lv2_added, lv2_scaff_added_path);
	//test_sort_read(&read_sorted_path, &g_scaff_lv2_added);
	set_log_stage("Scaffolding");
	char out_name[MAX_PATH];
	sprintf(out_name, "%s/scaffolds.fasta", opt->out_dir);
	log_info("Construct the scaffolds using barcode information.");
	FILE *out_file = fopen(out_name, "w");
	scaffolding(out_file, &g_scaff_lv2_added, opt);
	fclose(out_file);
	log_info("Done scaffolding. Please see the file: %s/scaffolds.fasta", opt->out_dir);
	asm_graph_destroy(&g_scaff_lv2_added);


	/**
	 * Build barcode for local assembly
	 */
	set_log_stage("BWAIndex");
	load_asm_graph(&g_lv2, lv2_path);
	char asm_path[1024];
	sprintf(fasta_path, "%s/barcode_build_dir_local", opt->out_dir);
	mkdir(fasta_path, 0755);
	sprintf(fasta_path, "%s/barcode_build_dir_local/contigs_tmp.fasta", opt->out_dir);
	write_fasta_seq(&g_lv2, fasta_path);
	set_log_stage("MapReads");
	construct_aux_info(opt, &g_lv2, &read_sorted_path, fasta_path,
			ASM_BUILD_BARCODE, NOT_FOR_SCAFF);
	save_graph_info(opt->out_dir, &g_lv2, "local_added_barcode");
	asm_graph_destroy(&g_lv2);



	/**
	* Local assembly
	*/
	char lv2_local_added_path[1024];
	sprintf(lv2_local_added_path, "%s/graph_k_%d_local_added_barcode.bin",
			opt->out_dir, opt->k0);
	set_log_stage("LocalAssembly");
	log_info("Start local assembly process with kmer size %d", opt->lk);
	char in_fasta[1024];
	char local_fasta[1024];
	sprintf(in_fasta, "%s/local_assembly_scaffold_path.txt", opt->out_dir);
	sprintf(local_fasta, "%s/scaffold.full.fasta", opt->out_dir);
	opt->in_file = lv2_local_added_path;
	opt->in_fasta = in_fasta;
	opt->lc = local_fasta;
	char local_path[1024];
	sprintf(local_path, "%s/local", opt->out_dir);
	mkdir(local_path, 0755);
	opt->out_dir = local_path;
	build_bridge_process(opt);
	log_info("Done local assembly. Please see the file: %s", opt->lc);

	close_logger();
	free(log_file);
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
	do_some_resolve_bridge(&g1);
	do_something_local(opt, &g1);
	test_asm_graph(&g1);
	save_graph_info(opt->out_dir, &g1, "level_pro");
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

void build_barcode_process_fastg(struct opt_proc_t *opt)
{
	struct asm_graph_t g, g1;
	load_asm_graph_fastg(&g, opt->in_fastg, opt->k0);
	test_asm_graph(&g);
	build_barcode_read(opt, &g);
	build_3_4(&g, &g1);
	save_graph_info(opt->out_dir, &g1, "level_4");
	asm_graph_destroy(&g);
	asm_graph_destroy(&g1);
}

void resolve_1_2_process(struct opt_proc_t *opt)
{
	char *log_file = str_concate(opt->out_dir, "/resolve_1_2.log");
	init_logger(opt->log_level, log_file);
	init_clock();
	struct asm_graph_t *g0 = create_and_load_graph(opt);
	log_info("resolve 1_2 process");
	resolve_1_2(g0,opt);
	save_graph_info(opt->out_dir, g0, "level_2_huu");
	asm_graph_destroy(g0);
	free(g0);
	free(log_file);
	close_logger();
}
