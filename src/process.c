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
	init_logger(opt->log_level, "./build_scaffolding_1_2.log");
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

void graph_query_process(struct opt_proc_t *opt)
{
	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_file);
	fprintf(stderr, "bin size = %d\n", g0->bin_size);
	test_asm_graph(g0);
	log_info("Graph query kmer size: %d", g0->ksize);
	log_info("Querying pair-edge barcode information");
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
			log_info("Building local graph");
			fscanf(fp, "%ld\n", &v);
			struct asm_graph_t lg = test_local_assembly(opt, g0,
							g0->edges[u].rc_id, v);
			asm_graph_destroy(&lg);
		} else if (c == 'S') {
			fscanf(fp, "%ld %ld\n", &v, &v2);
			print_test_barcode_superior(g0, u, v, v2);
		}
	}
}

void build_bridge_process(struct opt_proc_t *opt)
{
	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g0, opt->in_file);
	log_debug("bin size = %d", g0->bin_size);
	test_asm_graph(g0);
	log_info("Building bridges with local kmer size: %ld", g0->ksize);
	log_info("Building bridges on scaffold:");
	FILE *f = xfopen(opt->lc, "w");
	FILE *fp = xfopen(opt->in_fasta, "r");
	int n_paths;
	fscanf(fp, "%d\n", &n_paths);
	int *mark = (int *) calloc(g0->n_e, sizeof(int));
	for (int i = 0; i < n_paths; ++i){
		log_info("SCAFFOLD PATH: Processing %d on %d paths", i + 1, n_paths);
		int path_len;
		fscanf(fp, "%d\n", &path_len);
		int *path = (int *) calloc(path_len, sizeof(int));
		for (int i = 0; i < path_len; ++i){
			fscanf(fp, "%d", path + i);
			mark[path[i]] = 1;
			mark[g0->edges[path[i]].rc_id] = 1;
		}
		char *contig;
		get_contig_from_scaffold_path(opt, g0, path, path_len, &contig);
		fprintf(f, ">contig_path_%d\n", i);
		fprintf(f, "%s\n", contig);
		free(contig);
		free(path);
	}
	for (int i = 0; i < g0->n_e; ++i){
		if (g0->edges[i].seq_len < MIN_OUTPUT_CONTIG_LEN)
			continue;
		if (mark[i] == 0){
			int rc = g0->edges[i].rc_id;
			char *tmp;
			decode_seq(&tmp, g0->edges[i].seq, g0->edges[i].seq_len);
			fprintf(f, ">%d_%d\n", i, rc);
			fprintf(f, "%s\n", tmp);
			free(tmp);
			mark[rc] = 1;
		}
	}
	fclose(f);
	fclose(fp);
	free(mark);
	asm_graph_destroy(g0);
	free(g0);
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
	init_logger(opt->log_level, "./build_barcode_scaffold.log");
	struct asm_graph_t g;
	struct read_path_t read_sorted_path;

	load_asm_graph(&g, opt->in_file);
	char fasta_path[MAX_PATH];
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
		log_info("Read library is not sorted (type %d). Rearranging reads by barcodes", opt->lib_type);
		sort_read(opt, &read_sorted_path);
	}
	sprintf(fasta_path, "%s/barcode_build_dir", opt->out_dir);
	mkdir(fasta_path, 0755);
	sprintf(fasta_path, "%s/barcode_build_dir/contigs_tmp.fasta", opt->out_dir);
	write_fasta_seq(&g, fasta_path);
	log_info("Aligning reads on two heads of each contigs using BWA");
	construct_aux_info(opt, &g, &read_sorted_path, fasta_path, ASM_BUILD_BARCODE, FOR_SCAFFOLD);
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
	init_logger(opt->log_level, "./assembly3.log");
	struct read_path_t read_sorted_path; /* To sort fastq file if needed */

	/**
	 * Build Assembly graph from scratch
	 */
	struct asm_graph_t g1, g2;
	build_0_KMC(opt, opt->k0, &g1); /* Do kmer counting using KMC */
	save_graph_info(opt->out_dir, &g1, "level_0");

	build_0_1(&g1, &g2); /* Simplify g1 to g2 */
	save_graph_info(opt->out_dir, &g2, "level_1");

	/**
	 * Rearrange reads in fastq files. Reads from the same barcodes are grouped together
	 */
	char fasta_path[MAX_PATH];
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
		log_info("Read library is not sorted (type %d). Rearranging reads by barcodes", opt->lib_type);
		sort_read(opt, &read_sorted_path);
	}

	/**
	 * Use bwa to align reads on two ends of each contigs, barcode awared.
	 */
	sprintf(fasta_path, "%s/barcode_build_dir", opt->out_dir); /* Store temporary contigs for indexing two heads */
	mkdir(fasta_path, 0755);
	sprintf(fasta_path, "%s/barcode_build_dir/contigs_tmp.fasta", opt->out_dir);
	log_info("Write down temporary contigs for indexing.");
	write_fasta_seq(&g2, fasta_path);
	log_info("Aligning reads on two heads of each contigs using BWA");
	construct_aux_info(opt, &g2, &read_sorted_path, fasta_path, ASM_BUILD_BARCODE, FOR_SCAFFOLD);
	log_info("Done alignment. Serializing assembly graph to disk.");
	save_graph_info(opt->out_dir, &g2, "added_barcode"); /* Barcode-added assembly graph */

	/**
	 * Scaffolding
	 */
	char out_name[MAX_PATH];
	sprintf(out_name, "%s/scaffolds.fasta", opt->out_dir);
	log_info("Construct the scaffolds using barcode information.");
	FILE *out_file = fopen(out_name, "w");
	scaffolding(out_file, &g2, opt);
	fclose(out_file);
	log_info("Done scaffolding. The assembly process is completed.");

	asm_graph_destroy(&g2);
	close_logger();
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
