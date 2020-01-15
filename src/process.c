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
#include "minimizers/minimizers.h"
#include "minimizers/smart_load.h"
#include "minimizers/count_barcodes.h"
#include "minimizers/get_buffer.h"
#include "coverage/kmer_count.h"
#include "cluster_molecules.h"
#include "split_molecules.h"
#include "barcode_graph.h"
#include "read_pairs_resolve.h"

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
	asm_resolve_simple_bulges_ite(g);
	asm_resolve_complex_bulges_ite(g);
	test_asm_graph(g);
	log_info("Build graph level 1 time: %.3f", sec_from_prev_time());
}

int64_t dfs(struct asm_graph_t *g, int x, int *mark, int *hd)
{
	int64_t res = g->edges[x].seq_len-g->ksize;
	int l = 0, r = 1;
	hd[0] = x;
	while (l < r) {
		int x = hd[l++];
		res += g->edges[x].seq_len-g->ksize;
		int target = g->edges[x].target;
		assert(target >=0 && target < g->n_v);
		int src_rc = g->nodes[target].rc_id;
		for (int i = 0; i < g->nodes[src_rc].deg; i++) {
			int next = g->nodes[src_rc].adj[i];
			if (mark[next] == 0) {
				mark[next] = 1;
				hd[r++] =  next;
			}
		}
		for (int i = 0; i < g->nodes[target].deg; i++) {
			int next = g->nodes[target].adj[i];
			if (mark[next] == 0) {
				mark[next] = 1;
				hd[r++] =  next;
			}
		}
	}
	return res;
}

void count_cc(struct asm_graph_t *g)
{
	int *mark = calloc(g->n_e, sizeof(int));
	int MAX = 35;
	int *res = calloc(MAX, sizeof(int)), big_res = 0;
	int *hd = calloc(100000000, 4), total_cc = 0, total_10k_cc = 0;
	for (int i = 0; i < g->n_e; i++) if (mark[i] == 0){
		if (g->edges[i].source == -1)
			continue;
		int len = dfs(g, i, mark, hd);
		if (len > MAX*10000) {
			log_info("This component is big: %d", len);
		} else {
			res[len/10000]++;
		}
		if (len > 10000)
			total_10k_cc++;
		total_cc++;
	}
	for (int i = 0 ;  i < MAX; i++) {
		log_info("Number of CC of len %d is %d", i, res[i]);
	}
	log_info("CC bigger than 10k is: %d", total_10k_cc);
	log_info("Total cc %d", total_cc);
}


void build_0_1(struct opt_proc_t *opt, struct asm_graph_t *g0, struct asm_graph_t *g)
{
	char *log_file = str_concate(opt->out_dir, "/build_0_1.log");
	init_logger(opt->log_level, log_file);
	init_clock();
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

void write_neo4j_create(struct asm_graph_t *g)
{
	FILE *out = fopen("/home/che/bioturing/data/yeast/neo4j/listnode.tsv", "w");
	fprintf(out, "ID\n");
	for (int i = 0; i < g->n_e; i++){
		fprintf(out, "%d\n", i);
	}
	fclose(out);
	out = fopen("/home/che/bioturing/data/yeast/neo4j/reverse_edge.tsv", "w");
	fprintf(out, "i,i_rc\n");
	for (int i = 0; i < g->n_e; i++){
		int rc_id = g->edges[i].rc_id;
		if (rc_id < i)
			continue;
		fprintf(out, "%d,%d\n", i, rc_id);
	}
	fclose(out);
}

void dirty_process(struct opt_proc_t *opt)
{
	struct asm_graph_t *g = create_and_load_graph(opt);
	struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
	asm_condense(g, g0);
//	count_cc(g);
//	dirty(g, opt);
//	write_neo4j_create(g);
}

void resolve_212_cov_process(struct opt_proc_t *opt)
{
	struct asm_graph_t *g = create_and_load_graph(opt);
	resolve_212_by_cov(g, opt);
}

void resolve_molecule_process(struct opt_proc_t *opt)
{
	struct asm_graph_t *g = create_and_load_graph(opt);
	get_long_contig(g, opt);
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

void split_molecules_wrapper(struct opt_proc_t *opt)
{
	char path[1024];
	sprintf(path, "%s/debug.log", opt->out_dir);
	init_logger(opt->log_level, path);
	set_log_stage("Split molecules");
	FILE *f = fopen(opt->in_fasta, "r");
	char bx[19];
	int count;
	fclose(fopen(opt->lc, "w"));
	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	struct mm_db_edge_t *mm_edges = mm_index_edges(&g, MINIMIZERS_KMER, MINIMIZERS_WINDOW);
	khash_t(bcpos) *bx_pos_dict = kh_init(bcpos);
	struct read_path_t read_sorted_path;
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
		log_info("Reads are not sorted. Sort reads by barcode sequence...");
		sort_read(opt, &read_sorted_path);
	}
	smart_construct_read_index(&read_sorted_path, bx_pos_dict); //load the barcode indices
	int C = 0;
	while (fscanf(f, "%s\t%d\n", bx, &count)){
		if (C == 50000)
			break;
		if ((++C) % 1000 == 0)
			log_debug("Processing %d-th barcode", C);
		opt->bx_str = bx;
		split_molecules_process(opt, &g, mm_edges, bx_pos_dict);
	}
	fclose(f);
}

/**
 *
 * @param opt
 * @param g
 * @param mm_edges
 * @param bx_pos_dict
 */
void split_molecules_process(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct mm_db_edge_t *mm_edges, khash_t(bcpos) *bx_pos_dict)
{

	struct read_path_t read_sorted_path;
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
		log_info("Reads are not sorted. Sort reads by barcode sequence...");
		sort_read(opt, &read_sorted_path);
	}
	uint64_t bx_encoded = barcode_hash_mini(opt->bx_str);
	//log_info("Hashed barcode: %lu", bx_encoded);
	uint64_t bx[1] = {bx_encoded}; //43 15 mock barcode pseudo hash id here

	khint_t k = kh_get(bcpos, bx_pos_dict, bx_encoded);          // query the hash table
	if (k == kh_end(bx_pos_dict)) {
		log_error("Barcode does not exists!");
	} else {
		//log_info("Barcode does exist. Getting reads");
	}

	char *buf1, *buf2;
	uint64_t m_buf1, m_buf2;
	stream_filter_read(&read_sorted_path, bx_pos_dict, bx, 1, &buf1, &buf2, &m_buf1, &m_buf2);

	struct read_t r1, r2;
	int pos1 = 0, pos2 = 0;
	int n_reads = 0;
	struct mm_db_t *db1, *db2;
	struct mm_hits_t *hits;
	hits = mm_hits_init();


	while (get_read_from_fq(&r1, buf1, &pos1) == READ_SUCCESS && get_read_from_fq(&r2, buf2, &pos2) == READ_SUCCESS) {
		n_reads++;
		db1 = mm_index_char_str(r1.seq, MINIMIZERS_KMER, MINIMIZERS_WINDOW, r1.len);
		db2 = mm_index_char_str(r2.seq, MINIMIZERS_KMER, MINIMIZERS_WINDOW, r2.len);

		mm_hits_cmp(db1, mm_edges, hits, g);
		mm_hits_cmp(db2, mm_edges, hits, g);
	}
	//log_info("Number of read-pairs in barcode %s: %d", opt->bx_str, n_reads);
	//log_info("Number of singleton hits: %d", kh_size(hits->edges));
	mm_hits_print(hits, "barcode_hits.csv");
	free(buf1);
	free(buf2);

	order_edges(opt, g, hits);
}

void debug_process(struct opt_proc_t *opt)
{
	char path[1024];
	sprintf(path, "%s/debug.log", opt->out_dir);
	init_logger(opt->log_level, path);
	set_log_stage("Debug process");
	get_long_contigs_by_readpairs(opt);
	log_info("get long contig by readpair done");
}

void read_pairs_count_process(struct opt_proc_t *opt)
{
	char path[1024];
	sprintf(path, "%s/debug.log", opt->out_dir);
	init_logger(opt->log_level, path);
	set_log_stage("Debug process");
	khash_t(long_int) *rp_count = kh_init(long_int);
	get_all_read_pairs_count(opt, rp_count);
	sprintf(path, "%s/rp_counts.txt", opt->out_dir);
	FILE *f = fopen(path, "w");
	for (khiter_t it = kh_begin(rp_count); it != kh_end(rp_count); ++it){
		if (!kh_exist(rp_count, it))
			continue;
		uint64_t code = kh_key(rp_count, it);
		int v = code >> 32;
		int u = code & -1;
		assert(v >= 0);
		assert(u >= 0);
		fprintf(f, "%d %d %d\n", v, u, kh_val(rp_count, it));
	}
	fclose(f);
	opt->in_fasta = calloc(1024, 1);
	COPY_ARR(path, opt->in_fasta, strlen(path)+1);
}
void print_barcode_graph_process(struct opt_proc_t *opt)
{
	print_barcode_graph(opt);
}

//void cluster_molecules_process(struct opt_proc_t *opt)
//{
//	char path[1024];
//	sprintf(path, "%s/cluster_molecules.log", opt->out_dir);
//	init_logger(opt->log_level, path);
//	set_log_stage("Cluster molecules");
//
//	count_edge_links_bc(opt);
//}

void resolve_multi_kmer_process(struct opt_proc_t *opt)
{
	char path[1024];
	sprintf(path, "%s/resolve_ite_kmer.log", opt->out_dir);
	init_logger(opt->log_level, path);
	set_log_stage("Resolve iterative kmers");
	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	resolve_multi_kmer(opt, &g, opt->k1);
}

void resolve_complex_bulges_process(struct opt_proc_t *opt)
{
	char path[1024];
	sprintf(path, "%s/resolve_bulges.log", opt->out_dir);
	init_logger(opt->log_level, path);
	set_log_stage("Resolve complex bulges");
	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	asm_resolve_complex_bulges_ite(&g);

	save_graph_info(opt->out_dir, &g, "no_complex_bulges");
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
	asm_resolve_simple_bulges_ite(&g);

	save_graph_info(opt->out_dir, &g, "no_simple_bulges");
	asm_graph_destroy(&g);
}

void index_mm_process(struct opt_proc_t *opt)
{
	__VERBOSE("Index minimizers for an example string\n");
	set_log_stage("Minimizers Index");
	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	mm_index_edges(&g, MINIMIZERS_KMER, MINIMIZERS_WINDOW);
	__VERBOSE("Done indexing\n");
}

void hits_barcode_process(struct opt_proc_t *opt)
{
	__VERBOSE("Hitting reads from one barcode to edges\n");
	set_log_stage("Barcode hits");
	smart_load_barcode(opt);
	__VERBOSE("Done hitting\n");
}

void count_bx_process(struct opt_proc_t *opt)
{
	struct read_path_t read_path;
	__VERBOSE("Counting barcode frequencies\n");
	count_bx_freq(opt);
	__VERBOSE("Done counting\n");
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
	//resolve_1_2(&g0, opt);
	asm_resolve_simple_bulges_ite(&g0);
	asm_resolve_complex_bulges_ite(&g0);
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
	construct_aux_info(opt, &g, &read_sorted_path, fasta_path, ASM_BUILD_BARCODE);
	save_graph_info(opt->out_dir, &g, "added_barcode");
	asm_graph_destroy(&g);
}

void build_barcode_coverage_info(struct opt_proc_t *opt)
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
	construct_aux_info(opt, &g, &read_sorted_path, fasta_path, ASM_BUILD_BARCODE|ASM_BUILD_COVERAGE);
	save_graph_info(opt->out_dir, &g, "added_barcode");
	asm_graph_destroy(&g);
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
	build_0_1(opt, &g_lv0, &g_lv1); /* Simplify graph level 0 to graph level 1 */
	save_graph_info(opt->out_dir, &g_lv1, "level_1");
	if (g_lv1.n_e == 0) {
		log_error("graph after lv1 has 0 edges");
	}
	asm_graph_destroy(&g_lv1);

	/**
	  * Resolve process (dump function, not use yet)
	  */
	set_log_stage("ResolveProcess");
	log_info("Start resolve process with kmer size %d", opt->k0);
	char lv1_path[1024];
	sprintf(lv1_path, "%s/graph_k_%d_level_1.bin", opt->out_dir, opt->k0);
	opt->in_file = lv1_path;
	resolve_local_process(opt);
	char lv2_path[1024];
	sprintf(lv2_path, "%s/graph_k_%d_level_2.bin", opt->out_dir, opt->k0);
	opt->in_file = lv2_path;

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
	 * Get readpair count
	 */
	read_pairs_count_process(opt);

	/**
	 * resolve by readpairs
	 */
	set_log_stage("Get long contig by readpairs");
	get_long_contigs_by_readpairs(opt);

	/**
	 * Build barcodes
	 */
//	char lv3_path[1024];
//	sprintf(lv3_path, "%s/graph_k_%d_level_3.bin", opt->out_dir, opt->k0);
	set_log_stage("BWAIndex");
	build_barcode_process_fasta(opt);

	/**
	* Scaffolding
	*/
	sprintf(opt->in_file, "%s/graph_k_%d_added_barcode.bin", opt->out_dir,
		opt->k0);
	struct asm_graph_t g_lv3_added;
	load_asm_graph(&g_lv3_added, opt->in_file);
	set_log_stage("Scaffolding");
	char out_name[MAX_PATH];
	sprintf(out_name, "%s/scaffolds.fasta", opt->out_dir);
	log_info("Construct the scaffolds using barcode information.");
	FILE *out_file = fopen(out_name, "w");
	scaffolding(out_file, &g_lv3_added, opt);
	fclose(out_file);
	log_info("Done scaffolding. Please see the file: %s/scaffolds.fasta", opt->out_dir);
	asm_graph_destroy(&g_lv3_added);


	/**
	* Local assembly
	*/
	set_log_stage("LocalAssembly");
	log_info("Start local assembly process with kmer size %d", opt->lk);
	char in_fasta[1024];
	char local_fasta[1024];
	sprintf(in_fasta, "%s/local_assembly_scaffold_path.txt", opt->out_dir);
	sprintf(local_fasta, "%s/scaffold.full.fasta", opt->out_dir);
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
	build_0_1(opt, &g1, &g2);
	save_graph_info(opt->out_dir, &g2, "level_1");
	count_cc(&g2);
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
	gint_t le_idx = get_longest_edge(&g);
	if (le_idx != -1) {
		log_info("Longest edge %ld_%ld, length %u",
		         le_idx, g.edges[le_idx].rc_id, get_edge_len(g.edges + le_idx));
	}
	save_graph_info(opt->out_dir, &g, "from_fasta");
	char fasta_path[MAX_PATH];
	struct read_path_t read_path;
	read_path.R1_path = opt->files_1[0];
	read_path.R2_path = opt->files_2[0];
	sprintf(fasta_path, "%s/barcode_build_dir", opt->out_dir);
	mkdir(fasta_path, 0755);
	sprintf(fasta_path, "%s/barcode_build_dir/contigs_tmp.fasta", opt->out_dir);
	write_fasta_seq(&g, fasta_path);
	log_info("Aligning reads on two heads of each contigs using BWA");
	construct_aux_info(opt, &g, &read_path, fasta_path, ASM_BUILD_BARCODE|ASM_BUILD_COVERAGE);
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

void build_coverage_process(struct opt_proc_t *opt)
{
	char *log_file = str_concate(opt->out_dir, "/index_kmer_edges.log");
	init_logger(opt->log_level, log_file);
	set_log_stage("Calculate coverage");
	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	struct mini_hash_t *table = kmer_count_on_edges(opt, &g);
	add_cnt_to_graph(&g, table);
	save_graph_info(opt->out_dir, &g, "coverage_built");
	asm_graph_destroy(&g);
}
