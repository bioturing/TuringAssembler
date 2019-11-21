#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "attribute.h"
#include "verbose.h"
#include "utils.h"
#include "smart_load.h"
#include "count_barcodes.h"
#include "minimizers.h"
#include "assembly_graph.h"

void print_usage(const char *prog)
{
	__VERBOSE("Skipping sort reads\n");
	__VERBOSE("Usage: %s\n", prog);
	__VERBOSE("          sort [options] -1 read_1.fq -2 read_2.fq -I index.fq -l ust -o out_dir\n");
	__VERBOSE("          count_bx [options] -1 read_1_added.fq -2 read_2_added.fq -l bioturing -o out_dir\n");
	__VERBOSE("          search_bx [options] -1 read_1_added.fq -2 read_2_added.fq -l sorted -I barcode.idx -bx AAAAAA\n");
}

void sort_read_process(struct opt_proc_t *opt)
{
	struct read_path_t read_sorted_path;
	__VERBOSE("Sort reads by barcodes\n");
	sort_read(opt, &read_sorted_path);
	__VERBOSE("Done sorted\n");
}

void count_bx_process(struct opt_proc_t *opt)
{
	struct read_path_t read_path;
	__VERBOSE("Counting barcode frequencies\n");
	count_bx_freq(opt, &read_path);
	__VERBOSE("Done counting\n");
}

void search_bx_process(struct opt_proc_t *opt)
{
	struct read_path_t read_path;
	__VERBOSE("Searching reads from a barcode\n");
	smart_load_barcode(opt);
	__VERBOSE("Done searching\n");
}

void index_mm_process(struct opt_proc_t *opt)
{
	__VERBOSE("Index minimizers for an example string\n");
	uint32_t *s;
	set_log_stage("Minimizers Index");
	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	mm_index_edges(&g, MINIMIZERS_KMER, MINIMIZERS_WINDOW);
	__VERBOSE("Done indexing\n");
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		print_usage(argv[0]);
		return -1;
	}
	if (!strcmp(argv[1], "sort"))
		build_opt_process(argc, argv, &sort_read_process);
	else if (!strcmp(argv[1], "count_bx"))
		build_opt_process(argc, argv, &count_bx_process);
	else if (!strcmp(argv[1], "search_bx"))
		build_opt_process(argc, argv, &search_bx_process);
	else if (!strcmp(argv[1], "index"))
		build_opt_process(argc, argv, &index_mm_process);
	else
		return -1;
	return 0;
}
