#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "attribute.h"
#include "verbose.h"

void print_usage(const char *prog)
{
	__VERBOSE("Skipping sort reads\n");
	__VERBOSE("Usage: %s\n", prog);
	__VERBOSE("          sort [options] -1 read_1.fq -2 read_2.fq -I index.fq -l ust -o out_dir\n");
	__VERBOSE("          sort [options] -1 read_1_added.fq -2 read_2_added.fq -l bioturing -o out_dir\n");
}

void sort_read_process(struct opt_proc_t *opt)
{
	struct read_path_t read_sorted_path;
	__VERBOSE("Sort reads by barcodes");
	sort_read(opt, &read_sorted_path);
	__VERBOSE("Done sorted");
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		print_usage(argv[0]);
		return -1;
	}
	if (!strcmp(argv[1], "sort"))
		build_opt_process(argc, argv, &sort_read_process);
	else
		return -1;
	return 0;
}
