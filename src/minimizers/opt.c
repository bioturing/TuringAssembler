//
// Created by BioTuring on 2019-09-21.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "attribute.h"
#include "verbose.h"

char *lib_str[] = {"sorted", "bioturing", "ust", "10x"};
int n_lib = 4;

struct opt_proc_t *init_opt_proc()
{
	struct opt_proc_t *opt;
	opt = calloc(1, sizeof(struct opt_proc_t));
	opt->n_threads = 1;
	opt->hash_size = 1 << 24;
	opt->k0 = 17;
	opt->k1 = 31;
	opt->k2 = 55;
	opt->n_files = 0;
	opt->split_len = 1000;
	opt->files_1 = opt->files_2 = NULL;
	opt->in_file = NULL;
	opt->in_fasta = NULL;
	opt->in_fastg = NULL;
	opt->out_dir = ".";
	opt->lib_type = 0;
	opt->mmem = 32;
	return opt;
}

int get_library_index(const char *str)
{
	int i;
	for (i = 0; i < n_lib; ++i)
		if (!strcmp(str, lib_str[i]))
			return i;
	return -1;
}

void print_info(int argc, char *argv[])
{
	int cmd_len = 0, i;
	for (i = 0; i < argc; ++i)
		cmd_len += strlen(argv[i]) + 1;
	char *cmd = malloc(cmd_len);
	cmd_len = 0;
	for (i = 0; i < argc; ++i)
		cmd_len += sprintf(cmd + cmd_len, i + 1 == argc ? "%s" : "%s ", argv[i]);
	__VERBOSE_LOG("INFO", "Skipping assembly\n");
	__VERBOSE_LOG("INFO", "Version: %s%s\n", VERSION_STRING, GIT_SHA);
	__VERBOSE_LOG("INFO", "command: \"%s\"\n", cmd);
	free(cmd);
}

int opt_count_list(int argc, char *argv[])
{
	int n;
	for (n = 0; n < argc - 1; ++n) {
		if (argv[n + 1][0] == '-')
			break;
	}
	if (n == 0)
		log_error("Emtpy list %s", argv[0]);
	return n;
}

struct opt_proc_t *parse_proc_option(int argc, char *argv[])
{
	int pos = 0, n;
	struct opt_proc_t *opt = init_opt_proc();
	while (pos < argc) {
		if (!strcmp(argv[pos], "-t")) {
			opt->n_threads = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-s")) {
			opt->hash_size = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-k0")) {
			opt->k0 = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-k1")) {
			opt->k1 = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-k2")) {
			opt->k2 = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-o")) {
			opt->out_dir = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-i")) {
			opt->in_file = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-f")) {
			opt->in_fasta = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-fg")) {
			opt->in_fastg = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-bx")) {
			opt->bx_str = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-cf")) {
			opt->in_contig_file = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-metagenomics")){
			opt->metagenomics = 1;
			pos += 1;
		} else if (!strcmp(argv[pos], "-l")) {
			opt->lib_type = get_library_index(argv[pos + 1]);
			if (opt->lib_type == -1)
				log_error("Unknown library %s", argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-sm")) {
			opt->mmem = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-1")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				log_error("Inconsistent number of files");
			opt->n_files = n;
			opt->files_1 = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-2")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				log_error("Inconsistent number of files");
			opt->n_files = n;
			opt->files_2 = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-I")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				log_error("Inconsistent number of files");
			opt->n_files = n;
			opt->files_I = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-sl")) {
			opt->split_len = atoi(argv[pos + 1]);
			pos += 2;
		} else if (argv[pos][0] != '-') {
			if (opt->n_files != 0)
				log_error("Unknown %s", argv[pos]);
			opt->files_1 = argv + pos;
			while (pos < argc && argv[pos][0] != '-') {
				++pos;
				++opt->n_files;
			}
		} else if (!strcmp(argv[pos], "-lc")){
			opt->lc = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-lk")){
			opt->lk = atoi(argv[pos + 1]);
			pos += 2;
		} else {
			log_error("Unknown option %s", argv[pos]);
		}
	}
	mkdir(opt->out_dir, 0755);
	return opt;
}

void build_opt_process(int argc, char *argv[], void (*build_process)(struct opt_proc_t *))
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (opt == NULL) {
		log_error("Error parsing arguments");
	}
	char tmp_dir[1024];
	snprintf(tmp_dir, 1024, "%s/assembly.log", opt->out_dir);
	init_log(tmp_dir);
	init_clock();
	print_info(argc, argv);
	build_process(opt);
}
