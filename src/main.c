#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "assembly_graph.h"
#include "get_buffer.h"
#include "io_utils.h"
#include "khash.h"
#include "kmer_count.h"
#include "process.h"
#include "time_utils.h"
#include "utils.h"
#include "verbose.h"

char *lib_str[] = {"ust", "10x"};
int n_lib = 2;

int get_library_index(const char *str)
{
	int i;
	for (i = 0; i < n_lib; ++i)
		if (!strcmp(str, lib_str[i]))
			return i;
	return -1;
}

void print_usage_assembly(const char *prog)
{
	__VERBOSE("Usage: %s assembly [options] -1 read_1.fq -2 read_2.fq\n", prog);
	__VERBOSE("Options: -t                     <number of threads>\n");
	__VERBOSE("         -s                     <pre-alloc size>\n");
	__VERBOSE("         -sl                    <split length of edge sequence>\n");
	__VERBOSE("         -o                     <output directory>\n");
	__VERBOSE("         -k0                    <1st kmer size>\n");
	__VERBOSE("         -k1                    <2nd kmer size>\n");
	__VERBOSE("         -k2                    <3rd kmer size>\n");
	__VERBOSE("         -l                     <lib type [ust, 10x]>\n");
}

void print_usage_build0(const char *prog)
{
	__VERBOSE("Usage: %s build0 [options] -1 read_1.fq -2 read_2.fq\n", prog);
	__VERBOSE("Options: -t                     <number of threads>\n");
	__VERBOSE("         -s                     <pre-alloc size>\n");
	__VERBOSE("         -sl                    <split length of edge sequence>\n");
	__VERBOSE("         -o                     <output directory>\n");
	__VERBOSE("         -k0                    <1st kmer size>\n");
	__VERBOSE("         -k1                    <2nd kmer size>\n");
	__VERBOSE("         -k2                    <3rd kmer size>\n");
	__VERBOSE("         -l                     <lib type [ust, 10x]>\n");
}

void print_usage_build(const char *prog)
{
	__VERBOSE("Usage: %s build_x_y [options] -1 read_1.fq -2 read_2.fq\n", prog);
	__VERBOSE("Options: -t                     <number of threads>\n");
	__VERBOSE("         -s                     <pre-alloc size>\n");
	__VERBOSE("         -o                     <output directory>\n");
	__VERBOSE("         -i                     <input graph>\n");
	__VERBOSE("         -sl                    <split length of edge sequence>\n");
}

void print_usage(const char *prog)
{
	__VERBOSE("Skipping assembly\n");
	__VERBOSE("Version: %s%s\n", VERSION_STRING, GIT_SHA);
	__VERBOSE("Usage: %s\n", prog);
	__VERBOSE("          assembly [options] -1 read_1.fq -2 read_2.fq\n");
	__VERBOSE("          build_x_y [options]\n");
	__VERBOSE("          bin2text [options]\n");
	__VERBOSE("          build_barcode [options] -1 read_1.fq -2 read_2.fq\n");
	__VERBOSE("          build0 [options] -1 read_1.fq -2 read_2.fq\n");
	__VERBOSE("          query -i <input_graph> -f <list edge>\n");
}

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
	opt->out_dir = ".";
	opt->lib_type = 0;
	return opt;
}

struct opt_count_t *init_opt_count()
{
	struct opt_count_t *opt;
	opt = calloc(1, sizeof(struct opt_count_t));
	opt->n_threads = 1;
	opt->hash_size = 1 << 24;
	opt->k0 = 17;
	opt->k1 = 31;
	opt->k2 = 55;
	opt->n_files = 0;
	opt->split_len = 1000;
	opt->files_1 = opt->files_2 = NULL;
	opt->out_dir = ".";
	return opt;
}

struct opt_build_t *init_opt_build()
{
	struct opt_build_t *opt;
	opt = calloc(1, sizeof(struct opt_build_t));
	opt->n_threads = 1;
	opt->hash_size = 1 << 24;
	opt->out_dir = ".";
	opt->in_path = NULL;
	opt->in_file = NULL;
	opt->split_len = 1000;
	opt->files_1 = opt->files_2 = NULL;
	opt->out_dir = ".";
	return opt;
}

int opt_count_list(int argc, char *argv[])
{
	int n;
	for (n = 0; n < argc - 1; ++n) {
		if (argv[n + 1][0] == '-')
			break;
	}
	if (n == 0)
		__ERROR("Emtpy list %s", argv[0]);
	return n;
}

struct opt_build_t *parse_build_option(int argc, char *argv[])
{
	int pos = 0, n;
	struct opt_build_t *opt = init_opt_build();
	while (pos < argc) {
		if (!strcmp(argv[pos], "-t")) {
			opt->n_threads = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-o")) {
			opt->out_dir = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-i")) {
			opt->in_file = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-s")) {
			opt->hash_size = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-sl")) {
			opt->split_len = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-1")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__ERROR("Inconsistent number of files");
			opt->n_files = n;
			opt->files_1 = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-2")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__ERROR("Inconsistent number of files");
			opt->n_files = n;
			opt->files_2 = argv + pos + 1;
			pos += (n + 1);
		} else {
			__ERROR("Unknown option %s", argv[pos]);
		}
	}
	if (opt->in_path == NULL) {
		free(opt);
		return NULL;
	}
	mkdir(opt->out_dir, 0755);
	return opt;
}

struct opt_count_t *parse_count_option(int argc, char *argv[])
{
	int pos = 0, n;
	struct opt_count_t *opt = init_opt_count();
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
		} else if (!strcmp(argv[pos], "-1")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__ERROR("Inconsistent number of files");
			opt->n_files = n;
			opt->files_1 = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-2")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__ERROR("Inconsistent number of files");
			opt->n_files = n;
			opt->files_2 = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-sl")) {
			opt->split_len = atoi(argv[pos + 1]);
			pos += 2;
		} else if (argv[pos][0] != '-') {
			if (opt->n_files != 0)
				__ERROR("Unknown %s", argv[pos]);
			opt->files_1 = argv + pos;
			while (pos < argc && argv[pos][0] != '-') {
				++pos;
				++opt->n_files;
			}
		} else {
			__ERROR("Unknown option %s", argv[pos]);
		}
	}
	if (opt->n_files == 0) {
		free(opt);
		return NULL;
	}
	mkdir(opt->out_dir, 0755);
	return opt;
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
		} else if (!strcmp(argv[pos], "-l")) {
			opt->lib_type = get_library_index(argv[pos + 1]);
			if (opt->lib_type == -1)
				__ERROR("Unknown library %s", argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-1")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__ERROR("Inconsistent number of files");
			opt->n_files = n;
			opt->files_1 = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-2")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__ERROR("Inconsistent number of files");
			opt->n_files = n;
			opt->files_2 = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-sl")) {
			opt->split_len = atoi(argv[pos + 1]);
			pos += 2;
		} else if (argv[pos][0] != '-') {
			if (opt->n_files != 0)
				__ERROR("Unknown %s", argv[pos]);
			opt->files_1 = argv + pos;
			while (pos < argc && argv[pos][0] != '-') {
				++pos;
				++opt->n_files;
			}
		} else {
			__ERROR("Unknown option %s", argv[pos]);
		}
	}
	mkdir(opt->out_dir, 0755);
	return opt;
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

void assembly_opt_process(int argc, char *argv[])
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_assembly(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char log_dir[1024];
	snprintf(log_dir, 1024, "%s/assembly.log", opt->out_dir);
	init_log(log_dir);
	init_clock();
	print_info(argc, argv);
	assembly_process(opt);
}

void build_opt_process(int argc, char *argv[], void (*build_process)(struct opt_proc_t *))
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_build(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	snprintf(tmp_dir, 1024, "%s/assembly.log", opt->out_dir);
	init_log(tmp_dir);
	init_clock();
	print_info(argc, argv);
	build_process(opt);
}

void build_0_opt_process(int argc, char *argv[])
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_build0(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	snprintf(tmp_dir, 1024, "%s/build.log", opt->out_dir);
	init_log(tmp_dir);
	init_clock();
	print_info(argc, argv);
	build_0_process(opt);
}

void graph_query_opt_process(int argc, char *argv[])
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_build(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/query.log");
	init_log(tmp_dir);
	init_clock();
	graph_query_process(opt);
}

void graph_convert_opt_process(int argc, char *argv[])
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_build(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	snprintf(tmp_dir, 1024, "%s/convert.log", opt->out_dir);
	init_log(tmp_dir);
	init_clock();
	graph_convert_process(opt);
}

// ./skipping asm_graph0 -g <graph.bin> -o <output_folder>

int main(int argc, char *argv[])
{
	if (argc < 2) {
		print_usage(argv[0]);
		return -1;
	}
	if (!strcmp(argv[1], "assembly"))
		assembly_opt_process(argc, argv);
	else if (!strcmp(argv[1], "assembly2"))
		build_opt_process(argc, argv, &assembly2_process);
	else if (!strcmp(argv[1], "assembly_precount"))
		build_opt_process(argc, argv, &assembly_precount_process);
	else if (!strcmp(argv[1], "build_0"))
		build_0_opt_process(argc, argv);
	else if (!strcmp(argv[1], "build_barcode"))
		build_opt_process(argc, argv, &build_barcode_process);
	else if (!strcmp(argv[1], "build_barcode_fasta"))
		build_opt_process(argc, argv, &build_barcode_process_fasta);
	else if (!strcmp(argv[1], "build_0_1"))
		build_opt_process(argc, argv, &build_0_1_process);
	else if (!strcmp(argv[1], "build_1_2"))
		build_opt_process(argc, argv, &build_1_2_process);
	else if (!strcmp(argv[1], "build_2_3"))
		build_opt_process(argc, argv, &build_2_3_process);
	else if (!strcmp(argv[1], "build_3_4"))
		build_opt_process(argc, argv, &build_3_4_process);
	else if (!strcmp(argv[1], "bin2text"))
		graph_convert_opt_process(argc, argv);
	else if (!strcmp(argv[1], "query"))
		graph_query_opt_process(argc, argv);
	else
		print_usage(argv[0]);
	return 0;
}

