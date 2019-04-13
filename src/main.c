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
#include "k31_count.h"
#include "k63_count.h"
#include "process.h"
#include "utils.h"
#include "verbose.h"

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
}

void print_usage_count(const char *prog)
{
	__VERBOSE("Usage: %s count [options] read.[fq|fq]\n", prog);
	__VERBOSE("Options: -t                     <number of threads>\n");
	__VERBOSE("         -s                     <pre-alloc size>\n");
	__VERBOSE("         -o                     <output directory>\n");
	__VERBOSE("         --kmer                 <small kmer size>\n");
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
	__VERBOSE("Usage: %s\n", prog);
	__VERBOSE("          assembly [options] -1 read_1.fq -2 read_2.fq\n");
	__VERBOSE("          build_x_y [options]\n");
	__VERBOSE("          bin2text [options]\n");
	__VERBOSE("          build_barcode [options] -1 read_1.fq -2 read_2.fq\n");
	__VERBOSE("          build0 [options] -1 read_1.fq -2 read_2.fq\n");
	__VERBOSE("          query -i <input_graph> -f <list edge>\n");
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
			opt->in_path = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-f")) {
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

void print_opt_count_info(struct opt_count_t *opt, int argc, char *argv[])
{
	int cmd_len = 0, i;
	for (i = 0; i < argc; ++i)
		cmd_len += strlen(argv[i]) + 1;
	char *cmd = malloc(cmd_len);
	cmd_len = 0;
	for (i = 0; i < argc; ++i)
		cmd_len += sprintf(cmd + cmd_len, i + 1 == argc ? "%s" : "%s ", argv[i]);
	__VERBOSE_LOG("INFO", "command: \"%s\"\n", cmd);
	free(cmd);

	__VERBOSE_LOG("INFO", "k0: %d\n", opt->k0);
	__VERBOSE_LOG("INFO", "k1: %d\n", opt->k1);
	__VERBOSE_LOG("INFO", "k2: %d\n", opt->k2);
	__VERBOSE_LOG("INFO", "pre-allocated hash table size: %d\n", opt->hash_size);
	__VERBOSE_LOG("INFO", "number of threads: %d\n", opt->n_threads);
	__VERBOSE_LOG("INFO", "split bucket length: %d\n", opt->split_len);
	if (opt->n_files == 0) {
		__VERBOSE_LOG("INFO", "input: { stdin }\n");
	} else {
		if (opt->files_2 == NULL) {
			int len = 10, i;
			for (i = 0; i < opt->n_files; ++i)
				len += strlen(opt->files_1[i]) + 2;
			char *list_files = malloc(len);
			len = 0;
			len += sprintf(list_files, "{ ");
			for (i = 0; i < opt->n_files; ++i)
				len += sprintf(list_files + len,
						i + 1 == opt->n_files ? "%s" : "%s, ",
						opt->files_1[i]);
			sprintf(list_files + len, " }");
			__VERBOSE_LOG("INFO", "input: %s\n", list_files);
			free(list_files);
		} else {
			int len = 10, i;
			for (i = 0; i < opt->n_files; ++i)
				len += strlen(opt->files_1[i]) + strlen(opt->files_2[i]) + 6;
			char *list_files = malloc(len);
			len = 0;
			len += sprintf(list_files, "{ ");
			for (i = 0; i < opt->n_files; ++i)
				len += sprintf(list_files + len,
						i + 1 == opt->n_files ? "(%s, %s)" : "(%s, %s), ",
						opt->files_1[i], opt->files_2[i]);
			sprintf(list_files + len, " }");
			__VERBOSE_LOG("INFO", "input: %s\n", list_files);
			free(list_files);
		}
	}
}

void print_opt_build_info(struct opt_build_t *opt, int argc, char *argv[])
{
	int cmd_len = 0, i;
	for (i = 0; i < argc; ++i)
		cmd_len += strlen(argv[i]) + 1;
	char *cmd = malloc(cmd_len);
	cmd_len = 0;
	for (i = 0; i < argc; ++i)
		cmd_len += sprintf(cmd + cmd_len, i + 1 == argc ? "%s" : "%s ", argv[i]);
	__VERBOSE_LOG("INFO", "command: \"%s\"\n", cmd);
	free(cmd);
}

void assembly_opt_process63(int argc, char *argv[])
{
	struct opt_count_t *opt;
	opt = parse_count_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_assembly(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char log_dir[1024];
	strcpy(log_dir, opt->out_dir); strcat(log_dir, "/assembly.log");
	init_log(log_dir);
	print_opt_count_info(opt, argc, argv);
	// k63_process(opt);
}

void assembly_opt_process31(int argc, char *argv[])
{
	struct opt_count_t *opt;
	opt = parse_count_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_assembly(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char log_dir[1024];
	strcpy(log_dir, opt->out_dir); strcat(log_dir, "/assembly.log");
	init_log(log_dir);
	print_opt_count_info(opt, argc, argv);
	// k31_process(opt);
}

void assembly_opt_process(int argc, char *argv[])
{
	struct opt_count_t *opt;
	opt = parse_count_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_assembly(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char log_dir[1024];
	strcpy(log_dir, opt->out_dir); strcat(log_dir, "/assembly.log");
	init_log(log_dir);
	print_opt_count_info(opt, argc, argv);
	assembly_process(opt);
}

void test_opt_process(int argc, char *argv[])
{
	struct opt_count_t *opt;
	opt = parse_count_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_assembly(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/count.log");
	init_log(tmp_dir);
	print_opt_count_info(opt, argc, argv);
	// kmer_test_process(opt);
}

void test63_opt_process(int argc, char *argv[])
{
	struct opt_count_t *opt;
	opt = parse_count_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_assembly(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/count.log");
	init_log(tmp_dir);
	print_opt_count_info(opt, argc, argv);
	k63_test_process(opt);
}

void test31_opt_process(int argc, char *argv[])
{
	struct opt_count_t *opt;
	opt = parse_count_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_assembly(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/count.log");
	init_log(tmp_dir);
	print_opt_count_info(opt, argc, argv);
	k31_test_process(opt);
}

void build_opt_process(int argc, char *argv[], void (*build_process)(struct opt_build_t *))
{
	struct opt_build_t *opt;
	opt = parse_build_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_build(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/build.log");
	print_opt_build_info(opt, argc, argv);
	build_process(opt);
}

void clean_opt_process(int argc, char *argv[])
{
	struct opt_build_t *opt;
	opt = parse_build_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_build(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/build.log");
	init_log(tmp_dir);
	clean_process(opt);
}

void build0_opt_process(int argc, char *argv[])
{
	struct opt_count_t *opt;
	opt = parse_count_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_build0(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/build0.log");
	init_log(tmp_dir);
	print_opt_count_info(opt, argc, argv);
	build0_process(opt);
}

void graph_query_opt_process(int argc, char *argv[])
{
	struct opt_build_t *opt;
	opt = parse_build_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_build(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/query.log");
	init_log(tmp_dir);
	graph_query_process(opt);
}

void graph_convert_opt_process(int argc, char *argv[])
{
	struct opt_build_t *opt;
	opt = parse_build_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_build(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/convert.log");
	init_log(tmp_dir);
	graph_convert_process(opt);
}

// ./skipping asm_graph0 -g <graph.bin> -o <output_folder>

int main(int argc, char *argv[])
{
	if (argc < 2) {
		print_usage(argv[0]);
		return -1;
	}

	if (!strcmp(argv[1], "assembly31"))
		assembly_opt_process31(argc, argv);
	else if (!strcmp(argv[1], "assembly63"))
		assembly_opt_process63(argc, argv);
	else if (!strcmp(argv[1], "assembly"))
		assembly_opt_process(argc, argv);
	else if (!strcmp(argv[1], "test63"))
		test63_opt_process(argc, argv);
	else if (!strcmp(argv[1], "test31"))
		test31_opt_process(argc, argv);
	else if (!strcmp(argv[1], "build0"))
		build0_opt_process(argc, argv);
	else if (!strcmp(argv[1], "build_barcode"))
		// build_barcode_opt_process(argc, argv);
		build_opt_process(argc, argv, &build_barcode_process);
	else if (!strcmp(argv[1], "build_0_1"))
		// build0_1_opt_process(argc, argv);
		build_opt_process(argc, argv, &build0_1_process);
	else if (!strcmp(argv[1], "build_1_2"))
		// build1_2_opt_process(argc, argv);
		build_opt_process(argc, argv, &build1_2_process);
	else if (!strcmp(argv[1], "build_2_3"))
		// build2_3_opt_process(argc, argv);
		build_opt_process(argc, argv, &build2_3_process);
	else if (!strcmp(argv[1], "build_3_4"))
		// build3_4_opt_process(argc, argv);
		build_opt_process(argc, argv, &build3_4_process);
	else if (!strcmp(argv[1], "build_4_5"))
		// build4_5_opt_process(argc, argv);
		build_opt_process(argc, argv, &build4_5_process);
	else if (!strcmp(argv[1], "build_5_6"))
		// build5_6_opt_process(argc, argv);
		build_opt_process(argc, argv, &build5_6_process);
	else if (!strcmp(argv[1], "build_6_7"))
		build_opt_process(argc, argv, &build6_7_process);
	else if (!strcmp(argv[1], "bin2text"))
		graph_convert_opt_process(argc, argv);
	else if (!strcmp(argv[1], "query"))
		graph_query_opt_process(argc, argv);
	else if (!strcmp(argv[1], "clean"))
		clean_opt_process(argc, argv);
	else
		print_usage(argv[0]);

	return 0;
}

