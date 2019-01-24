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
#include "kmhash.h"
#include "utils.h"
#include "verbose.h"

__KHASH_IMPL(kvert, kh_inline, kmkey_t, struct kvert_info_t, 1, kh_int64_hash_func, kh_int64_hash_equal)

void print_usage_assembly(const char *prog)
{
	__VERBOSE("Usage: %s assembly [options] -1 read_1.fq -2 read_2.fq\n", prog);
	__VERBOSE("Options: -t                     <number of threads>\n");
	__VERBOSE("         -s                     <pre-alloc size>\n");
	__VERBOSE("         -o                     <output directory>\n");
	__VERBOSE("         --kmer-small           <small kmer size>\n");
	__VERBOSE("         --kmer-large           <large kmer size>\n");
	__VERBOSE("         --filter-threshold     <kmer count cut off>\n");
	__VERBOSE("         --help                 show help\n");
}

void print_usage_count(const char *prog)
{
	__VERBOSE("Usage: %s count [options] read.[fq|fq]\n", prog);
	__VERBOSE("Options: -t                     <number of threads>\n");
	__VERBOSE("         -s                     <pre-alloc size>\n");
	__VERBOSE("         -o                     <output directory>\n");
	__VERBOSE("         --kmer                 <small kmer size>\n");
	__VERBOSE("         --filter-threshold     <kmer count cut off>\n");
}

void print_usage(const char *prog)
{
	__VERBOSE("Usage: %s\n", prog);
	__VERBOSE("          assembly [options] -1 read_1.fq -2 read_2.fq\n");
	__VERBOSE("          count [options] read.[fq|fa]\n");
}

struct opt_count_t *init_opt_count()
{
	struct opt_count_t *opt;
	opt = calloc(1, sizeof(struct opt_count_t));
	opt->n_threads = 1;
	opt->hash_size = (1 << 24);
	opt->kmer_master = 31;
	opt->kmer_slave = 17;
	opt->n_files = 0;
	opt->filter_thres = 0;
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
		} else if (!strcmp(argv[pos], "--kmer-small")) {
			opt->kmer_slave = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "--kmer-large")) {
			opt->kmer_master = atoi(argv[pos + 1]);
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
		} else if (!strcmp(argv[pos], "--filter-threshold")) {
			opt->filter_thres = atoi(argv[pos + 1]);
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

void assembly_opt_process(int argc, char *argv[])
{
	struct opt_count_t *opt;
	opt = parse_count_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage_assembly(argv[0]);
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/assembly.log");
	init_log(tmp_dir);

	int cmd_len = 0, i;
	for (i = 0; i < argc; ++i)
		cmd_len += strlen(argv[i]) + 1;
	char *cmd = malloc(cmd_len);
	cmd_len = 0;
	for (i = 0; i < argc; ++i)
		cmd_len += sprintf(cmd + cmd_len, i + 1 == argc ? "%s" : "%s ", argv[i]);
	__VERBOSE_LOG("INFO", "command: \"%s\"\n", cmd);
	free(cmd);

	__VERBOSE_LOG("INFO", "large kmer size: %d\n", opt->kmer_master);
	__VERBOSE_LOG("INFO", "small kmer size: %d\n", opt->kmer_slave);
	__VERBOSE_LOG("INFO", "pre-allocated hash table size: %d\n", opt->hash_size);
	__VERBOSE_LOG("INFO", "number of threads: %d\n", opt->n_threads);
	// __VERBOSE_LOG("INFO", "cut off with kmer count less or equal: %d\n", opt->filter_thres);
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

	int cmd_len = 0, i;
	for (i = 0; i < argc; ++i)
		cmd_len += strlen(argv[i]) + 1;
	char *cmd = malloc(cmd_len);
	cmd_len = 0;
	for (i = 0; i < argc; ++i)
		cmd_len += sprintf(cmd + cmd_len, i + 1 == argc ? "%s" : "%s ", argv[i]);
	__VERBOSE_LOG("INFO", "command: \"%s\"\n", cmd);
	free(cmd);

	__VERBOSE_LOG("INFO", "large kmer size: %d\n", opt->kmer_master);
	__VERBOSE_LOG("INFO", "small kmer size: %d\n", opt->kmer_slave);
	__VERBOSE_LOG("INFO", "pre-allocated hash table size: %d\n", opt->hash_size);
	__VERBOSE_LOG("INFO", "number of threads: %d\n", opt->n_threads);
	// __VERBOSE_LOG("INFO", "cut off with kmer count less or equal: %d\n", opt->filter_thres);
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
	kmer_test_process(opt);
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		print_usage(argv[0]);
		return -1;
	}

	if (!strcmp(argv[1], "assembly"))
		assembly_opt_process(argc, argv);
	else if (!strcmp(argv[1], "test"))
		test_opt_process(argc, argv);
	else
		print_usage(argv[0]);

	return 0;
}

