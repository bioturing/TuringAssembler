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
#include "process.h"
#include "time_utils.h"
#include "utils.h"
#include "verbose.h"
#include "log.h"

char *lib_str[] = {"sorted", "bioturing", "ust", "10x"};
int n_lib = 4;

int get_library_index(const char *str)
{
	int i;
	for (i = 0; i < n_lib; ++i)
		if (!strcmp(str, lib_str[i]))
			return i;
	return -1;
}

void print_usage()
{
	__VERBOSE("=====================================================================================================\n");
	__VERBOSE("TuringAssembler is a program developed by BioTuring for doing genome assembly with read-cloud technology\n");
	__VERBOSE("Please contact support@bioturing.com if you need further support.\n");
	__VERBOSE("=====================================================================================================\n");
	__VERBOSE("TuringAssembler - A genome assembler for read-cloud technology\n");
	__VERBOSE("Version: %s%s\n", VERSION_STRING, GIT_SHA);
	__VERBOSE("Usage:\n");
	__VERBOSE("          assembly3 [options] -1 read_1.fq -2 read_2.fq -l ust/bioturing/sorted\n");
	__VERBOSE("          local_assembly [options] -i graph.bin -1 R1_sorted.fq -2 R2_sorted.fq -l sorted -lk 31 -lc scaffold.full.fasta\n");
	__VERBOSE("Required parameters:\n");
	__VERBOSE("          -1            Forward reads. e.g: R1.fq or R1_lane1.fq R1_lane2.fq ...\n");
	__VERBOSE("          -2            Reverse reads. e.g: R2.fq or R2_lane1.fq R2_lane2.fq ...\n");
	__VERBOSE("          -l            Type of reads library:\n");
	__VERBOSE("                           bioturing: has BX:Z:<barcode> in the comment of the read names\n");
	__VERBOSE("                           sorted   : Sorted reads produced by TuringAssembler. Must be accompanied with -I barcode.idx\n");
	__VERBOSE("                           ust      : reads are generated by TELL-Seq protocols. Must be accompanied with -I I1.fq (I1_lane1.fq I1_lane2.fq)\n");
	__VERBOSE("Optional parameters:\n");
	__VERBOSE("          -t            Number of threads [4]\n");
	__VERBOSE("          -k0           Kmer size for global assembly process [45]\n");
	__VERBOSE("          -lk           Kmer size for local assembly [31]\n");
	__VERBOSE("          -lc           Output file after local assembly step [scaffold.full.fasta]\n");
	__VERBOSE("          -metagenomics Doing assembly for metagenomics dataset [no]\n");
	__VERBOSE("          -o            Output directory [./]\n");
	__VERBOSE("          -i            Input graph binary file. Only use for sub-processes [./graph_k_xx_level_x.bin]\n");
	__VERBOSE("          -sm           Maximum memory size for read sorting (GB) [32]\n");
	__VERBOSE("          -v            Verbose mode. Print log trace and log debug [no]\n");
}

void check_process_assembly3(int argc, struct opt_proc_t *opt)
{
	int stop = 0;
	__VERBOSE("=====================================================================================================\n");
	__VERBOSE("TuringAssembler is a program developed by BioTuring for doing genome assembly with read-cloud technology\n");
	__VERBOSE("Please contact support@bioturing.com if you need further support.\n");
	__VERBOSE("=====================================================================================================\n");
	__VERBOSE("TuringAssembler - A genome assembler for read-cloud technology\n");
	if (argc == 2) {
		__VERBOSE("Version: %s%s\n", VERSION_STRING, GIT_SHA);
		__VERBOSE("TuringAssembler assembly3 usage:\n");
	}
	if (opt->n_files <= 0 || opt->lib_type == -1) {
		stop = 1;
		__VERBOSE("Missing required argument(s) for assembly3:\n");
	}

	if (opt->n_files <= 0) {
		__VERBOSE("          -1            Forward reads. e.g: R1.fq or R1_lane1.fq R1_lane2.fq ...\n");
		__VERBOSE("          -2            Reverse reads. e.g: R2.fq or R2_lane1.fq R2_lane2.fq ...\n");
	}
	if (opt->lib_type == -1) {
		__VERBOSE("          -l            Type of reads library:\n");
		__VERBOSE("                           bioturing: has BX:Z:<barcode> in the comment of the read names\n");
		__VERBOSE("                           sorted   : Sorted reads produced by TuringAssembler. Must be accompanied with -I barcode.idx\n");
		__VERBOSE("                           ust      : reads are generated by TELL-Seq protocol. Must be accompanied with -I I1.fq (I1_lane1.fq I1_lane2.fq)\n");
	}
	if (stop)
		exit(1);

}

void check_process_local_ass(int argc, struct opt_proc_t *opt)
{
	int stop = 0;
	__VERBOSE("=====================================================================================================\n");
	__VERBOSE("TuringAssembler is a program developed by BioTuring for doing genome assembly with read-cloud technology\n");
	__VERBOSE("Please contact support@bioturing.com if you need further support.\n");
	__VERBOSE("=====================================================================================================\n");
	__VERBOSE("TuringAssembler - A genome assembler for read-cloud technology\n");
	__VERBOSE("Suggested arguments for local_assembly:\n");
	__VERBOSE("          -lk           Kmer size for local assembly\n");
	__VERBOSE("                        For read length 100-148bp: 69/71/73/75\n");
	__VERBOSE("                        For read length 71-100bp : 27/29/31/33\n");
	__VERBOSE("          -sm           Maximum memory size for read sorting (GB) [32]\n");
	if (argc == 2) {
		__VERBOSE("Version: %s%s\n", VERSION_STRING, GIT_SHA);
		__VERBOSE("TuringAssembler local_assembly usage:\n");
	}
	if (opt->n_files <= 0 || opt->lib_type == -1 || opt->in_file == NULL) {
		stop = 1;
		__VERBOSE("Missing required argument(s) for local_assembly:\n");
	}

	if (opt->n_files <= 0) {
		__VERBOSE("          -1            Forward reads. e.g: R1.fq or R1_lane1.fq R1_lane2.fq ...\n");
		__VERBOSE("          -2            Reverse reads. e.g: R2.fq or R2_lane1.fq R2_lane2.fq ...\n");
	}
	if (opt->lib_type == -1) {
		__VERBOSE("          -l            Type of reads library:\n");
		__VERBOSE("                           bioturing: has BX:Z:<barcode> in the comment of the read names\n");
		__VERBOSE("                           sorted   : Sorted reads produced by TuringAssembler. Must be accompanied with -I barcode.idx\n");
		__VERBOSE("                           ust      : reads are generated by TELL-Seq protocol. Must be accompanied with -I I1.fq (I1_lane1.fq I1_lane2.fq)\n");
	}
	if (opt->in_file == NULL) {
		__VERBOSE("          -i            Input graph binary file. Only use for sub-processes [./graph_k_xx_added_barcode.bin]\n");
	}
	if (stop)
		exit(1);

}

struct opt_proc_t *init_opt_proc()
{
	struct opt_proc_t *opt;
	opt = calloc(1, sizeof(struct opt_proc_t));
	opt->log_level = LOG_INFO;
	opt->n_threads = 4;
	opt->hash_size = 1 << 24;
	opt->k0 = 45; /* Default kmer size */
	opt->n_files = 0;
	opt->split_len = 1000;
	opt->files_1 = opt->files_2 = NULL;
	opt->in_file = NULL;
	opt->in_fasta = NULL;
	opt->in_fastg = NULL;
	opt->files_1 = NULL;
	opt->files_2 = NULL;
	opt->files_I = NULL;
	opt->out_dir = "."; /* Default output folder */
	opt->lib_type = -1;/* Default read library */
	opt->mmem = 32; /*Default memory limit for reads sorting */
	opt->lk = 31; /* Default local kmer size */
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

void check_proc_opt(struct opt_proc_t *opt)
{
	if (opt->k0 < 17)
		log_error("Kmer size (k0) must be greater than or equal to 17");
	if (opt->lk < 17)
		log_error("Local kmer size (lk) must be greater than or equal to 17");
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
		} else if (!strcmp(argv[pos], "-v")) {
			opt->log_level = LOG_DEBUG;
			pos += 1;
		} else if (!strcmp(argv[pos], "-vv")) {
			opt->log_level = LOG_DEBUG_TECH;
			pos += 1;
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
		} else if (!strcmp(argv[pos], "-bx")) {
			opt->bx_str = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-f")) {
			opt->in_fasta = argv[pos + 1];
			pos += 2;
                } else if (!strcmp(argv[pos], "-fg")) {
                        opt->in_fastg = argv[pos + 1];
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
				__ERROR("Unknown library %s", argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-sm")) {
			opt->mmem = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-1")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__ERROR("Inconsistent number of files");
			opt->n_files = n;
			opt->files_1 = argv + pos + 1;
			pos += (n + 1);
		}else if (!strcmp(argv[pos], "-var")) {
			n = opt_count_list(argc - pos, argv + pos);
			opt->var = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-2")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__ERROR("Inconsistent number of files");
			opt->n_files = n;
			opt->files_2 = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-I")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__ERROR("Inconsistent number of files");
			opt->n_files = n;
			opt->files_I = argv + pos + 1;
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
		} else if (!strcmp(argv[pos], "-lc")){
			opt->lc = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-lk")){
			opt->lk = atoi(argv[pos + 1]);
			pos += 2;
		} else {
			__ERROR("Unknown option %s", argv[pos]);
		}
	}
	check_proc_opt(opt);
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

void build_opt_process(int argc, char *argv[], void (*build_process)(struct opt_proc_t *))
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (build_process == &(assembly3_process)) {
		check_process_assembly3(argc, opt);
	} else if (build_process == &(build_bridge_process)){
		check_process_local_ass(argc, opt);
	};
	if (!opt->k0 & 1) {
		__ERROR("Kmer size must be an odd number! e.g -k0 45\n");
	}
	char tmp_dir[1024];
	snprintf(tmp_dir, 1024, "%s/assembly.log", opt->out_dir);
	init_log(tmp_dir);
	init_clock();
	print_info(argc, argv);
	build_process(opt);
	free(opt);
	//destroy_opt_proc_t(opt); /* No need */
}

void build_0_opt_process(int argc, char *argv[])
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage();
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	snprintf(tmp_dir, 1024, "%s/build.log", opt->out_dir);
	init_log(tmp_dir);
	init_clock();
	print_info(argc, argv);
	build_0_process(opt);
}

//void graph_query_opt_process(int argc, char *argv[])
//{
//	struct opt_proc_t *opt;
//	opt = parse_proc_option(argc - 2, argv + 2);
//	if (opt == NULL) {
//		print_usage();
//		__ERROR("Error parsing arguments");
//	}
//	char tmp_dir[1024];
//	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/query.log");
//	init_log(tmp_dir);
//	init_clock();
//	graph_query_process(opt);
//}

void build_bridge_opt_process(int argc, char *argv[])
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (opt == NULL){
		print_usage();
		__ERROR("Error parsing arguments");
	}
	char path[1024];
	sprintf(path, "%s/build_bridge.log", opt->out_dir);
	init_logger(LOG_DEBUG_TECH, path);
	set_log_stage("Local assembly");
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/build_bridge.log");
	init_log(tmp_dir);
	init_clock();
	build_bridge_process(opt);
	free(opt);
}

void reduce_read_opt_process(int argc, char *argv[])
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (opt == NULL){
		print_usage();
		__ERROR("Error parsing arguments");
	}
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/reduce_read.log");
	init_log(tmp_dir);
	init_clock();
	reduce_read_process(opt);
	free(opt);
}

void resolve_local_opt_process(int argc, char *argv[])
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (opt == NULL){
		print_usage();
		log_error("Error parsing arguments");
	}
	resolve_local_process(opt);
	free(opt);
}

void graph_convert_opt_process(int argc, char *argv[])
{
	struct opt_proc_t *opt;
	opt = parse_proc_option(argc - 2, argv + 2);
	if (opt == NULL) {
		print_usage();
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
		print_usage();
		return -1;
	}
	if (!strcmp(argv[1], "assembly"))
		build_opt_process(argc, argv, &assembly_process);
	else if (!strcmp(argv[1], "assembly3"))
		build_opt_process(argc, argv, &assembly3_process);
	else if (!strcmp(argv[1], "build_0"))
		build_0_opt_process(argc, argv);
	else if (!strcmp(argv[1], "build_barcode"))
		build_opt_process(argc, argv, &build_barcode_info);
	else if (!strcmp(argv[1], "build_barcode_scaffold"))
		build_opt_process(argc, argv, &build_barcode_scaffold);
	else if (!strcmp(argv[1], "build_barcode_fasta"))
		build_opt_process(argc, argv, &build_barcode_process_fasta);
	else if (!strcmp(argv[1], "build_barcode_fastg"))
		build_opt_process(argc, argv, &build_barcode_process_fastg);
	else if (!strcmp(argv[1], "build_0_1"))
		build_opt_process(argc, argv, &build_0_1_process);
	else if (!strcmp(argv[1], "build_1_2"))
		build_opt_process(argc, argv, &build_1_2_process);
	else if (!strcmp(argv[1], "build_2_3"))
		build_opt_process(argc, argv, &build_2_3_process);
	else if (!strcmp(argv[1], "build_3_4"))
		build_opt_process(argc, argv, &build_3_4_process);
	else if (!strcmp(argv[1], "build_3_4_nobc"))
		build_opt_process(argc, argv, &build_3_4_no_bc_rebuild_process);
	else if (!strcmp(argv[1], "build_4_5"))
		build_opt_process(argc, argv, &build_4_5_process);
	else if (!strcmp(argv[1], "bin2text"))
		graph_convert_opt_process(argc, argv);
	/*else if (!strcmp(argv[1], "query"))
		graph_query_opt_process(argc, argv);*/
	else if (!strcmp(argv[1], "build_bridge"))
		build_bridge_opt_process(argc, argv);
	else if (!strcmp(argv[1], "local_assembly"))
		build_opt_process(argc, argv, &build_bridge_process);
	else if (!strcmp(argv[1], "resolve_bulges"))
		build_opt_process(argc, argv, &resolve_bulges_process);
	else if (!strcmp(argv[1], "resolve_complex_bulges"))
		build_opt_process(argc, argv, &resolve_complex_bulges_process);
	else if (!strcmp(argv[1], "reduce_reads"))
		reduce_read_opt_process(argc, argv);
	else if (!strcmp(argv[1], "resolve_local"))
		resolve_local_opt_process(argc, argv);
	else if (!strcmp(argv[1], "build_scaffolding_1_2"))
		build_opt_process(argc, argv, &build_scaffolding_1_2_process); 
	else if (!strcmp(argv[1], "build_scaffolding_test"))
		build_opt_process(argc, argv, &build_scaffolding_test_process);
	else if(!strcmp(argv[1], "resolve_1_2"))
		build_opt_process(argc, argv, &resolve_1_2_process);
	else if (!strcmp(argv[1], "resolve_complex_bulges"))
		build_opt_process(argc, argv, &resolve_complex_bulges_process);
	else if (!strcmp(argv[1], "mm_index"))
		build_opt_process(argc, argv, &index_mm_process);
	else if (!strcmp(argv[1], "barcode_hit"))
		build_opt_process(argc, argv, &hits_barcode_process);
	else if (!strcmp(argv[1], "count_bx"))
		build_opt_process(argc, argv, &count_bx_process);
	else
		print_usage();
	return 0;
}

