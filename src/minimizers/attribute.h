#ifndef _ATTRIBUTE_H_
#define _ATTRIBUTE_H_

#include <pthread.h>
#include <stdint.h>

#include "pthread_barrier.h"

#define SIZE_1MB		1048576
#define SIZE_2MB		2097152
#define SIZE_4MB		4194304
#define SIZE_16MB		16777216
#define SIZE_128MB		134217728
#define BARCODE_LEN_10X 16
#define UMI_LEN_10X 7

#define KMCP			"KMCP"
#define KMCS			"KMCS"

#define MINIMIZERS_WINDOW       21
#define MINIMIZERS_KMER       21

#define LIB_TYPE_SORTED		0
#define LIB_TYPE_BIOT		1
#define LIB_TYPE_UST		2
#define LIB_TYPE_10X		3

#define MAX_PATH 4096

#if !defined(GIT_SHA)
#define GIT_SHA			"unknown"
#endif

#define VERSION_STRING		"0.9-"

#define EPS			1e-6

typedef int64_t gint_t;

struct read_t {
	char *seq;
	char *name;
	char *qual;
	char *note;
	char *info;
	int len;
};

struct opt_proc_t {
	int n_threads;
	int hash_size;
	int k0;
	int k1;
	int k2;
	int split_len;
	int lib_type;
	int n_files;
	char **files_1, **files_2, **files_I;
	int metagenomics;
	char *out_dir;
	char *in_file;
	char *in_fasta;
	char *in_fastg;
	char *in_contig_file;
	int mmem;
	char *lc;
	int lk;
	char *bx_str;
};

struct read_path_t {
	char *R1_path;
	char *R2_path;
	char *idx_path;
};

#endif
