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

#define KMCP			"KMCP"
#define KMCS			"KMCS"

#define LIB_TYPE_SORTED		0
#define LIB_TYPE_BIOT		1
#define LIB_TYPE_UST		2
#define LIB_TYPE_10X		3

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
	char *out_dir;
	char *in_file;
	char *in_fasta;
	int mmem;
};

#endif
