#ifndef __READ_LIST__
#define __READ_LIST__
#define MAX_READ_BUF 1000
#include "attribute.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "verbose.h"
struct read_list_t{
	int n;
	int m;
	struct read_t *reads;
};

void read_list_init(struct read_list_t *rlst);
void load_reads_from_fastq(char *r_path, struct read_list_t *rlst);
void load_read_pairs_from_fastq(char *r1_path, char *r2_path,
		struct read_list_t *list1, struct read_list_t *list2);
void push_read_to_list(struct read_t read, struct read_list_t *rlst);
#endif
