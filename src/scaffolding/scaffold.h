#ifndef __SCAFFOLDING_SCAFFOLD_H__
#define __SCAFFOLDING_SCAFFOLD_H__
#include "scaffolding/edge.h"
#include <stdio.h>

struct scaffold_path {
	int n_right_half, n_left_half;
	int *right_half, *left_half;
};

struct scaffold_type {
	int n_path;
	struct scaffold_path *path;
};
void append_i_contig(struct scaffold_path *path, int i_contig);
void prepend_i_contig(struct scaffold_path *path, int i_contig);

struct scaffold_type *new_scaffold_type();
void print_scaffold(struct asm_graph_t *g, FILE *out_file, struct scaffold_type *scaffold);
void add_path(struct scaffold_type *scaffold, struct scaffold_path *path);
void print_scaffold_contig(struct opt_proc_t *opt, struct scaffold_type *scaffold);
int get_last_n(struct scaffold_path *path, int is_left, int pos);
void reverse_n_th(struct asm_graph_t *g, struct scaffold_path *path, int is_left, int pos);
#endif /* __SCAFFOLDING_SCAFFOLD_H__*/
