#ifndef __SCAFFOLDING_SCAFFOLD_H__
#define __SCAFFOLDING_SCAFFOLD_H__
#include "scaffolding/edge.h"
#include <stdio.h>
struct scaffold_type {
	int n_connected_component;
	int *n_edge;
	int **list_i_contig;
};
struct scaffold_type *new_scaffold_type();
void print_scaffold(struct asm_graph_t *g, FILE *out_file, struct scaffold_type *scaffold);
void add_edge(struct scaffold_type *scaffold, struct scaffold_edge *edge);
void add_connected_component(struct scaffold_type *scaffold, int start_contig);
#endif /* __SCAFFOLDING_SCAFFOLD_H__*/
