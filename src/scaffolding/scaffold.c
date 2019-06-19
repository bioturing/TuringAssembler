#include "scaffolding/scaffold.h"
#include "scaffolding/edge.h"
#include "scaffolding/output.h"
#include <stdlib.h>
#include <stdio.h>

void print_scaffold(struct asm_graph_t *g, FILE *out_file, struct scaffold_type *scaffold)
{
	for (int i = 0; i < scaffold->n_connected_component; i++) 
		print_contig(g, out_file, i, scaffold->n_edge[i], scaffold->list_i_contig[i]);
}

void add_edge(struct scaffold_type *scaffold, struct scaffold_edge *edge)
{
	int v = scaffold->n_connected_component-1;
	scaffold->n_edge[v] += 1;
	scaffold->list_i_contig[v] = realloc(scaffold->list_i_contig[v], scaffold->n_edge[v] * sizeof(int));
	scaffold->list_i_contig[v][scaffold->n_edge[v] - 1] = edge->des;
}

void add_i_contig(struct scaffold_type *scaffold, int i_contig)
{
	int v = scaffold->n_connected_component-1;
	scaffold->n_edge[v] += 1;
	scaffold->list_i_contig[v] = realloc(scaffold->list_i_contig[v], scaffold->n_edge[v] * sizeof(int));
	scaffold->list_i_contig[v][scaffold->n_edge[v] - 1] = i_contig;
}

void add_connected_component(struct scaffold_type *scaffold, int start_contig)
{
	scaffold->n_connected_component++;
	int n = scaffold->n_connected_component;
	scaffold->n_edge = realloc(scaffold->n_edge, n*sizeof(int));
	scaffold->list_i_contig = realloc(scaffold->list_i_contig, 
			n * sizeof(int*));
	add_i_contig(scaffold, start_contig);
}

struct scaffold_type *new_scaffold_type()
{
	struct scaffold_type *this = calloc(1, sizeof(struct scaffold_type));
	return this;
}

