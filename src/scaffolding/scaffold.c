#include "scaffolding/scaffold.h"
#include "scaffolding/edge.h"
#include "scaffolding/output.h"
#include <stdlib.h>
#include <stdio.h>
#include "verbose.h"
#include "scaffolding/global_params.h"

void print_scaffold(struct asm_graph_t *g, FILE *out_file, struct scaffold_type *scaffold)
{
	for (int i = 0; i < scaffold->n_path; i++) 
		print_contig(g, out_file, i, scaffold->path[i].n_contig, scaffold->path[i].list_i_contig);
}

void add_edge(struct scaffold_type *scaffold, struct scaffold_edge *edge)
{
	int v = scaffold->n_path-1;
	scaffold->path[v].n_contig += 1;
	scaffold->path[v].list_i_contig = realloc(scaffold->path[v].list_i_contig, 
			scaffold->path[v].n_contig * sizeof(int));
	scaffold->path[v].list_i_contig[scaffold->path[v].n_contig - 1] = edge->des;
}

void add_i_contig(struct scaffold_path *path, int i_contig)
{
	path->n_contig += 1;
	path->list_i_contig = realloc(path->list_i_contig, 
			path->n_contig* sizeof(int));
	path->list_i_contig[path->n_contig - 1] = i_contig;
}

struct scaffold_type *new_scaffold_type()
{
	struct scaffold_type *this = calloc(1, sizeof(struct scaffold_type));
	return this;
}

void add_path(struct scaffold_type *scaffold, struct scaffold_path *path)
{
	scaffold->n_path++;
	int n = scaffold->n_path;
	scaffold->path = realloc(scaffold->path, n*sizeof(struct scaffold_path));
	scaffold->path[n-1] = *path;
}

void destroy_path(struct scaffold_path *path)
{
	free(path->list_i_contig);
	free(path);
}

void deepcopy_path(struct scaffold_path *src, struct scaffold_path *des)
{
	if (des != NULL) 
		destroy_path(des);
	des = calloc(1, sizeof(struct scaffold_path));
	des->n_contig = src->n_contig;
	des->list_i_contig = calloc(des->n_contig, sizeof(int));
	for (int i = 0; i < des->n_contig; i++) 
		des->list_i_contig[i] = src->list_i_contig[i];
}

void print_scaffold_contig(struct scaffold_type *scaffold) 
{
	for (int i = 0 ; i < scaffold->n_path; i++) {
		VERBOSE_FLAG(1, "path\n");
		for (int j = 0; j < scaffold->path[i].n_contig; j++) {
			VERBOSE_FLAG(1, "contig %d ", scaffold->path[i].list_i_contig[j]);
		}
	}
}
