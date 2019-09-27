#include "scaffolding/scaffold.h"
#include "scaffolding/edge.h"
#include "scaffolding/output.h"
#include <stdlib.h>
#include <stdio.h>
#include "verbose.h"
#include "utils.h"
#include "scaffolding/global_params.h"

void print_scaffold(struct asm_graph_t *g, FILE *out_file, struct scaffold_type *scaffold)
{
	for (int i = 0; i < scaffold->n_path; i++)  {
		struct scaffold_path *path = &scaffold->path[i];
		int sum_n = path->n_left_half + path->n_right_half;
		int *list_contig = calloc(sum_n, sizeof(int));
		for(int i = 0; i < path->n_left_half; i++)
			list_contig[path->n_left_half - 1 - i] = path->left_half[i];
		COPY_ARR(path->right_half, list_contig + path->n_left_half, path->n_right_half);
		print_contig(g, out_file, i, sum_n, list_contig);
		free(list_contig);
	}
}

void append_i_contig(struct scaffold_path *path, int i_contig)
{
	path->n_right_half += 1;
	path->right_half = realloc(path->right_half, path->n_right_half *sizeof(int));
	path->right_half[path->n_right_half - 1] = i_contig;
}

void prepend_i_contig(struct scaffold_path *path, int i_contig)
{
	path->n_left_half += 1;
	path->left_half = realloc(path->left_half, path->n_left_half *sizeof(int));
	path->left_half[path->n_left_half - 1] = i_contig;
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
	scaffold->path[n-1].n_right_half = path->n_right_half;
	scaffold->path[n-1].n_left_half = path->n_left_half;
	scaffold->path[n-1].right_half = path->right_half;
	scaffold->path[n-1].left_half = path->left_half;
}

void destroy_path(struct scaffold_path *path)
{
	free(path->left_half);
	free(path->right_half);
	free(path);
}

void print_scaffold_contig(struct opt_proc_t *opt, struct scaffold_type *scaffold) 
{
	int n_paths = 0;
	for (int i = 0; i < scaffold->n_path; ++i)
		n_paths += scaffold->path[i].n_left_half
			+ scaffold->path[i].n_right_half > 1;
	char tmp_dir[1024];
	snprintf(tmp_dir, 1024, "%s/local_assembly_scaffold_path.txt", opt->out_dir);
	FILE *f = fopen(tmp_dir, "w");
	fprintf(f, "%d\n", n_paths);
	for (int i = 0 ; i < scaffold->n_path; i++) {
		struct scaffold_path *path = &scaffold->path[i];
		int sum_n = path->n_left_half + path->n_right_half;
		if (sum_n <= 1)
			continue;
		int *list_contig = calloc(sum_n, sizeof(int));
		for(int i = 0; i < path->n_left_half; i++)
			list_contig[path->n_left_half - 1 - i] = path->left_half[i];
		COPY_ARR(path->right_half, list_contig + path->n_left_half,
				path->n_right_half);
		fprintf(f, "%d\n", sum_n);
		for (int j = 0; j < sum_n; ++j)
			fprintf(f, "%d ", list_contig[j]);
		fprintf(f, "\n");
		free(list_contig);
	}
}

int get_last_n(struct scaffold_path *path, int is_left, int pos)
{
	if (is_left) {
		if (pos + 1 <= path->n_left_half)
			return path->left_half[path->n_left_half-1-pos];
		pos -= path->n_left_half;
		if (pos >= path->n_right_half)
			return -1;
		return path->right_half[pos];
	} else {
		if (pos + 1 <= path->n_right_half)
			return path->right_half[path->n_right_half-1-pos];
		pos -= path->n_right_half;
		if (pos >= path->n_left_half)
			return -1;
		return path->left_half[pos];
	}
}

void reverse_contig(struct asm_graph_t *g, int *i_contig)
{
	*i_contig = g->edges[*i_contig].rc_id;
}

void reverse_n_th(struct asm_graph_t *g, struct scaffold_path *path, int is_left, int pos)
{
	if (is_left) {
		if (pos + 1 <= path->n_left_half)
			reverse_contig(g, &path->left_half[path->n_left_half-1-pos]);
		else {
			pos -= path->n_left_half;
			if (pos >= path->n_right_half)
				assert(0);
			int x = path->right_half[pos];
			reverse_contig(g, &path->right_half[pos]);
			log_debug("reverse edge:%d\n", path->right_half[pos]);
		}
	} else {
		if (pos + 1 <= path->n_right_half)
			reverse_contig(g, &path->right_half[path->n_right_half-1-pos]);
		else {
			pos -= path->n_right_half;
			if (pos >= path->n_left_half)
				assert(0);
			reverse_contig(g, &path->left_half[pos]);
            log_debug("reverse edge:%d\n", path->left_half[pos]);
		}
	}
}

