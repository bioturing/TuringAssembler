#include "scaffolding/scaffold.h"
#include "scaffolding/edge.h"
#include "scaffolding/output.h"
#include <stdlib.h>
#include <stdio.h>
#include "verbose.h"
#include "utils.h"
#include "scaffolding/global_params.h"

int greater_int(const void *a, const void *b) {
	int x = *(int *) a;
	int y = *(int *) b;
	if (x < y) {
		return 1;
	} else if (x == y) {
		return 0;
	} else {
		return -1;
	}
}

void print_scaffold(struct asm_graph_t *g, FILE *out_file, struct scaffold_type *scaffold)
{
	int total_length = 0, total_N = 0, largest_contig=0, N50 = 0, Gaps= 0, L50 = 0;
	int *list_length = calloc(scaffold->n_path, sizeof(int));
	for (int i = 0; i < scaffold->n_path; i++)  {
		struct scaffold_path *path = &scaffold->path[i];
		int sum_n = path->n_left_half + path->n_right_half;
		int *list_contig = calloc(sum_n, sizeof(int));
		for(int i = 0; i < path->n_left_half; i++)
			list_contig[path->n_left_half - 1 - i] = path->left_half[i];
		COPY_ARR(path->right_half, list_contig + path->n_left_half, path->n_right_half);
		struct contig_statictis sta = print_contig(g, out_file, i, sum_n, list_contig);
		total_length += sta.contig_length;
		total_N += sta.number_N;
		if (sta.contig_length > largest_contig )
			largest_contig = sta.contig_length;
		list_length[i] =  sta.contig_length;
		free(list_contig);
	}
	qsort(list_length, scaffold->n_path, sizeof(int), greater_int);
	int cur_sum_len = 0;
	for (int i = 0; i < scaffold->n_path; i++) {
		cur_sum_len += list_length[i];
		if (cur_sum_len >= total_length/2) {
			L50 = i+1;
			N50 = list_length[i];
			break;
		}
	}
	Gaps = total_N/ global_number_n;
	log_info("------- INFO SUMMARY ------ ");
	log_info("Total length: %d", total_length);
	log_info("Total N: %d", total_N);
	log_info("Largest contig: %d", largest_contig);
	log_info("N50 :%d", N50);
	log_info("L50 :%d", L50);
	log_info("");
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
		n_paths += (scaffold->path[i].n_left_half
			+ scaffold->path[i].n_right_half) > 1;
	char tmp_file[1024];
	snprintf(tmp_file, 1024, "%s/local_assembly_scaffold_path.txt", opt->out_dir);
	FILE *f = fopen(tmp_file, "w");
	log_debug("Number of scaffold paths: %d", n_paths);
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
	fclose(f);
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
			log_debug("reverse edge:%d", path->right_half[pos]);
		}
	} else {
		if (pos + 1 <= path->n_right_half)
			reverse_contig(g, &path->right_half[path->n_right_half-1-pos]);
		else {
			pos -= path->n_right_half;
			if (pos >= path->n_left_half)
				assert(0);
			reverse_contig(g, &path->left_half[pos]);
            log_debug("reverse edge:%d", path->left_half[pos]);
		}
	}
}

