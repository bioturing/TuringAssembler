#ifndef SCAFFOLDING_ALGORITHM_H
#define SCAFFOLDING_ALGORITHM_H

#include "assembly_graph.h"

void dfs_hamiltonian(int x, int depth, int n_adj, int *list_adj, int n_v, float *E, int *remain_unvisited, int *cur_path, int *best_n_res, int *best_res, int *listV);
void algo_find_hamiltonian(FILE *out_file, struct asm_graph_t *g, float *E, int n_v, int *res, int *n_res, int *listV, float avg_bin_hash);
int get_new_size(int size);

#endif
