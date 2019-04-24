#include <stdint.h>
#include <assert.h>
#include "verbose.h"
#include <stdlib.h>
void dfs_hamiltonian(int x, int depth, int n_adj, int *list_adj, int n_v, float *E, int *remain_unvisited, 
		int *cur_path, int *best_n_res, int *best_res, int *listV)
{
	int cmp(const void *i, const void *j)
	{
		float zzz = E[x * n_v + *(int*)j] - E[x *n_v + *(int*)i];
		if (zzz < 0) 
			return -1;
		if (zzz == 0)
			return 0;
		return 1;
	}
	__VERBOSE("dfs %d %d\n", x, depth);
	remain_unvisited[x]--;
	remain_unvisited[x^1]--;
	assert(remain_unvisited[x]>=0 && remain_unvisited[x^1]>=0);
	cur_path[depth-1] = x;
	if (depth > *best_n_res) {
		*best_n_res = depth;
		__VERBOSE("new path depth %d\n", depth);
		for (int i = 0; i < depth; i++){
			best_res[i] = cur_path[i];
			__VERBOSE("%d ", cur_path[i]);
		}
		__VERBOSE("\n");
	} 

	int *list_sort = calloc(n_adj, sizeof(int));
	for (int i = 0; i < n_adj; i++) {
		list_sort[i] = list_adj[i];
	}
	qsort(list_sort, n_adj, sizeof(int), cmp);
	for (int i = 0; i < n_adj; i++) {
		int adj = list_sort[i];
		if (remain_unvisited[adj] && E[x * n_v + adj]) {
			dfs_hamiltonian(adj, depth+1, n_adj, list_adj, n_v, E, remain_unvisited, 
				cur_path, best_n_res, best_res, listV);
		}
	}
	remain_unvisited[x]++;
	remain_unvisited[x^1]++;
	free(list_sort);
}

