#include <stdint.h>
#include <assert.h>
#include "verbose.h"
void dfs_hamiltonian(int x, int depth, int n_adj, int *list_adj, int n_v, int *E, int *remain_unvisited, 
		int *cur_path, int *best_n_res, int *best_res, int *listV)
{
	int a = best_res;
	__VERBOSE("dfs %d %d %d\n", x, depth, a);
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
	for (int i = 0; i < n_adj; i++) {
		int adj = list_adj[i];
		if (remain_unvisited[adj] && E[x * n_v + adj]) {
			dfs_hamiltonian(adj, depth+1, n_adj, list_adj, n_v, E, remain_unvisited, 
				cur_path, best_n_res, best_res, listV);
		}
	}
	remain_unvisited[x]++;
	remain_unvisited[x^1]++;
}

