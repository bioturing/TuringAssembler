#include <stdint.h>
#include <assert.h>
#include "verbose.h"
void dfs_hamiltonian(int x, uint32_t depth, uint32_t *E, int *head, uint32_t *next, int *mark, uint32_t n_v, 
		uint32_t *res, uint32_t *best_res, uint32_t *best_n_res, uint32_t *listV)
{
	mark[x]--;
	mark[x^1]--;
	res[depth-1] = x;
	if (depth > *best_n_res) {
		*best_n_res = depth;
		for (uint32_t i = 0; i < depth; i++){
			best_res[i] = res[i];
			__VERBOSE("%d ", listV[res[i]]);
		}
		__VERBOSE("\n");
	} 
	for (int i = head[x]; i != -1; i = next[i]) {
		assert(mark[E[i]] >= 0);
		if (mark[E[i]] > 0) {
			dfs_hamiltonian(E[i], depth + 1, E, head, next, mark, n_v, res, best_res, best_n_res, listV);
		}
	}
	mark[x]++;
	mark[x^1]++;
}

