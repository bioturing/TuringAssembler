#include <stdint.h>
#include <assert.h>
#include "verbose.h"
#include <stdlib.h>
#include <math.h>
#include "assembly_graph.h"
#include "scaffolding/output.h"
#include "scaffolding/global_params.h"
#include "scaffolding/edge.h"

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
// #########################      BEGIN FUNCTION      ########################
	remain_unvisited[x]--;
	remain_unvisited[x^1]--;
	assert(remain_unvisited[x]>=0 && remain_unvisited[x^1]>=0);
	cur_path[depth-1] = x;
	if (depth > *best_n_res) {
		*best_n_res = depth;
		for (int i = 0; i < depth; i++){
			best_res[i] = cur_path[i];
		}
	} 

	int *list_sort = calloc(n_adj, sizeof(int));
	for (int i = 0; i < n_adj; i++) {
		list_sort[i] = list_adj[i];
	}
	qsort(list_sort, n_adj, sizeof(int), cmp);
	for (int i = 0; i < n_adj; i++) {
		int adj = list_sort[i];
		if (remain_unvisited[adj] && E[x * n_v + adj] > 0) {
			dfs_hamiltonian(adj, depth+1, n_adj, list_adj, n_v, E, remain_unvisited, 
				cur_path, best_n_res, best_res, listV);
		}
	}
	remain_unvisited[x]++;
	remain_unvisited[x^1]++;
	free(list_sort);
}

void algo_find_hamiltonian(FILE *out_file, struct asm_graph_t *g, float *E, int n_v, int *res, int *n_res, int *listV, float avg_bin_hash)
{
	// deprecated due to low performance
	void insert_short_contigs(int n_big_contigs, int *big_contigs, int n_insert, int *arr_insert, 
			int n_short, int *arr_i_short, int *mark_short)
	{
		//todo @huu insert to the beginning
		VERBOSE_FLAG(1, "big contig length %d:", n_big_contigs);
		for (int i = 0; i < n_big_contigs; i++){
			VERBOSE_FLAG(3, "%d ", big_contigs[i]);
		}
		for (int i = 1; i < n_insert - 1; i++) {
			VERBOSE_FLAG(3, "i insert: %d\n", i);
			int max_score = -1, pos = -1;
			for (int j = 0; j < n_short; j++) if (mark_short[j] == 0) {
				int score0 = get_score_big_small(big_contigs[i-1], arr_i_short[j], g, avg_bin_hash);
				int score1 = get_score_big_small(g->edges[big_contigs[i]].rc_id, g->edges[arr_i_short[j]].rc_id, g, avg_bin_hash);
				int score = score0 + score1;
				if (score > max_score){
					max_score = score;
					pos = j;
				}
			}
			VERBOSE_FLAG(3, "max score: %d \n",max_score);
			//todo @huu check if normal score large enough
			if (max_score > 15) {
				arr_insert[i] = arr_i_short[pos];
				mark_short[pos] = 1;
			}
		}
		//todo @huu insert to the end
	}

	void add_insert_except_m1(int *n_arr, int **new_arr, int ele)
	{
		if (ele == -1)
			return;
		*new_arr = realloc(*new_arr, (*n_arr+1) * sizeof(int));
		(*new_arr)[*n_arr] = ele;
		++(*n_arr);
	}

	void merge_big_and_small(int *best_n_res, int **best_res, int n_insert, int *arr_insert)
	{
		int *new_arr = NULL;
		int n_arr = 0;
		add_insert_except_m1(&n_arr, &new_arr, arr_insert[0]);
		for (int i = 0; i < *best_n_res; i++) {
			add_insert_except_m1(&n_arr, &new_arr, (*best_res)[i]); 
			add_insert_except_m1(&n_arr, &new_arr, arr_insert[i]); 
		}
		free(*best_res);
		*best_n_res = n_arr;
		*best_res = new_arr;
		VERBOSE_FLAG(3, "after merge contig length %d:", *best_n_res);
		for (int i = 0; i < *best_n_res; i++){
			VERBOSE_FLAG(3, "%d ", (*best_res)[i]);
		}
	}

	void find_longest_path_from_node(int x, int *best_n_hamiltonian_path, int *best_hamiltonian_path, int *remain_unvisited)
	{
		VERBOSE_FLAG(1, "find longest path from node");
		remain_unvisited[x]--;
		remain_unvisited[x^1]--;
		assert(*best_n_hamiltonian_path == 0);
		*best_n_hamiltonian_path = 1;
		best_hamiltonian_path[0] = x;
		while (1) {
			int *list_adj = NULL, count_adj = 0;
			int last_pos = best_hamiltonian_path[*best_n_hamiltonian_path-1];
			for (int i = 0; i < n_v; i++) if (E[last_pos * n_v + i] != 0 && remain_unvisited[i]){
				assert(remain_unvisited[i] >0);
				list_adj = realloc(list_adj, (count_adj+1) * sizeof(int));
				list_adj[count_adj] = i;
				count_adj++;
			}
			if (count_adj == 0)
				break;
			VERBOSE_FLAG(3, "node %d count_adj %d\n", listV[last_pos], count_adj);
			int best_n_local_path = 0 , *best_local_path = calloc(n_v, sizeof(int)), best_add_len = 0;
			for (int i_adj = 0; i_adj < count_adj; i_adj++) {
				int adj = list_adj[i_adj];
				VERBOSE_FLAG(3, "dfs from %d %d %d\n", i_adj, adj, remain_unvisited[adj]);
				int *cur_path = calloc(n_v, sizeof(int));
				dfs_hamiltonian(
					adj, 1, count_adj, list_adj, n_v, E, remain_unvisited,
					cur_path, &best_add_len, best_hamiltonian_path + *best_n_hamiltonian_path, listV
				);
				free(cur_path);
			}
			VERBOSE_FLAG(3, "best n from node n_v best add len %d %d %d\n", *best_n_hamiltonian_path, n_v, best_add_len);
			for (int i = 0; i < *best_n_hamiltonian_path + best_add_len ; i++) {
				VERBOSE_FLAG(3, "%d ", best_hamiltonian_path[i]);
			}
			for (int i_path = *best_n_hamiltonian_path; i_path < *best_n_hamiltonian_path + best_add_len; i_path++){
				remain_unvisited[best_hamiltonian_path[i_path]]--;
				remain_unvisited[best_hamiltonian_path[i_path]^1]--;
			}
			*best_n_hamiltonian_path += best_add_len;
			free(best_local_path);
			free(list_adj);
		}
		for (int i = 0; i < *best_n_hamiltonian_path; i++){
			remain_unvisited[best_hamiltonian_path[i]]++;
			remain_unvisited[best_hamiltonian_path[i]^1]++;
		}
	}

	void dfs_find_connected_component(int root, int *remain_unvisited, int
			*connected_component)
	{
		connected_component[root] = 1;
		for (int i = 0; i < n_v; i++) if (remain_unvisited[i] && connected_component[i] == 0 
				&& ((E[root * n_v + i]) || (E[i * n_v + root]))) {
			dfs_find_connected_component(i, remain_unvisited, connected_component);
		}
	}

	void iter_find_longest_path(int *remain_unvisited, int *best_n_hamiltonian_path, int *best_hamiltonian_path)
	{
		int *save_remain_unvisited = calloc(n_v, sizeof(int));
		COPY_ARR(remain_unvisited, save_remain_unvisited, n_v);
		int *connected_component = calloc(n_v, sizeof(int));
		for (int i = 0; i < n_v; i++) if (remain_unvisited[i]) {
			dfs_find_connected_component(i, remain_unvisited, connected_component);
			break;
		}
		int count_connected_compoent = 0;
		for (int i = 0; i < n_v; i++) if (remain_unvisited[i] && connected_component[i]) {
			count_connected_compoent++;
			assert(remain_unvisited[i] > 0);
			int *hamiltonian_path = calloc(n_v, sizeof(int)), n_hamiltonian_path = 0 ;
			
			find_longest_path_from_node(i, &n_hamiltonian_path, hamiltonian_path, remain_unvisited);
			if (n_hamiltonian_path > *best_n_hamiltonian_path) {
				*best_n_hamiltonian_path = n_hamiltonian_path;
				VERBOSE_FLAG(3, "new best hamilton iter %d %d\n", *best_n_hamiltonian_path, n_hamiltonian_path);
				for(int i_path = 0; i_path < n_hamiltonian_path; i_path++) {
					best_hamiltonian_path[i_path] = hamiltonian_path[i_path];
					VERBOSE_FLAG(3, "%d ", best_hamiltonian_path[i_path]);
				}
			}
			free(hamiltonian_path);
		}
		VERBOSE_FLAG(3, "count connected component %d\n", count_connected_compoent);
		VERBOSE_FLAG(1, "\n");
		for (int i = 0 ; i < n_v; i++)
			assert(save_remain_unvisited[i] ==  remain_unvisited[i]);
		free(save_remain_unvisited);
	}
	
//#############################  BEGIN OF FUNCTION  ###########################
	int thres_len = global_thres_length, thres_len_min = global_thres_length_min;
	VERBOSE_FLAG(1, "thres len MIN %d ", thres_len_min);

	int *arr_i_short = NULL, n_arr_short = 0;
	int *mark_short = NULL;
	VERBOSE_FLAG(1, "g->n_e %ld", g->n_e);
	for (int e = 0; e < g->n_e; ++e) {
		int len = get_edge_len(&g->edges[e]);
		VERBOSE_FLAG(3, "len %d\n", len);
		if (len < thres_len && len > thres_len_min) {
			++n_arr_short;
			arr_i_short = realloc(arr_i_short, n_arr_short * sizeof(int));
			mark_short = realloc(mark_short, n_arr_short * sizeof(int));
			arr_i_short[n_arr_short-1] = e;
			mark_short[n_arr_short-1] = 0;
		}
	}

	int *remain_unvisited = calloc(n_v, sizeof(int));
	for (int i = 0; i < n_v; i++) {
		float cvr = global_genome_coverage;
		VERBOSE_FLAG(3, "list v %d", listV[i]);
		float cov_times = (__get_edge_cov(&g->edges[listV[i]], g->ksize)/cvr) ;
		remain_unvisited[i] = lround(cov_times);
		VERBOSE_FLAG(3, "%d " , remain_unvisited[i]);
	}
	int count = 0;
	while (1){
		int *best_hamiltonian_path = calloc(n_v, sizeof(int)), best_n_hamiltonian_path = 0 ;
		iter_find_longest_path(remain_unvisited, &best_n_hamiltonian_path, best_hamiltonian_path);
		VERBOSE_FLAG(3, "best_n_hamiltonian_path %d\n", best_n_hamiltonian_path);
		for(int i = 0; i < best_n_hamiltonian_path; i++) {
			VERBOSE_FLAG(3, "remain %d %d\n", best_hamiltonian_path[i], remain_unvisited[best_hamiltonian_path[i]]);
			remain_unvisited[best_hamiltonian_path[i]]--;
			remain_unvisited[best_hamiltonian_path[i]^1]--;
			assert(remain_unvisited[best_hamiltonian_path[i]]>=0);
			assert(remain_unvisited[best_hamiltonian_path[i]^1]>=0);
		}
		if (best_n_hamiltonian_path == 0){
			free(best_hamiltonian_path);
			break;
		}

		int *arr_insert = NULL, n_insert = best_n_hamiltonian_path + 1;
		arr_insert = realloc(arr_insert, n_insert * sizeof(int));
		for (int i = 0; i < n_insert; i++) arr_insert[i] = -1;
		for (int i = 0; i < best_n_hamiltonian_path; i++) 
			best_hamiltonian_path[i] = listV[best_hamiltonian_path[i]];
		VERBOSE_FLAG(1, "best n hamiltonian %d \n", best_n_hamiltonian_path);
		for(int i = 0 ; i < best_n_hamiltonian_path; i++) 
			VERBOSE_FLAG(1, "best ham path %d ", best_hamiltonian_path[i]);
//		insert_short_contigs(best_n_hamiltonian_path, best_hamiltonian_path, n_insert, arr_insert, n_arr_short, arr_i_short, mark_short);
		merge_big_and_small(&best_n_hamiltonian_path, &best_hamiltonian_path, n_insert, arr_insert);

		print_contig(g, out_file, count, best_n_hamiltonian_path, best_hamiltonian_path);
		VERBOSE_FLAG(3, "contig path ");
		for (int i = 0; i < best_n_hamiltonian_path; i++) 
			VERBOSE_FLAG(3, "%d " , best_hamiltonian_path[i]);
		VERBOSE_FLAG(3, "best n %d\n", best_n_hamiltonian_path);
		count++;
		free(best_hamiltonian_path);
		free(arr_insert);
	}

	VERBOSE_FLAG(1, "n arr short %d\n", n_arr_short);
	for (int i = 0; i < n_arr_short; i++) if (remain_unvisited[i] == 0) {
		VERBOSE_FLAG(3, "printf short edge \n");
		int e = arr_i_short[i]; 
		uint32_t seq_len = 0;
		char *seq = NULL;
		int len = dump_edge_seq_reduce_N(&seq, &seq_len, &g->edges[e]);
		print_seq(out_file, count, seq, seq_len, 1); 
		count++;
		free(seq);
	}

	for (int e = 0; e < g->n_e; e++){
		int len  = get_edge_len(&g->edges[e]);
		VERBOSE_FLAG(3, "len very short %d %d\n", len, thres_len_min); 
		if (len < thres_len_min && len > 1000) {
			VERBOSE_FLAG(3, "printf very short edge \n");
			uint32_t seq_len = 0;
			char *seq = NULL;
			int len = dump_edge_seq_reduce_N(&seq, &seq_len, &g->edges[e]);
			print_seq(out_file, count, seq, seq_len, 1); 
			count++;
			free(seq);
		}
	}

	free(remain_unvisited);
	free(arr_i_short);
	free(mark_short);
}

