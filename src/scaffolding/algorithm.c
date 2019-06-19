#include <stdint.h>
#include <assert.h>
#include "verbose.h"
#include <stdlib.h>
#include <math.h>
#include "assembly_graph.h"
#include "scaffolding/output.h"
#include "scaffolding/global_params.h"
#include "scaffolding/edge.h"
#include "scaffolding/contig_graph.h"

int get_relative_cov(struct asm_graph_t *g, int i_e)
{
	float cvr = global_genome_coverage;
	return lround(__get_edge_cov(&g->edges[i_e], g->ksize)/cvr);
}


int hround(float x)
{
	return lround(x);
}

int get_rc_id(struct asm_graph_t *g, int i_e)
{
	return g->edges[i_e].rc_id;
}

int get_rc_id_V(struct asm_graph_t *g, int n_v, int *listV, int i_e)
{
	int rc = g->edges[listV[i_e]].rc_id;
	if ((i_e^1) < n_v && listV[i_e^1] == rc) {
		return i_e^1;
	} else{
		for (int i = 0; i < n_v; i++) {
			if (listV[i] == rc) {
				VERBOSE_FLAG(0, "rcidV %d %d\n", i_e, i);
				return i;
			}
		}
	}
}

int get_new_size(int size)
{
	if (size == 0)
		return 1;
	return size << 1;
}

float get_component_cov(struct asm_graph_t *g, int n_component_node, int *connected_component, int *listV)
{
	uint64_t sum_count = 0, sum_len = 0, sum_n_holes = 0;
	for (int i = 0 ; i < n_component_node; i++) {
		struct asm_edge_t *e = &g->edges[listV[connected_component[i]]];
		sum_count += e->count;
		sum_len += e->seq_len;
		sum_n_holes += e->n_holes;
	}
	return 1.0*sum_count/(sum_len - g->ksize * sum_n_holes - g->ksize * n_component_node);
}

void dfs_hamiltonian(int x, int depth, int n_adj, int *list_adj, int n_v, float *E, int *remain_unvisited, 
		int *cur_path, int *best_n_res, int *best_res, struct asm_graph_t *g, 
		int *listV)
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
	remain_unvisited[get_rc_id_V(g, n_v, listV, x)]--;
	assert(remain_unvisited[x]>=0 && remain_unvisited[get_rc_id_V(g, n_v, listV, x)]>=0);
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
				cur_path, best_n_res, best_res, g, listV);
		}
	}
	remain_unvisited[x]++;
	remain_unvisited[get_rc_id_V(g, n_v, listV, x)]++;
	free(list_sort);
}

void algo_find_hamiltonian(FILE *out_file, struct asm_graph_t *g, float *E, int n_v, 
		int *res, int *n_res, int *listV, float avg_bin_hash, struct opt_proc_t *opt)
{
	void insert_short_contigs(struct asm_graph_t *g,int *n_big_contigs, int **big_contigs, 
			struct contig_graph *contig_g , int *mark_short)
	{
		int *new_array = calloc(1, sizeof(int)), n_new_array = 1;
		new_array[0] = (*big_contigs)[0];
		for (int i = 1 ; i < *n_big_contigs; i++) {
			int *short_path  = NULL, n_short_path = 0;
			VERBOSE_FLAG(0, "find path short edge %d %d\n", (*big_contigs)[i-1], (*big_contigs)[i]);
			int src_edge = (*big_contigs)[i-1], des_edge = (*big_contigs)[i];
			struct barcode_hash_t *hl = &(g->edges[g->edges[src_edge].rc_id].bucks[0]);
			struct barcode_hash_t *hr = &(g->edges[des_edge].bucks[0]);
			find_path_short(g->edges[src_edge].target, g->edges[des_edge].source, g, contig_g, 
					&n_short_path, &short_path, mark_short, hl, hr);
			new_array = realloc(new_array, (n_new_array + n_short_path+1)* sizeof(int));
			for (int j = 0; j < n_short_path; j++){
				new_array[n_new_array + j] = short_path[j];
			}
			n_new_array += n_short_path;
			new_array[n_new_array++] = (*big_contigs)[i];
		}
		VERBOSE_FLAG(1, "big contig length %d:", *n_big_contigs);
		for (int i = 0; i < *n_big_contigs; i++){
			VERBOSE_FLAG(3, "%d ", (*big_contigs)[i]);
		}
		*n_big_contigs = n_new_array;
		*big_contigs = new_array;
//				int score0 = get_score_big_small(big_contigs[i-1], arr_i_short[j], g, avg_bin_hash);
//				int score1 = get_score_big_small(g->edges[big_contigs[i]].rc_id, g->edges[arr_i_short[j]].rc_id, g, avg_bin_hash);
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
		remain_unvisited[get_rc_id_V(g, n_v, listV, x)]--;
		assert(*best_n_hamiltonian_path == 0);
		*best_n_hamiltonian_path = 1;
		best_hamiltonian_path[0] = x;
		while (1) {
			int *list_adj = NULL, count_adj = 0;
			int last_pos = best_hamiltonian_path[*best_n_hamiltonian_path-1];
			for (int i = 0; i < n_v; i++) if (E[last_pos * n_v + i] != 0 && remain_unvisited[i]){
				assert(remain_unvisited[i] > 0);
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
					cur_path, &best_add_len, best_hamiltonian_path + *best_n_hamiltonian_path, g,
					listV
				);
				free(cur_path);
			}
			VERBOSE_FLAG(3, "best n from node n_v best add len %d %d %d\n", *best_n_hamiltonian_path, n_v, best_add_len);
			for (int i = 0; i < *best_n_hamiltonian_path + best_add_len ; i++) {
				VERBOSE_FLAG(3, "%d ", best_hamiltonian_path[i]);
			}
			for (int i_path = *best_n_hamiltonian_path; i_path < *best_n_hamiltonian_path + best_add_len; i_path++){
				remain_unvisited[best_hamiltonian_path[i_path]]--;
				remain_unvisited[get_rc_id_V(g, n_v, listV, best_hamiltonian_path[i_path])]--;
			}
			*best_n_hamiltonian_path += best_add_len;
			free(best_local_path);
			free(list_adj);
		}
		for (int i = 0; i < *best_n_hamiltonian_path; i++){
			remain_unvisited[best_hamiltonian_path[i]]++;
			remain_unvisited[get_rc_id_V(g, n_v, listV, best_hamiltonian_path[i])]++;
		}
	}

	void dfs_find_connected_component(int root, int *mark, int
			*n_component_node, int *connected_component)
	{
		connected_component[(*n_component_node)++] = root;
		mark[root] = 0;
		mark[get_rc_id_V(g, n_v, listV, root)] = 0;
		for (int i = 0; i < n_v; i++) if (mark[i] == 1 
				&& ((E[root * n_v + i]) || (E[i * n_v + root]))) {
			dfs_find_connected_component(i, mark, n_component_node, connected_component);
		}
	}

	void iter_find_longest_path(int *remain_unvisited, int *best_n_hamiltonian_path, int *best_hamiltonian_path)
	{
		for (int i = 0; i < n_v; i++) if (remain_unvisited[i]) {
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
		VERBOSE_FLAG(1, "\n");
	}
	
	void find_path_in_component(int n_component_node, int *connected_component, 
			int *count, struct contig_graph *contig_g, int *mark_short){
		int *remain_unvisited = calloc(n_v, sizeof(int));
		for (int i = 0 ; i < n_component_node; i++) {
			float cvr = 0;
			if (opt->metagenomics) 
				cvr = get_component_cov(g, n_component_node, connected_component, listV);
			else 
				cvr = global_genome_coverage;
			int i_edge = connected_component[i];
			int i_contig = listV[i_edge];

			float cov_times = (__get_edge_cov(&g->edges[i_contig], g->ksize)/cvr);
			remain_unvisited[i_edge] = MIN(hround(cov_times), 1);
			remain_unvisited[get_rc_id_V(g, n_v, listV, i_edge)] = MIN(hround(cov_times), 1);
		}
		// todo @huu calculate remain_unvisited
		while (1) {
			int *best_hamiltonian_path = calloc(n_v, sizeof(int)), best_n_hamiltonian_path = 0;
			iter_find_longest_path(remain_unvisited, &best_n_hamiltonian_path, best_hamiltonian_path);
			if (best_n_hamiltonian_path == 0){
				free(best_hamiltonian_path);
				break;
			}
			VERBOSE_FLAG(1, "calculate seq %d", *count);
			for(int i = 0; i < best_n_hamiltonian_path; i++) {
				VERBOSE_FLAG(3, "remain %d %d\n", best_hamiltonian_path[i], remain_unvisited[best_hamiltonian_path[i]]);
				remain_unvisited[best_hamiltonian_path[i]]--;
				remain_unvisited[get_rc_id_V(g, n_v, listV, best_hamiltonian_path[i])]--;
			}

			for (int i = 0; i < best_n_hamiltonian_path; i++) 
				best_hamiltonian_path[i] = listV[best_hamiltonian_path[i]];
			VERBOSE_FLAG(1, "best n hamiltonian %d \n", best_n_hamiltonian_path);
			for(int i = 0 ; i < best_n_hamiltonian_path; i++) 
				VERBOSE_FLAG(1, "best ham path %d ", best_hamiltonian_path[i]);
			VERBOSE_FLAG(1, "\n");
			insert_short_contigs(g, &best_n_hamiltonian_path, &best_hamiltonian_path, contig_g, mark_short);
			VERBOSE_FLAG(1, "best n hamiltonian %d \n", best_n_hamiltonian_path);
			for(int i = 0 ; i < best_n_hamiltonian_path; i++) 
				VERBOSE_FLAG(1, "best ham path %d ", best_hamiltonian_path[i]);
			VERBOSE_FLAG(1, "\n");

			print_contig(g, out_file, *count, best_n_hamiltonian_path, best_hamiltonian_path);
			VERBOSE_FLAG(3, "contig path ");
			for (int i = 0; i < best_n_hamiltonian_path; i++) 
				VERBOSE_FLAG(3, "%d " , best_hamiltonian_path[i]);
			(*count)++;
			free(best_hamiltonian_path);
		}
	}
//#############################  BEGIN OF FUNCTION  ###########################
	int thres_len_min = global_thres_length;
	VERBOSE_FLAG(1, "thres len MIN %d ", thres_len_min);

	int *mark = calloc(n_v, sizeof(int));
	for (int i = 0; i < n_v; i++) {
		mark[i] = 1;
	}
	int count = 0;
	int *mark_short =  calloc(g->n_e, sizeof(int));
	for (int e = 0; e < g->n_e; e++){
		int len  = get_edge_len(&g->edges[e]);
		if (len < global_thres_length ) {
			float cvr = global_genome_coverage;
			mark_short[e] = get_relative_cov(g, e);
		}
	}
	struct contig_graph *contig_g = build_graph_contig(g, mark_short);
	while (1){
		int *connected_component = calloc(n_v, sizeof(int)), n_component_node = 0;
		for (int i = 0; i < n_v; i++) if (mark[i]) {
			dfs_find_connected_component(i, mark, &n_component_node, connected_component);
			find_path_in_component(n_component_node, connected_component, &count
					, contig_g, mark_short);
			break;
		}
		if (n_component_node == 0)
			break;
	}

	for (int e = 0; e < g->n_e; e++){
		if (mark_short[e] == 1) {
			mark_short[e] = 0;
			mark_short[g->edges[e].rc_id] = 0;
			VERBOSE_FLAG(3, "printf very short edge \n");
			uint32_t seq_len = 0;
			char *seq = NULL;
			int len = dump_edge_seq_reduce_N(&seq, &seq_len, &g->edges[e]);
			print_seq(out_file, count, seq, seq_len, 1); 
			count++;
			free(seq);
		}
	}
}

