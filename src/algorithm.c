void dfs_hamiltonian(int x, uint32_t depth, uint32_t *E, int *head, uint32_t *next, int *mark, uint32_t n_v, 
		uint32_t *res, uint32_t *best_res, uint32_t *best_n_res)
{
	mark[x]--;
	mark[x^1]--;
	res[depth-1] = x;
	if (depth > *best_n_res) {
		*best_n_res = depth;
		for (uint32_t i = 0; i < depth; i++){
			best_res[i] = res[i];
		}
	} 
	for (int i = head[x]; i != -1; i = next[i]) {
		assert(mark[E[i]] >= 0);
		if (mark[E[i]] > 0) {
			dfs_hamiltonian(E[i], depth + 1, E, head, next, mark, n_v, res, best_res, best_n_res);
		}
	}
	mark[x]++;
	mark[x^1]++;
}

void algo_find_hamiltonian(FILE *out_file, struct asm_graph_t *g, uint32_t *E, uint32_t n_e, int *head, uint32_t *next, uint32_t n_v, uint32_t *res, uint32_t *n_res, uint32_t *listV)
{
	void print_contig(int index, uint32_t n_contig, uint32_t *list_contig)
	{
		char *seq = NULL, *total_seq = NULL, *NNN = NULL;
		NNN = calloc(1000, sizeof(char));
		for (uint32_t i = 0; i < 1000; i++) 
			NNN[i] = 'N';
		uint32_t seq_len = 0, total_len = 0;
		for(uint32_t i = 0; i < n_contig; i++) {
			uint32_t e = listV[list_contig[i]];
			uint32_t len = dump_edge_seq(&seq, &seq_len, &g->edges[e]);
			total_seq = realloc(total_seq, (total_len + len) * sizeof(char));
			memcpy(total_seq+total_len, seq, len);
			total_len += len;
		}
		print_seq(out_file, index, total_seq, total_len, 1);
		for(uint32_t i = 0; i < n_contig; i++) {
			uint32_t e = listV[list_contig[i]];
			__VERBOSE("%d ", e);
		}
	}

	int *mark = calloc(n_v, sizeof(uint32_t));
	for (uint32_t i = 0; i < n_v; i++) {
		float cvr = get_genome_coverage(g);
		float cov_times = (__get_edge_cov(&g->edges[listV[i]], g->ksize)/cvr) ;
		mark[i] = roundint(cov_times);
		__VERBOSE("%d " , mark[i]);
	}
	uint32_t count = 0;
	for(uint32_t ii = 0; ii < n_v; ii++){
		uint32_t *best_res = calloc(n_v, sizeof(uint32_t)), best_n_res = 0 ;
		for (uint32_t i = 0; i < n_v; i++) if (mark[i] > 0) {
			__VERBOSE("dfs from %d %d\n", i, mark[i]);
			dfs_hamiltonian(i,  1, E, head, next, mark, n_v, res, best_res, &best_n_res);
		}
		if (best_n_res == 0) 
			break;
		print_contig(ii, best_n_res, best_res);
		__VERBOSE("best n %d\n", best_n_res);
//		__VERBOSE("\nmark\n");
		for (uint32_t i = 0; i < best_n_res; i++) {
			mark[best_res[i]]--;
			mark[best_res[i]^1]--;
			//__VERBOSE("%d %d \n", listV[best_res[i]], mark[listV[best_res[i]]] );
		}
		count++;
	}
	int thres_len = global_thres_length;
	for (uint32_t e = 0; e < g->n_e; ++e) {
		int len = get_edge_len(&g->edges[e]);
		if (len < thres_len) {
			uint32_t seq_len = 0;
			char *seq = NULL;
			uint32_t len = dump_edge_seq(&seq, &seq_len, &g->edges[e]);
			print_seq(out_file, count, seq, seq_len, 1); 
			count++;
		}
	}
	fclose(out_file);
}

