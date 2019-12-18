#include "read_pairs_resolve.h"

void get_read_pairs_count(struct asm_graph_t *g, char *path,
		struct read_pair_cand_t *rp_cand)
{
	FILE *f = fopen(path, "r");
	int v, u_rc, count;
	khash_t(long_int) *rp_count = kh_init(long_int);
	while (fscanf(f, "%d %d %d\n", &v, &u_rc, &count) == 3){
		int u = g->edges[u_rc].rc_id;
		int v_rc = g->edges[v].rc_id;
		uint64_t code = GET_CODE(v, u);
		uint64_t code_rc = GET_CODE(u_rc, v_rc);
		kh_long_int_set(rp_count, code,
				kh_long_int_try_get(rp_count, code, 0) + count);
		kh_long_int_set(rp_count, code_rc,
				kh_long_int_try_get(rp_count, code_rc, 0) + count);
	}

	for (khiter_t it = kh_begin(rp_count); it != kh_end(rp_count); ++it){
		if (!kh_exist(rp_count, it))
			continue;
		uint64_t code = kh_key(rp_count, it);
		int v = code >> 32;
		int u = code & -1;
		int val = kh_val(rp_count, it);
		//if (val < MIN_READ_PAIR_MAPPED_SOFT){
		//	log_debug("Edges %d and %d share only %d read pairs",
		//			v, u, val);
		//	continue;
		//}

		rp_cand[v].cand = realloc(rp_cand[v].cand, (rp_cand[v].n + 1) * sizeof(int));
		rp_cand[v].cand[rp_cand[v].n] = u;
		rp_cand[v].score = realloc(rp_cand[v].score, (rp_cand[v].n + 1) * sizeof(int));
		rp_cand[v].score[rp_cand[v].n] = val;
		++rp_cand[v].n;
	}
	kh_destroy(long_int, rp_count);
	fclose(f);
}

int cmp_edge_length(const void *a, const void *b)
{
	int len_a = (*(uint64_t *) a) & -1;
	int len_b = (*(uint64_t *) b) & -1;
	if (len_a < len_b)
		return -1;
	if (len_a > len_b)
		return 1;
	return 0;
}

int get_next_cand(struct asm_graph_t *g, float unit_cov, struct read_pair_cand_t *rp_cand,
		int *path, int n_path)
{
	khash_t(set_int) *cand = kh_init(set_int);
	int last = path[n_path - 1];
	for (int i = 0; i < rp_cand[last].n; ++i){
		int v = rp_cand[last].cand[i];
		float cov = __get_edge_cov(g->edges + v, g->ksize);
		if (cov >= 0.25 * unit_cov
			&& rp_cand[last].score[i] >= MIN_READ_PAIR_MAPPED_HARD
			&& g->edges[v].seq_len >= 100)
			kh_set_int_insert(cand, rp_cand[last].cand[i]);
	}

	int p = rp_cand[last].n - 2;
	while (p >= 0 && kh_size(cand) > 1){
		int v = path[p];
		khash_t(set_int) *new_cand = kh_init(set_int);
		for (int i = 0; i < rp_cand[v].n; ++i){
			int u = rp_cand[v].cand[i];
			if (rp_cand[v].score[i] >= MIN_READ_PAIR_MAPPED_SOFT
				&& kh_set_int_exist(cand, u)
				&& g->edges[u].seq_len >= 100)
				kh_set_int_insert(new_cand, u);
		}
		kh_destroy(set_int, cand);
		cand = new_cand;
		--p;
	}
	int res = -1;
	if (kh_size(cand) == 1){
		int it = kh_begin(cand);
		while (kh_exist(cand, it) == 0)
			++it;
		res = kh_key(cand, it);
	}
	kh_destroy(set_int, cand);
	return res;
}

void extend_by_read_pairs(struct asm_graph_t *g, int s, float unit_cov,
		struct read_pair_cand_t *rp_cand, int **path, int *n_path)
{
	*path = calloc(1, sizeof(int));
	*n_path = 1;
	(*path)[0] = s;
	 while (1){
		int v = get_next_cand(g, unit_cov, rp_cand, *path, *n_path);
		if (v == -1)
			 break;
		int v_rc = g->edges[v].rc_id;
		*path = realloc(*path, ((*n_path) + 1) * sizeof(int));
		(*path)[(*n_path)++] = v;
	 }

	 for (int i = 0; i < *n_path; ++i){
		int v = (*path)[i];
		int v_rc = g->edges[v].rc_id;
		int count = min(unit_cov * (g->edges[v].seq_len - g->ksize + 1),
				g->edges[v].count);
		g->edges[v].count -= count;
		g->edges[v_rc].count -= count;
	 }
}

void get_long_contigs(struct opt_proc_t *opt)
{
	struct asm_graph_t *g = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g, opt->in_file);
	struct read_pair_cand_t *rp_cand = calloc(g->n_e, sizeof(struct read_pair_cand_t));
	get_read_pairs_count(g, opt->in_fasta, rp_cand);

	uint64_t *edge_sorted = calloc(g->n_e, sizeof(uint64_t));
	for (int i = 0; i < g->n_e; ++i)
		edge_sorted[i] = GET_CODE(i, g->edges[i].seq_len);
	qsort(edge_sorted, g->n_e, sizeof(uint64_t), &cmp_edge_length);

	float unit_cov = get_genome_coverage(g);
	int *visited = calloc(g->n_e, sizeof(int));
	struct asm_graph_t *g_new = calloc(1, sizeof(struct asm_graph_t));
	int n_e = 0;
	struct asm_edge_t *edges = NULL;
	for (int i = g->n_e - 1, j = 0; i >= 0 && j < 10; --i, ++j){
		int v = edge_sorted[i] >> 32;
		float cov = __get_edge_cov(g->edges + v, g->ksize);
		if (cov < 0.25 * unit_cov || g->edges[v].seq_len < 100)
			continue;
		int *path_fw, *path_rv;
		int n_path_fw, n_path_rv;
		extend_by_read_pairs(g, v, unit_cov, rp_cand, &path_fw, &n_path_fw);
		extend_by_read_pairs(g, g->edges[v].rc_id, unit_cov, rp_cand,
				&path_rv, &n_path_rv);
		int *path = calloc(n_path_fw + n_path_rv - 1, sizeof(int));
		int n_path = 0;
		for (int i = n_path_rv - 1; i >= 0; --i)
			path[n_path++] = g->edges[path_rv[i]].rc_id;
		for (int i = 1; i < n_path_fw; ++i)
			path[n_path++] = path_fw[i];
		free(path_fw);
		free(path_rv);
		/*__VERBOSE("Start from %d, fw: %d, rv: %d\n", v, n_path_fw,
				n_path_rv);
		for (int i = 0; i < n_path; ++i)
			__VERBOSE("%d ", path[i]);
		__VERBOSE("\n");*/
	}

	asm_graph_destroy(g);
	free(g);
}
