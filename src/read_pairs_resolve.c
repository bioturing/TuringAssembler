#include <stdlib.h>
#include "minimizers/count_barcodes.h"
#include "read_pairs_resolve.h"
#include "barcode_graph.h"

void read_pair_cand_destroy(struct read_pair_cand_t *rp_cand)
{
	free(rp_cand->cand);
	free(rp_cand->score);
}

void get_read_pairs_count(struct asm_graph_t *g, char *path,
		struct read_pair_cand_t *rp_cand)
{
	FILE *f = fopen(path, "r");
	int v, u_rc, count;
	khash_t(long_int) *rp_count = kh_init(long_int);
	while (fscanf(f, "%d %d %d\n", &v, &u_rc, &count) == 3){
		int u = g->edges[u_rc].rc_id;
		int v_rc = g->edges[v].rc_id;
		if (u == g->edges[v_rc].rc_id)
			continue;
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

int cmp_edge_length(const void *a, const void *b, void *arg)
{
	struct asm_graph_t *g = (struct asm_graph_t *) arg;
	int a_id = *(int *) a;
	int b_id = *(int *) b;
	if (g->edges[a_id].seq_len < g->edges[b_id].seq_len)
		return -1;
	if (g->edges[a_id].seq_len > g->edges[b_id].seq_len)
		return 1;
	return 0;
}

int get_next_cand(struct asm_graph_t *g, float unit_cov, struct read_pair_cand_t *rp_cand,
		khash_t(long_int) *share_bc, int *path, int n_path)
{
	int res = -1;
	khash_t(set_int) *cand = kh_init(set_int);
	int last = path[n_path - 1], count = 0;
	int second_score = 0, best_score = 0, best = -1;
	log_debug("Get next candidates of %d", last);
	for (int i = 0; i < rp_cand[last].n; ++i){
		int v = rp_cand[last].cand[i];
		int score = rp_cand[last].score[i];
		if (g->edges[v].rc_id == last)
			continue;
		float cov = __get_edge_cov(g->edges + v, g->ksize);
		int added = 0;
		if (cov >= 0.5 * unit_cov && g->edges[v].seq_len >= 100){
			if (rp_cand[last].score[i] > second_score) {
				second_score = rp_cand[last].score[i];
				if (second_score > best_score) {
					int tmp = best_score;
					best_score = second_score;
					second_score = tmp;
					best= rp_cand[last].cand[i];
				}
			}

			if (score >= MIN_READ_PAIR_MAPPED_HARD)
				added = 1;
		}
		if (added){
			log_debug("Add edge %d to list of candidates, score: %d",
					v, score);
			kh_set_int_insert(cand, v);
		} else {
			log_debug("%d does not sastisfy threshold, cov: %.3f, len: %.3f, score: %d, unit cov: %.3f",
					v, cov, g->edges[v].seq_len, score, unit_cov);
		}
	}
	if (best_score >  (second_score + 10)*1.3){
		res = best;
		log_debug("Edge %d has score ignificantly greater than the others, best score: second score:",
				best, best_score, second_score);
		goto end_function;
	}

	//for (khiter_t it = kh_begin(cand); it != kh_end(cand); ++it){
	//	if (!kh_exist(cand, it))
	//		continue;
	//	int v = kh_key(cand, it);
	//	int count = 0;
	//	for (int i = 0; i < rp_cand[v].n; ++i){
	//		int u = rp_cand[v].cand[i];
	//		/*if (v == 43252)
	//			__VERBOSE("%d ", u);
	//		__VERBOSE("\n");*/
	//		if (u != v && kh_set_int_exist(cand, u))
	//			++count;
	//	}
	//	if (count == kh_size(cand) - 1){
	//		res = v;
	//}
	//goto end_function;

	int p = n_path - 2;
	while (p >= 0 && kh_size(cand) > 1){
		int v = path[p];
		khash_t(set_int) *new_cand_rp = kh_init(set_int);
		for (int i = 0; i < rp_cand[v].n; ++i){
			int u = rp_cand[v].cand[i];
			if (rp_cand[v].score[i] >= MIN_READ_PAIR_MAPPED_SOFT
				&& kh_set_int_exist(cand, u)
				&& g->edges[u].seq_len >= 100)
				kh_set_int_insert(new_cand_rp, u);
		}

		khash_t(set_int) *new_cand_bc = kh_init(set_int);
		for (khiter_t it = kh_begin(cand); it != kh_end(cand); ++it){
			if (!kh_exist(cand, it))
				continue;
			int u = kh_key(cand, it);
			int shared_bc = get_share_bc(g, share_bc, v, u);
			if (shared_bc >= MIN_SHARE_BC_RP_EXTEND)
				kh_set_int_insert(new_cand_bc, u);
		}
		kh_destroy(set_int, cand);
		if (kh_size(new_cand_rp) != 0){
			cand = new_cand_rp;
			kh_destroy(set_int, new_cand_bc);
		} else {
			cand = new_cand_bc;
			kh_destroy(set_int, new_cand_rp);
		}
		--p;
	}
	if (kh_size(cand) == 1){
		khiter_t it = kh_begin(cand);
		while (kh_exist(cand, it) == 0)
			++it;
		res = kh_key(cand, it);
	}
end_function:
	kh_destroy(set_int, cand);
	return res;
}

int check_good_cand(struct asm_graph_t *g, int *path, int n_path,
		khash_t(long_int) *share_bc, int cand)
{
	int last_e = path[n_path - 1];
	int sum_share = 0;
	for (int i = max(0, n_path - 4); i < n_path - 1; ++i){
		int e = path[i];
		int val = get_share_bc(g, share_bc, e, last_e);
		sum_share += val;
	}
	return sum_share < get_share_bc(g, share_bc, last_e, cand);
}

void extend_by_read_pairs(struct asm_graph_t *g, int s, float unit_cov,
		struct read_pair_cand_t *rp_cand, khash_t(long_int) *share_bc,
		int **path, int *n_path)
{
	//__VERBOSE("s = %d\n", s);
	*path = calloc(1, sizeof(int));
	*n_path = 1;
	(*path)[0] = s;
	int count = min(unit_cov * (g->edges[s].seq_len - g->ksize + 1),
			g->edges[s].count);
	g->edges[s].count -= count;
	g->edges[g->edges[s].rc_id].count -= count;

	while (1){
		int v = get_next_cand(g, unit_cov, rp_cand, share_bc, *path,
				*n_path);
		if (__get_edge_cov(&g->edges[v], g->ksize) > REPEAT_COV_RATIO * unit_cov
			|| __get_edge_cov(&g->edges[s], g->ksize) > REPEAT_COV_RATIO * unit_cov)
			return;
		if (v == -1)
			return;
		if (__get_edge_cov(&g->edges[v], g->ksize) > REPEAT_COV_RATIO * unit_cov
			|| __get_edge_cov(&g->edges[s], g->ksize) > REPEAT_COV_RATIO * unit_cov){
			log_debug("Next cand %d is repeat, v cov: %.3f, s cov: %.3f, unit cov: %.3f",
					v, __get_edge_cov(g->edges + v, g->ksize),
					__get_edge_cov(g->edges + s, g->ksize), unit_cov);
			return;
		}
		log_debug("Next cand of %d: %d", (*path)[(*n_path) - 1], v);
		int v_rc = g->edges[v].rc_id;
		int count = min(unit_cov * (g->edges[v].seq_len - g->ksize + 1),
				g->edges[v].count);
		g->edges[v].count = 0;
		g->edges[v_rc].count = 0;
		*path = realloc(*path, ((*n_path) + 1) * sizeof(int));
		(*path)[(*n_path)++] = v;
	}
}

int get_share_bc(struct asm_graph_t *g, khash_t(long_int) *share_bc, int u, int v)
{
	uint64_t code = get_min_code(g, u, v);
	gint_t k = kh_get(long_int, share_bc, code);
	int ret = kh_value(share_bc,k);
	return ret;
}

void join_read_pair_path(struct asm_graph_t *g, int *path, int n_path,
		char **seq)
{
	decode_seq(seq, g->edges[path[0]].seq, g->edges[path[0]].seq_len);
	int seq_len = strlen(*seq);
	for (int i = 1; i < n_path; ++i){
		int e = path[i];
		//if (path[0] == 50913)
		//	__VERBOSE("%d\n", e);
		char *tmp_seq;
		decode_seq(&tmp_seq, g->edges[e].seq, g->edges[e].seq_len);
		char N[51] = {};
		for (int j = 0; j < 50; ++j)
			N[j] = 'N';
		int new_len = seq_len + 50 + g->edges[e].seq_len;
		*seq = realloc(*seq, (new_len + 1) * sizeof(char));
		strcat(*seq, N);
		strcat(*seq, tmp_seq);
		free(tmp_seq);
		seq_len = new_len;
	}
}

void get_long_contigs(struct opt_proc_t *opt)
{
	struct asm_graph_t *g = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g, opt->in_file);
	struct read_pair_cand_t *rp_cand = calloc(g->n_e, sizeof(struct read_pair_cand_t));
	get_read_pairs_count(g, opt->in_fasta, rp_cand);

	int *edge_sorted = calloc(g->n_e, sizeof(int));
	for (int i = 0; i < g->n_e; ++i)
		edge_sorted[i] = i;
	qsort_r(edge_sorted, g->n_e, sizeof(int), &cmp_edge_length, g);

	khash_t(long_int) *share_bc = load_khash(opt->var[0]);

	float unit_cov = get_genome_coverage_h(g);
	log_info("global cov: %lf", unit_cov);
	int *visited = calloc(g->n_e, sizeof(int));
	int n_e = 0;
	char out_path[1024];
	sprintf(out_path, "%s/graph_k_%d_extend.fasta", opt->out_dir, g->ksize);
	FILE *f = fopen(out_path, "w");
	for (int i = g->n_e - 1; i >= 0; --i){
		int p = g->n_e - i;
		if (p * 100 / g->n_e > (p - 1) * 100 / g->n_e)
			log_info("Processed %d/%d edges (%d\%)", p, g->n_e,
					p * 100 / g->n_e);
		int e = edge_sorted[i];
		log_debug("Trying to extend from %d", e);
		float cov = __get_edge_cov(g->edges + e, g->ksize);
		//if (cov < 0.5 * unit_cov || g->edges[e].seq_len < 100 || cov > 1.3 * unit_cov)
		//	continue;
		if (cov < 0.5 * unit_cov || g->edges[e].seq_len < 100){
			log_debug("Edge violates threshold, cov: %.3f, len: %d, unit cov: %.3f",
					cov, g->edges[e].seq_len, unit_cov);
			continue;
		}
		int *path_fw, *path_rv;
		int n_path_fw, n_path_rv;
		extend_by_read_pairs(g, e, unit_cov, rp_cand, share_bc, &path_fw,
				&n_path_fw);
		extend_by_read_pairs(g, g->edges[e].rc_id, unit_cov, rp_cand,
				share_bc, &path_rv, &n_path_rv);
		int *path = calloc(n_path_fw + n_path_rv - 1, sizeof(int));
		int n_path = 0;
		for (int i = n_path_rv - 1; i >= 0; --i)
			path[n_path++] = g->edges[path_rv[i]].rc_id;
		for (int i = 1; i < n_path_fw; ++i)
			path[n_path++] = path_fw[i];
		//for (int i = 0; i < n_path; ++i)
		//	path[i] = g->edges[path[i]].rc_id;
		//for (int i = 0; i < n_path; ++i){
		//	int e = path[i];
		//	int e_rc = g->edges[e].rc_id;
		//	int min_count = min(g->edges[e].count, g->edges[e_rc].count);
		//	g->edges[e].count = min_count;
		//	g->edges[e_rc].count = min_count;
		//}
		free(path_fw);
		free(path_rv);

		char *seq;
		join_read_pair_path(g, path, n_path, &seq);

		char path_str[2000] = {};
		for (int i = 0; i < n_path; ++i){
			visited[path[i]] = visited[g->edges[path[i]].rc_id] = 1;
			char tmp[1000];
			sprintf(tmp, "%d ", path[i]);
			strcat(path_str, tmp);
		}
		log_debug("Path number %d start from %d: %s", n_e, e, path_str);
		free(path);
		fprintf(f, ">SEQ_%d\n", n_e++);
		fprintf(f, "%s\n", seq);
		free(seq);
	}
	free(edge_sorted);

	for (int i = 0; i < g->n_e; ++i)
		read_pair_cand_destroy(rp_cand + i);

	for (int e = 0; e < g->n_e; ++e){
		int e_rc = g->edges[e].rc_id;
		if (e > e_rc || g->edges[e].seq_len < MIN_NOTICE_LEN)
			continue;
		float cov = __get_edge_cov(g->edges + e, g->ksize);
		if (visited[e] == 0 || cov >= 0.2 * unit_cov){
			char *seq;
			decode_seq(&seq, g->edges[e].seq, g->edges[e].seq_len);
			fprintf(f, ">SEQ_%d\n%s\n", n_e++, seq);
			free(seq);
		}
	}
	fclose(f);
	asm_graph_destroy(g);
	free(g);
}


