#include "read_pairs_resolve.h"
#define GET_CODE(v, u) ((((uint64_t) (v)) << 32) | (u))

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
		if (val < MIN_READ_PAIR_MAPPED){
			log_debug("Edges %d and %d share only %d read pairs",
					v, u, val);
			continue;
		}

		rp_cand[v].cand = realloc(rp_cand[v].cand, (rp_cand[v].n + 1) * sizeof(int));
		rp_cand[v].cand[rp_cand[v].n++] = u;
	}
	kh_destroy(long_int, rp_count);
	fclose(f);
}

void get_long_contigs(struct opt_proc_t *opt)
{
	struct asm_graph_t *g = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g, opt->in_file);
	struct read_pair_cand_t *rp_cand = calloc(g->n_e, sizeof(struct read_pair_cand_t));
	get_read_pairs_count(g, opt->in_fasta, rp_cand);

	asm_graph_destroy(g);
	free(g);
}
