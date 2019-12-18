#include "read_pairs_resolve.h"
#define GET_CODE(v, u) ((((uint64_t) (v)) << 32) | (u))

void get_read_pairs_count(struct asm_graph_t *g, char *path, khash_t(long_int) *rp_count)
{
	FILE *f = fopen(path, "r");
	int v, u_rc, count;
	while (fscanf(f, "%d %d %d\n", &v, &u_rc, &count) == 3){
		if (count < MIN_READ_PAIR_MAPPED)
			continue;
		int u = g->edges[u_rc].rc_id;
		int v_rc = g->edges[v].rc_id;
		uint64_t code = GET_CODE(v, u);
		uint64_t code_rc = GET_CODE(u_rc, v_rc);
		kh_long_int_set(rp_count, code,
				kh_long_int_try_get(rp_count, code, 0) + count);
		kh_long_int_set(rp_count, code_rc,
				kh_long_int_try_get(rp_count, code_rc, 0) + count);
	}
	fclose(f);
}

void get_long_contigs(struct opt_proc_t *opt)
{
	struct asm_graph_t *g = calloc(1, sizeof(struct asm_graph_t));
	load_asm_graph(g, opt->in_file);
	khash_t(long_int) *rp_count = kh_init(long_int);
	get_read_pairs_count(g, opt->in_fasta, rp_count);

	asm_graph_destroy(g);
	free(g);
}
