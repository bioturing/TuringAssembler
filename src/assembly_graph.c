#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "io_utils.h"
#include "k31hash.h"
#include "k63hash.h"
#include "k31_count.h"
#include "k63_count.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

#define TIPS_THRESHOLD			0.1

void remove_tips(struct asm_graph_t *g0, struct asm_graph_t *g);
void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void write_gfa(struct asm_graph_t *g, const char *path);

void k63_process(struct opt_count_t *opt)
{
	char path[1024];
	init_clock();

	struct k63hash_t *kmer_hash;
	kmer_hash = calloc(1, sizeof(struct k63hash_t));
	__VERBOSE("Estimating kmer\n");
	build_k63_table_lazy(opt, kmer_hash, opt->kmer_slave);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Estimating kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nBuilding assembly graph\n");
	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	build_asm_graph_from_k63(opt, opt->kmer_slave, kmer_hash, g0);
	__VERBOSE_LOG("Graph #1", "Number of nodes: %lld\n", (long long)g0->n_v);
	__VERBOSE_LOG("Graph #1", "Number of edges: %lld\n", (long long)g0->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_0.gfa");
	write_gfa(g0, path);

	__VERBOSE("\nRemoving tips\n");
	struct asm_graph_t *g1;
	g1 = calloc(1, sizeof(struct asm_graph_t));
	remove_tips(g0, g1);
	__VERBOSE_LOG("Graph #2", "Number of nodes: %lld\n", (long long)g1->n_v);
	__VERBOSE_LOG("Graph #2", "Number of edges: %lld\n", (long long)g1->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k63_1.gfa");
	write_gfa(g1, path);

}

void k31_process(struct opt_count_t *opt)
{
	char path[1024];
	init_clock();

	struct k31hash_t *kmer_hash;
	kmer_hash = calloc(1, sizeof(struct k31hash_t));
	__VERBOSE("Estimating kmer\n");
	build_k31_table_lazy(opt, kmer_hash, opt->kmer_slave);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Estimating kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nBuilding assembly graph\n");
	struct asm_graph_t *g0;
	g0 = calloc(1, sizeof(struct asm_graph_t));
	build_asm_graph_from_k31(opt, opt->kmer_slave, kmer_hash, g0);
	__VERBOSE_LOG("Graph #1", "Number of nodes: %lld\n", (long long)g0->n_v);
	__VERBOSE_LOG("Graph #1", "Number of edges: %lld\n", (long long)g0->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_0.gfa");
	write_gfa(g0, path);

	__VERBOSE("\nRemoving tips\n");
	struct asm_graph_t *g1;
	g1 = calloc(1, sizeof(struct asm_graph_t));
	remove_tips(g0, g1);
	__VERBOSE_LOG("Graph #2", "Number of nodes: %lld\n", (long long)g1->n_v);
	__VERBOSE_LOG("Graph #2", "Number of edges: %lld\n", (long long)g1->n_e);
	strcpy(path, opt->out_dir); strcat(path, "/graph_k31_1.gfa");
	write_gfa(g1, path);
}

static uint32_t *append_bin_seq(uint32_t *dst, gint_t dlen, uint32_t *src,
					gint_t slen, int skip)
{
	uint32_t *new_ptr;
	gint_t cur_m, m, i, k;
	cur_m = (dlen + 15) >> 4;
	m = (dlen + slen + 15) >> 4;
	if (m > cur_m) {
		new_ptr = realloc(dst, m * sizeof(uint32_t));
		memset(new_ptr + cur_m, 0, (m - cur_m) * sizeof(uint32_t));
	} else {
		new_ptr = dst;
	}
	for (i = skip; i < skip + slen; ++i) {
		k = dlen + (i - skip);
		new_ptr[k >> 4] |= ((src[i >> 4] >> ((i & 15) << 1)) & (uint32_t)3) << ((k & 15) << 1);
	}
	return new_ptr;
}

static uint32_t *asm_clone_binseq(uint32_t *src, gint_t len)
{
	uint32_t *ret;
	ret = calloc((len + 15) >> 4, sizeof(uint32_t));
	memcpy(ret, src, ((len + 15) >> 4) * sizeof(uint32_t));
	return ret;
}

static int asm_find_edge_index(gint_t *adj, int deg, gint_t id)
{
	int i, ret;
	ret = -1;
	for (i = 0; i < deg; ++i) {
		if (adj[i] == id)
			ret = i;
	}
	return ret;
}

static gint_t asm_only_positive_edge(gint_t *adj, int deg)
{
	int i, ret;
	ret = -1;
	for (i = 0; i < deg; ++i) {
		if (adj[i] != -1) {
			if (ret == -1)
				ret = adj[i];
			else
				ret = -2;
		}
	}
	return ret;
}

static void dump_bin_seq(char *seq, uint32_t *bin, gint_t len)
{
	gint_t i;
	for (i = 0; i < len; ++i)
		seq[i] = nt4_char[(bin[i >> 4] >> ((i & 15) << 1)) & 3];
	seq[len] = '\0';
}

#define __bin_seq_get_char(seq, l) (((seq)[(l) >> 4] >> (((l) & 15) << 1)) & (uint32_t)0x3)

int asm_is_edge_rc(uint32_t *seq1, gint_t l1, uint32_t *seq2, gint_t l2)
{
	if (l1 != l2)
		return 0;
	gint_t i, k;
	uint32_t c1, c2;
	for (i = 0; i < l1; ++i) {
		k = l1 - i - 1;
		c1 = __bin_seq_get_char(seq1, i);
		c2 = __bin_seq_get_char(seq2, k);
		if (c1 != (c2 ^ 3)) {
			return 0;
		}
	}
	return 1;
}

void write_gfa(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "w");
	char *seq = NULL;
	int seq_len = 0;
	gint_t e, e_rc;
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		int deg = g->nodes[g->edges[e].target].deg +
			g->nodes[g->edges[e_rc].target].deg;
		if (deg == 0)
			continue;
		if (seq_len < g->edges[e].seq_len + 1) {
			seq_len = g->edges[e].seq_len + 1;
			seq = realloc(seq, seq_len);
		}
		dump_bin_seq(seq, g->edges[e].seq, g->edges[e].seq_len);
		fprintf(fp, "S\t%lld\t%s\tKC:i:%llu\n", (long long)e, seq,
					(long long unsigned)g->edges[e].count);
	}
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		gint_t pe, next_pe;
		char ce, next_ce;
		if (e > e_rc) {
			pe = e_rc;
			ce = '-';
		} else {
			pe = e;
			ce = '+';
		}
		gint_t n = g->edges[e].target;
		int k;
		for (k = 0; k < g->nodes[n].deg; ++k) {
			gint_t next_e, next_e_rc;
			next_e = g->nodes[n].adj[k];
			next_e_rc = g->edges[next_e].rc_id;
			if (next_e > next_e_rc) {
				next_pe = next_e_rc;
				next_ce = '-';
			} else {
				next_pe = next_e;
				next_ce = '+';
			}
			fprintf(fp, "L\t%lld\t%c\t%lld\t%c\t%dM\n",
				(long long)pe, ce, (long long)next_pe, next_ce,
				g->ksize);
		}
	}
	fclose(fp);
}

void remove_tips(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	gint_t u;
	for (u = 0; u < g0->n_v; ++u) {
		double max_cov = 0.0, cov;
		int deg = g0->nodes[u].deg, c;
		gint_t e_id, e_rc;
		for (c = 0; c < deg; ++c) {
			e_id = g0->nodes[u].adj[c];
			if (e_id == -1)
				continue;
			cov = g0->edges[e_id].count * 1.0 /
					(g0->edges[e_id].seq_len - g0->ksize);
			max_cov = __max(max_cov, cov);
		}
		for (c = 0; c < deg; ++c) {
			e_id = g0->nodes[u].adj[c];
			if (e_id == -1)
				continue;
			cov = g0->edges[e_id].count * 1.0 /
					(g0->edges[e_id].seq_len - g0->ksize);
			if (cov / max_cov < TIPS_THRESHOLD) {
				gint_t v, v_rc;
				v = g0->edges[e_id].target;
				v_rc = g0->nodes[v].rc_id;
				e_rc = g0->edges[e_id].rc_id;
				int j = asm_find_edge_index(g0->nodes[v_rc].adj,
						g0->nodes[v_rc].deg, e_rc);
				assert(j != -1);
				/* disconnect edge */
				g0->nodes[u].adj[c] = -1;
				g0->nodes[v_rc].adj[j] = -1;
			}
		}
		// g0->nodes[u].deg = 0;
		// for (c = 0; c < deg; ++c) {
		// 	gint_t e_id = g0->nodes[u].adj[c];
		// 	if (e_id != -1)
		// 		g0->nodes[u].adj[g0->nodes[u].deg++] = e_id;
		// }
		// g0->nodes[u].adj = realloc(g0->nodes[u].adj,
		// 		g0->nodes[u].deg * sizeof(gint_t));
	}
	for (u = 0; u < g0->n_v; ++u) {
		int c, deg;
		deg = g0->nodes[u].deg;
		g0->nodes[u].deg = 0;
		for (c = 0; c < deg; ++c) {
			gint_t e_id = g0->nodes[u].adj[c];
			if (e_id != -1)
				g0->nodes[u].adj[g0->nodes[u].deg++] = e_id;
		}
		g0->nodes[u].adj = realloc(g0->nodes[u].adj,
					g0->nodes[u].deg * sizeof(gint_t));
	}

	asm_condense(g0, g);
}

void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	gint_t *node_id, *edge_id;
	gint_t n_v, n_e;
	gint_t u;
	node_id = malloc(g0->n_v * sizeof(gint_t));
	edge_id = malloc(g0->n_e * sizeof(gint_t));
	memset(node_id, 255, g0->n_v * sizeof(gint_t));
	memset(edge_id, 255, g0->n_e * sizeof(gint_t));

	n_v = n_e = 0;
	for (u = 0; u < g0->n_v; ++u) {
		int deg_fw = g0->nodes[u].deg;
		int deg_rv = g0->nodes[g0->nodes[u].rc_id].deg;
		if ((deg_fw == 1 && deg_rv == 1) || deg_fw + deg_rv == 0)
			continue;
		node_id[u] = n_v++;
		n_e += deg_fw;
	}

	struct asm_node_t *nodes = calloc(n_v, sizeof(struct asm_node_t));
	struct asm_edge_t *edges = calloc(n_e, sizeof(struct asm_edge_t));
	n_e = 0;
	for (u = 0; u < g0->n_v; ++u) {
		gint_t new_id = node_id[u];
		if (new_id == -1)
			continue;
		nodes[new_id].adj = malloc(g0->nodes[u].deg * sizeof(gint_t));
		nodes[new_id].deg = 0;
		int c;
		for (c = 0; c < g0->nodes[u].deg; ++c) {
			gint_t e_id = g0->nodes[u].adj[c], e_rc, n_id;
			assert(e_id != -1);
			edges[n_e].seq = asm_clone_binseq(g0->edges[e_id].seq,
						g0->edges[e_id].seq_len);
			// edges[n_e].seq = g0->edges[e_id].seq;
			// g0->edges[e_id].seq = NULL;
			edges[n_e].seq_len = g0->edges[e_id].seq_len;
			edges[n_e].count = g0->edges[e_id].count;
			edges[n_e].source = new_id;

			do {
				edge_id[e_id] = n_e;
				n_id = g0->edges[e_id].target;
				// fprintf(stderr, "n_id = %lld; e_id = %lld\n", n_id, e_id);
				if (node_id[n_id] == -1) { /* middle node */
					// e_id = asm_positive_edge(g0->nodes[n_id].adj,
					// 		g0->nodes[n_id].deg);
					assert(g0->nodes[n_id].deg == 1);
					e_id = g0->nodes[n_id].adj[0];
					assert(e_id != -1);
					edges[n_e].seq =
						append_bin_seq(edges[n_e].seq,
							edges[n_e].seq_len,
							g0->edges[e_id].seq,
							g0->edges[e_id].seq_len
								- g0->ksize,
							g->ksize);
					edges[n_e].seq_len +=
						g0->edges[e_id].seq_len - g0->ksize;
					edges[n_e].count += g0->edges[e_id].count;
				} else {
					e_id = -1;
				}
			} while (e_id != -1);
			assert(node_id[n_id] != -1);
			edges[n_e].target = node_id[n_id];
			nodes[new_id].adj[nodes[new_id].deg++] = n_e;
			e_id = g0->nodes[u].adj[c];
			e_rc = g0->edges[e_id].rc_id;
			if (edge_id[e_rc] != -1) {
				edges[n_e].rc_id = edge_id[e_rc];
				edges[edge_id[e_rc]].rc_id = n_e;
			}
			++n_e;
		}
		nodes[new_id].adj = realloc(nodes[new_id].adj, nodes[new_id].deg * sizeof(gint_t));
	}

	fprintf(stderr, "Done condense edges\n");

	free(node_id);
	free(edge_id);
	g->ksize = g0->ksize;
	g->n_v = n_v;
	g->n_e = n_e;
	g->nodes = nodes;
	g->edges = edges;
}

void test_asm_graph(struct asm_graph_t *g)
{
}
