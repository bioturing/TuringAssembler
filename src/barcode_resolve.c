#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "io_utils.h"
#include "resolve.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

#define __diff(g, e1, e2) ((e1) != (e2) && (g)->edges[e1].rc_id != (e2))

static inline int test_bridge(struct asm_graph_t *g, gint_t e)
{
	if (g->edges[e].source == -1)
		return 0;
	gint_t e_rc, u, v, u_rc, v_rc;
	e_rc = g->edges[e].rc_id;
	u = g->edges[e_rc].target;
	v = g->edges[e].target;
	u_rc = g->nodes[u].rc_id;
	v_rc = g->nodes[v].rc_id;
	if (g->nodes[u].deg != 2 || g->nodes[v].deg != 2 ||
		g->nodes[u_rc].deg != 1 || g->nodes[v_rc].deg != 1)
		return 0;
	fprintf(stderr, "Bridge: %ld_%ld\n", e, e_rc);
	int *test_result = alloca(4 * sizeof(int));
	int positive = 0;
	int i, k;
	for (i = 0; i < 2; ++i) {
		for (k = 0; k < 2; ++k) {
			test_result[i * 2 + k] = test_edge_barcode(g,
				g->nodes[u].adj[i], g->nodes[v].adj[k]);
			if (test_result[i * 2 + k] == 1) {
				++positive;
				// fprintf(stderr, "Path [%ld] -> [%ld]: positive\n",
				// 	g->nodes[u].adj[i], g->nodes[v].adj[k]);
				// print_test_barcode_edge(g, g->nodes[u].adj[i],
				// 			g->nodes[v].adj[k]);
			} else if (test_result[i * 2 + k] == -1) {
				// fprintf(stderr, "Path [%ld] -> [%ld]: not-consider\n",
				// 	g->nodes[u].adj[i], g->nodes[v].adj[k]);
			} else {
				// fprintf(stderr, "Path [%ld] -> [%ld]: negative\n",
				// 	g->nodes[u].adj[i], g->nodes[v].adj[k]);
				// print_test_barcode_edge(g, g->nodes[u].adj[i],
				// 			g->nodes[v].adj[k]);
			}
		}
	}
	fprintf(stderr, "Number of positive path: %d\n", positive);
	if (positive == 0 || positive > 2)
		return 0;
	gint_t idx_i, idx_k;
	for (i = 0; i < 2; ++i) {
		for (k = 0; k < 2; ++k) {
			if (test_result[i * 2 + k] == 1) {
				idx_i = i;
				idx_k = k;
			}
		}
	}
	/* join edges */
	gint_t el1, el2, er1, er2, el_rc1, el_rc2, er_rc1, er_rc2,
	       el1_len, el2_len, er1_len, er2_len;
	float cov1, cov2;
	el_rc1 = g->nodes[u].adj[idx_i];
	el1 = g->edges[el_rc1].rc_id;
	er1 = g->nodes[v].adj[idx_k];
	er_rc1 = g->edges[er1].rc_id;
	el1_len = get_edge_len(g->edges + el1);
	er1_len = get_edge_len(g->edges + er1);
	cov1 = (el1_len * __get_edge_cov(g->edges + el1, g->ksize) +
		er1_len * __get_edge_cov(g->edges + er1, g->ksize)) /
		(el1_len + er1_len);

	el_rc2 = g->nodes[u].adj[idx_i ^ 1];
	el2 = g->edges[el_rc2].rc_id;
	er2 = g->nodes[v].adj[idx_k ^ 1];
	er_rc2 = g->edges[er2].rc_id;
	el2_len = get_edge_len(g->edges + el2);
	er2_len = get_edge_len(g->edges + er2);
	cov2 = (el2_len * __get_edge_cov(g->edges + el2, g->ksize) +
		er2_len * __get_edge_cov(g->edges + er2, g->ksize)) /
		(el2_len + er2_len);

	/* check for auto loop here */
	if (!__diff(g, el1, el2) || !__diff(g, er1, er2) || !__diff(g, el1, er1)
		|| !__diff(g, el1, er2) || !__diff(g, el2, er1) ||
		!__diff(g, el2, er2))
		return 0;

	fprintf(stderr, "joining edge %ld->%ld->%ld\n", el1, e, er1);
	asm_append_edge_seq(g->edges + el1, g->edges + e, g->ksize);
	asm_append_edge_seq(g->edges + el1, g->edges + er1, g->ksize);
	g->edges[el1].count += g->edges[er1].count;
	g->edges[el1].target = g->edges[er1].target;
	asm_append_edge_seq(g->edges + er_rc1, g->edges + e_rc, g->ksize);
	asm_append_edge_seq(g->edges + er_rc1, g->edges + el_rc1, g->ksize);
	g->edges[er_rc1].count += g->edges[el_rc1].count;
	g->edges[er_rc1].target += g->edges[el_rc1].target;
	g->edges[el1].rc_id = er_rc1;
	g->edges[er_rc1].rc_id = el1;

	fprintf(stderr, "joining edge %ld->%ld->%ld\n", el2, e, er2);
	asm_append_edge_seq(g->edges + el2, g->edges + e, g->ksize);
	asm_append_edge_seq(g->edges + el2, g->edges + er2, g->ksize);
	g->edges[el2].count += g->edges[er2].count;
	g->edges[el2].target = g->edges[er2].target;
	asm_append_edge_seq(g->edges + er_rc2, g->edges + e_rc, g->ksize);
	asm_append_edge_seq(g->edges + er_rc2, g->edges + el_rc2, g->ksize);
	g->edges[er_rc2].count += g->edges[el_rc2].count;
	g->edges[er_rc2].target = g->edges[el_rc2].target;
	g->edges[el2].rc_id = er_rc2;
	g->edges[er_rc2].rc_id = el2;

	uint64_t c = g->edges[e].count * cov1 / (cov1 + cov2);
	g->edges[el1].count += c;
	g->edges[er_rc1].count += c;
	g->edges[el2].count += (g->edges[e].count - c);
	g->edges[er_rc2].count += (g->edges[e].count - c);

	asm_remove_edge(g, e);
	asm_remove_edge(g, e_rc);
	asm_remove_edge(g, er1);
	asm_remove_edge(g, er2);
	asm_remove_edge(g, el_rc1);
	asm_remove_edge(g, el_rc2);

	return 1;
}

void resolve_bridge(struct asm_graph_t *g)
{
	gint_t e;
	int cnt = 0;
	for (e = 0; e < g->n_e; ++e) {
		gint_t e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		cnt += test_bridge(g, e);
	}
	fprintf(stderr, "Number of possible bridge to resolve: %d\n", cnt);
}
