#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "io_utils.h"
#include "khash.h"
#include "resolve.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

KHASH_SET_INIT_INT64(gint);

#define __positive_ratio(r)		((r) + (1e-6) >= 0.05)
#define LEN_VAR				1
#define MIN_SCAFFOLD_LEN			3000
#define MAX_EDGE_COUNT			1000

static inline gint_t find_adj_idx(gint_t *adj, gint_t deg, gint_t id)
{
	gint_t i, ret;
	ret = -1;
	for (i = 0; i < deg; ++i) {
		if (adj[i] == id)
			ret = i;
	}
	return ret;
}

void asm_condense2(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	gint_t *node_id;
	gint_t n_v, n_e, m_e;
	node_id = malloc(g0->n_v * sizeof(gint_t));
	memset(node_id, 255, g0->n_v * sizeof(gint_t));

	/* remove unused links */
	gint_t u;
	for (u = 0; u < g0->n_v; ++u) {
		gint_t c, deg;
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
	/* nodes on new graph only consist of branching nodes on old graph */
	n_v = 0;
	for (u = 0; u < g0->n_v; ++u) {
		gint_t deg_fw = g0->nodes[u].deg;
		gint_t deg_rv = g0->nodes[g0->nodes[u].rc_id].deg;
		/* non-branching node */
		if ((deg_fw == 1 && deg_rv == 1) || deg_fw + deg_rv == 0)
			continue;
		fprintf(stdout, "NODE: [%ld] -> [%ld]\n", u, n_v);
		node_id[u] = n_v++;
	}
	struct asm_node_t *nodes = calloc(n_v, sizeof(struct asm_node_t));
	/* set reverse complement link between nodes */
	for (u = 0; u < g0->n_v; ++u) {
		gint_t x, x_rc, u_rc;
		x = node_id[u];
		if (x == -1)
			continue;
		u_rc = g0->nodes[u].rc_id;
		x_rc = node_id[u_rc];
		assert(x_rc != -1);
		nodes[x].rc_id = x_rc;
		nodes[x_rc].rc_id = x;
		nodes[x].adj = NULL;
		nodes[x].deg = 0;
	}
	n_e = 0;
	m_e = 0x10000;
	struct asm_edge_t *edges = calloc(m_e, sizeof(struct asm_edge_t));
	/* construct new edges */
	for (u = 0; u < g0->n_v; ++u) {
		gint_t x, y_rc;
		x = node_id[u];
		if (x == -1)
			continue;
		gint_t c;
		for (c = 0; c < g0->nodes[u].deg; ++c) {
			if (n_e + 2 > m_e) {
				edges = realloc(edges, (m_e << 1) * sizeof(struct asm_edge_t));
				memset(edges + m_e, 0, m_e * sizeof(struct asm_edge_t));
				m_e <<= 1;
			}
			gint_t e = g0->nodes[u].adj[c], e_rc, v, v_rc, p, q;
			if (e == -1)
				continue;
			p = n_e; q = n_e + 1;
			edges[p].rc_id = q;
			edges[q].rc_id = p;
			asm_clone_edge(edges + p, g0->edges + e);
			gint_t *ea = NULL;
			gint_t n_ea = 0;

			do {
				ea = realloc(ea, (n_ea + 1) * sizeof(gint_t));
				ea[n_ea++] = e;
				v = g0->edges[e].target;
				if (node_id[v] == -1) { /* middle node */
					assert(g0->nodes[v].deg == 1);
					e = g0->nodes[v].adj[0];
					assert(e != -1);
					asm_append_edge_seq(edges + p,
						g0->edges + e, g0->ksize);
					edges[p].count += g0->edges[e].count;
				} else {
					break;
				}
			} while (1);
			edges[p].source = x;
			edges[p].target = node_id[v];
			asm_clone_reverse(edges + q, edges + p);
			v_rc = g0->nodes[v].rc_id;
			e_rc = g0->edges[e].rc_id;
			gint_t j = find_adj_idx(g0->nodes[v_rc].adj,
						g0->nodes[v_rc].deg, e_rc);
			assert(j >= 0);
			g0->nodes[v_rc].adj[j] = -1;
			y_rc = node_id[v_rc];
			edges[q].source = y_rc;
			edges[q].target = nodes[x].rc_id;

			nodes[x].adj = realloc(nodes[x].adj, (nodes[x].deg + 1) * sizeof(gint_t));
			nodes[x].adj[nodes[x].deg++] = p;
			nodes[y_rc].adj = realloc(nodes[y_rc].adj,
					(nodes[y_rc].deg + 1) * sizeof(gint_t));
			nodes[y_rc].adj[nodes[y_rc].deg++] = q;
			n_e += 2;
			fprintf(stdout, "EDGE: [%ld] -> ", p);
			for (j = 0; j < n_ea; ++j)
				fprintf(stdout, j + 1 == n_ea ?
						"[%ld]\n" : "[%ld], ", ea[j]);
			fprintf(stdout, "EDGE: [%ld] -> ", q);
			for (j = 0; j < n_ea; ++j)
				fprintf(stdout, j + 1 == n_ea ?
						"[%ld]\n" : "[%ld], ",
						g0->edges[ea[j]].rc_id);
			free(ea);
		}
	}
	free(node_id);
	edges = realloc(edges, n_e * sizeof(struct asm_edge_t));
	g->ksize = g0->ksize;
	g->n_v = n_v;
	g->n_e = n_e;
	g->nodes = nodes;
	g->edges = edges;
}

int is_dead_end(struct asm_graph_t *g, gint_t u)
{
	gint_t e, u_rc, v, v_rc;
	u_rc = g->nodes[u].rc_id;
	if (g->nodes[u].deg + g->nodes[u_rc].deg != 1)
		return 0;
	if (g->nodes[u].deg)
		e = g->nodes[u].adj[0];
	else
		e = g->nodes[u_rc].adj[0];
	v = g->edges[e].target;
	v_rc = g->nodes[v].rc_id;
	if (g->nodes[v].deg + g->nodes[v_rc].deg != 1)
		return 0;
	uint32_t len = get_edge_len(g->edges + e);
	return len >= 250 ? 0 : 1;
}

void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	gint_t *node_id;
	gint_t n_v, n_e, m_e;
	node_id = malloc(g0->n_v * sizeof(gint_t));
	memset(node_id, 255, g0->n_v * sizeof(gint_t));

	/* remove unused links */
	gint_t u;
	for (u = 0; u < g0->n_v; ++u) {
		gint_t c, deg;
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
	/* nodes on new graph only consist of branching nodes on old graph */
	n_v = 0;
	for (u = 0; u < g0->n_v; ++u) {
		gint_t deg_fw = g0->nodes[u].deg;
		gint_t deg_rv = g0->nodes[g0->nodes[u].rc_id].deg;
		/* non-branching node */
		if ((deg_fw == 1 && deg_rv == 1) || deg_fw + deg_rv == 0 || is_dead_end(g0, u))
			continue;
		node_id[u] = n_v++;
	}
	struct asm_node_t *nodes = calloc(n_v, sizeof(struct asm_node_t));
	/* set reverse complement link between nodes */
	for (u = 0; u < g0->n_v; ++u) {
		gint_t x, x_rc, u_rc;
		x = node_id[u];
		if (x == -1)
			continue;
		u_rc = g0->nodes[u].rc_id;
		x_rc = node_id[u_rc];
		assert(x_rc != -1);
		nodes[x].rc_id = x_rc;
		nodes[x_rc].rc_id = x;
		nodes[x].adj = NULL;
		nodes[x].deg = 0;
	}
	n_e = 0;
	m_e = 0x10000;
	struct asm_edge_t *edges = calloc(m_e, sizeof(struct asm_edge_t));
	/* construct new edges */
	for (u = 0; u < g0->n_v; ++u) {
		gint_t x, y_rc;
		x = node_id[u];
		if (x == -1)
			continue;
		gint_t c;
		for (c = 0; c < g0->nodes[u].deg; ++c) {
			if (n_e + 2 > m_e) {
				edges = realloc(edges, (m_e << 1) * sizeof(struct asm_edge_t));
				memset(edges + m_e, 0, m_e * sizeof(struct asm_edge_t));
				m_e <<= 1;
			}
			gint_t e = g0->nodes[u].adj[c], e_rc, v, v_rc, p, q;
			if (e == -1)
				continue;
			p = n_e; q = n_e + 1;
			edges[p].rc_id = q;
			edges[q].rc_id = p;
			asm_clone_edge(edges + p, g0->edges + e);

			do {
				v = g0->edges[e].target;
				if (node_id[v] == -1) { /* middle node */
					assert(g0->nodes[v].deg == 1);
					e = g0->nodes[v].adj[0];
					assert(e != -1);
					asm_append_edge_seq(edges + p,
						g0->edges + e, g0->ksize);
					edges[p].count += g0->edges[e].count;
				} else {
					break;
				}
			} while (1);
			edges[p].source = x;
			edges[p].target = node_id[v];
			asm_clone_reverse(edges + q, edges + p);
			v_rc = g0->nodes[v].rc_id;
			e_rc = g0->edges[e].rc_id;
			gint_t j = find_adj_idx(g0->nodes[v_rc].adj,
						g0->nodes[v_rc].deg, e_rc);
			assert(j >= 0);
			g0->nodes[v_rc].adj[j] = -1;
			y_rc = node_id[v_rc];
			edges[q].source = y_rc;
			edges[q].target = nodes[x].rc_id;

			nodes[x].adj = realloc(nodes[x].adj, (nodes[x].deg + 1) * sizeof(gint_t));
			nodes[x].adj[nodes[x].deg++] = p;
			nodes[y_rc].adj = realloc(nodes[y_rc].adj,
					(nodes[y_rc].deg + 1) * sizeof(gint_t));
			nodes[y_rc].adj[nodes[y_rc].deg++] = q;
			n_e += 2;
		}
	}
	free(node_id);
	edges = realloc(edges, n_e * sizeof(struct asm_edge_t));
	g->ksize = g0->ksize;
	g->n_v = n_v;
	g->n_e = n_e;
	g->nodes = nodes;
	g->edges = edges;
}

static inline int check_self_loop(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t e_rc, u, v, u_rc, v_rc;
	e_rc = g->edges[e].rc_id;
	u = g->edges[e].source;
	v = g->edges[e].target;
	u_rc = g->nodes[u].rc_id;
	v_rc = g->nodes[v].rc_id;
	if (u == v) {
		__VERBOSE("Self loop at edge %ld\n", e);
		return 1;
	}
	if (u == g->nodes[v].rc_id) {
		__VERBOSE("Self loop reverse at edge %ld\n", e);
		return 1;
	}
	return 0;
}

static inline int check_double_loop(struct asm_graph_t *g, gint_t ef, double uni_cov)
{
	gint_t ef_rc, u, v, u_rc, v_rc, er, er_rc, e1, e2;
	ef_rc = g->edges[ef].rc_id;
	u = g->edges[ef].source;
	v = g->edges[ef].target;
	u_rc = g->nodes[u].rc_id;
	v_rc = g->nodes[v].rc_id;
	/* Check topology */
	if (g->nodes[u].deg != 1 || g->nodes[v].deg != 2 ||
		g->nodes[u_rc].deg != 2 || g->nodes[v_rc].deg != 1)
		return 0;
	/* Check return edge */
	if (g->edges[g->nodes[v].adj[0]].target == u) {
		er = g->nodes[v].adj[0];
		e2 = g->nodes[v].adj[1];
	} else if (g->edges[g->nodes[v].adj[1]].target == u) {
		er = g->nodes[v].adj[1];
		e2 = g->nodes[v].adj[0];
	} else {
		return 0;
	}
	if (g->edges[g->nodes[u_rc].adj[0]].target == v_rc) {
		er_rc = g->nodes[u_rc].adj[0];
		e1 = g->nodes[u_rc].adj[1];
	} else if (g->edges[g->nodes[u_rc].adj[1]].target == v_rc) {
		er_rc = g->nodes[u_rc].adj[1];
		e1 = g->nodes[u_rc].adj[0];
	} else {
		return 0;
	}
	double cov1, cov2, cov_f, cov_r;
	cov1 = __get_edge_cov(g->edges + e1, g->ksize);
	cov2 = __get_edge_cov(g->edges + e2, g->ksize);
	cov_f = __get_edge_cov(g->edges + ef, g->ksize);
	cov_r = __get_edge_cov(g->edges + er, g->ksize);
	__VERBOSE("Loop: %ld-(%ld-%ld)-%ld\n%.6lf-(%.6lf-%.6lf)-%.6lf\n",
		e1, ef, er, e2, cov1, cov_f, cov_r, cov2);
	// asm_append_edge_seq(g->edges + e_return,
	// 			g->edges + e, g->ksize);
	// asm_append_edge_seq(g->edges + e,
	// 			g->edges + e_return, g->ksize);
	// g->edges[e].count += g->edges[e_return].count;
	// asm_append_edge_seq(g->edges + e_return_rc,
	// 			g->edges + e_rc, g->ksize);
	// asm_append_edge_seq(g->edges + e_rc,
	// 		g->edges + e_return_rc, g->ksize);
	// asm_remove_edge(g, e_return);
	// asm_remove_edge(g, e_return_rc);
	// return 1;

	return 1;
}

/* return 0: not at all
 * return 1: self - loop
 * return 2: double - loop
 * return 3: false loop, just go through
 */
static inline int check_simple_loop(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	// int ret = check_self_loop(g, e, uni_cov);
	// if (ret)
	// 	return ret;
	// ret = check_double_loop(g, e, uni_cov);
	// if (ret)
	// 	return ret;
	// return 0;
	gint_t u, v, u_rc, v_rc, e_rc, e_left, e_right, e_return, e_return_rc;
	e_rc = g->edges[e].rc_id;
	u = g->edges[e].source;
	v = g->edges[e].target;
	u_rc = g->nodes[u].rc_id;
	v_rc = g->nodes[v].rc_id;
	if (g->nodes[u].deg != 1 || g->nodes[v].deg != 2 ||
		g->nodes[u_rc].deg != 2 || g->nodes[v_rc].deg != 1)
		return 0;
	if (g->edges[g->nodes[v].adj[0]].target == u) {
		e_return = g->nodes[v].adj[0];
		e_right = g->nodes[v].adj[1];
	} else if (g->edges[g->nodes[v].adj[1]].target == u) {
		e_return = g->nodes[v].adj[1];
		e_right = g->nodes[v].adj[0];
	} else {
		return 0;
	}
	if (g->edges[g->nodes[u_rc].adj[0]].target == v_rc) {
		e_return_rc = g->nodes[u_rc].adj[0];
		e_left = g->nodes[u_rc].adj[1];
	} else if (g->edges[g->nodes[u_rc].adj[1]].target == v_rc) {
		e_return_rc = g->nodes[u_rc].adj[1];
		e_left = g->nodes[u_rc].adj[0];
	} else {
		return 0;
	}
	assert(g->edges[e_return].rc_id == e_return_rc);
	// float cov_left, cov_right, cov_ahead, cov_return;
	// cov_left = __get_edge_cov(g->edges + e_left, g->ksize);
	// cov_right = __get_edge_cov(g->edges + e_right, g->ksize);
	// cov_ahead = __get_edge_cov(g->edges + e, g->ksize);
	// cov_return = __get_edge_cov(g->edges + e_return, g->ksize);

	// if (__int_cov_ratio(cov_left, cov_right) != 1)
	// 	return 0;
	// if (__int_cov_ratio(cov_ahead, cov_return) < 2) {
		// if (__int_cov_ratio(cov_left, cov_ahead) == 1) {
		// 	asm_remove_edge(g, e_return);
		// 	asm_remove_edge(g, e_return_rc);
		// 	/* remove the return edge */
		// 	return 2;
		// } else {
		// 	/* still resolve */
		// 	asm_append_edge_seq(g->edges + e_return,
		// 				g->edges + e, g->ksize);
		// 	asm_append_edge_seq(g->edges + e,
		// 				g->edges + e_return, g->ksize);
		// 	g->edges[e].count += g->edges[e_return].count;
		// 	asm_append_edge_seq(g->edges + e_return_rc,
		// 				g->edges + e_rc, g->ksize);
		// 	asm_append_edge_seq(g->edges + e_rc,
		// 			g->edges + e_return_rc, g->ksize);
		// 	asm_remove_edge(g, e_return);
		// 	asm_remove_edge(g, e_return_rc);
		// 	return 1;
		// }
	// } else {
		/* resolve */
		asm_append_edge_seq(g->edges + e_return,
					g->edges + e, g->ksize);
		asm_append_edge_seq(g->edges + e,
					g->edges + e_return, g->ksize);
		g->edges[e].count += g->edges[e_return].count;
		asm_append_edge_seq(g->edges + e_return_rc,
					g->edges + e_rc, g->ksize);
		asm_append_edge_seq(g->edges + e_rc,
				g->edges + e_return_rc, g->ksize);
		g->edges[e_rc].count += g->edges[e_return_rc].count;
		asm_remove_edge(g, e_return);
		asm_remove_edge(g, e_return_rc);
		return 1;
	// }
}

int unroll_simple_loop(struct asm_graph_t *g, double uni_cov)
{
	/*
	 *                    1g
	 *              +-------------+
	 *             /               \
	 *      1g    v       2g        \       1g
	 * o--------->o=================>o-------------->o
	 *            u                  v
	 */
	gint_t e, cnt_loop, cnt_removed;
	int ret;
	cnt_loop = cnt_removed = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		ret = check_simple_loop(g, e, uni_cov);
		if (ret == 1)
			++cnt_loop;
		else if (ret == 2)
			++cnt_removed;
	}
	__VERBOSE("Number of unroll loops: %ld\n", cnt_loop);
	__VERBOSE("Number of false loops: %ld\n", cnt_removed);
	return cnt_loop + cnt_removed;
}

void remove_self_loop(struct asm_graph_t *g)
{
	/* Loop is bad
	 *               +---+
	 *              /    |
	 *             +     +  <---------- loop
	 *             |    /
	 * o==========>o<--+
	 *             \\
	 *               ++========>o
	 */
	gint_t e, e_rc, u, u_rc;
	gint_t cnt = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		e_rc = g->edges[e].rc_id;
		if (e != e_rc)
			continue;
		u = g->edges[e].source;
		u_rc = g->nodes[e].rc_id;
		assert(g->edges[e].target == u_rc);
		if (g->nodes[u_rc].deg == 1 && g->nodes[u].deg == 2) {
			gint_t e_true, e_prev;
			if (g->nodes[u].adj[0] == e)
				e_true = g->nodes[u].adj[1];
			else
				e_true = g->nodes[u].adj[0];
			e_prev = g->nodes[u_rc].adj[0];
			double cov1, cov2;
			cov1 = g->edges[e_prev].count / (g->edges[e_prev].seq_len - g->ksize);
			cov2 = g->edges[e_true].count / (g->edges[e_true].seq_len - g->ksize);
			if (__int_cov_ratio(cov1, cov2) != 1)
				continue;
			g->nodes[u].adj[0] = e_true;
			g->nodes[u].deg = 1;
			g->nodes[u].adj = realloc(g->nodes[u].adj, sizeof(gint_t));
			++cnt;
		}
	}
	__VERBOSE("Number of removed loops: %ld\n", cnt);
}

void remove_bubble_simple(struct asm_graph_t *g0, double uni_cov)
{
	/* "Don't you want a balloon?"
	 * 1g = 1 genome walk
	 *                   0.8g
	 *                +--------+
	 *               /          \
	 *       1g     /            v    1g
	 * o---------->U             V--------->o
	 *              \            ^
	 *               \   0.3g   /
	 *                +--------+
	 */
	gint_t cnt = 0;
	gint_t u;
	for (u = 0; u < g0->n_v; ++u) {
		gint_t ctg;
		double cov;
		/* check topology of node u */
		gint_t u_rc = g0->nodes[u].rc_id;
		if (g0->nodes[u_rc].deg != 1 || g0->nodes[u].deg != 2)
			continue;
		/* check if previous contig is uni genome walk */
		ctg = g0->nodes[u_rc].adj[0];
		cov = g0->edges[ctg].count * 1.0 /
			(g0->edges[ctg].seq_len - g0->ksize);
		if (__int_cov_ratio(cov, uni_cov) != 1)
			continue;

		/* check topology of 2 small contigs (end at same node) */
		gint_t e0, e1;
		e0 = g0->nodes[u].adj[0];
		e1 = g0->nodes[u].adj[1];
		if (e0 == g0->edges[e1].rc_id)
			continue;
		if (g0->edges[e0].target != g0->edges[e1].target)
			continue;
		/* check topology of node v */
		gint_t v, v_rc;
		v = g0->edges[e0].target;
		v_rc = g0->nodes[v].rc_id;
		if (g0->nodes[v_rc].deg != 2 || g0->nodes[v].deg != 1)
			continue;
		/* check if next contig is uni genome walk */
		ctg = g0->nodes[v].adj[0];
		cov = g0->edges[ctg].count * 1.0 /
			(g0->edges[ctg].seq_len - g0->ksize);
		if (__int_cov_ratio(cov, uni_cov) != 1)
			continue;

		gint_t e_rc0, e_rc1;
		e_rc0 = g0->nodes[v_rc].adj[0];
		e_rc1 = g0->nodes[v_rc].adj[1];
		assert(g0->edges[e_rc0].target == g0->edges[e_rc1].target);
		/* keep the longer edge */
		gint_t e;
		if (g0->edges[e0].seq_len > g0->edges[e1].seq_len)
			e = e0;
		else
			e = e1;
		g0->nodes[u].adj[0] = e;
		g0->nodes[u].deg = 1;
		g0->nodes[u].adj = realloc(g0->nodes[u].adj, sizeof(gint_t));
		g0->edges[e].count = g0->edges[e0].count + g0->edges[e1].count;

		assert(g0->edges[e].rc_id == e_rc0 || g0->edges[e].rc_id == e_rc1);
		g0->nodes[v_rc].adj[0] = g0->edges[e].rc_id;
		g0->nodes[v_rc].deg = 1;
		g0->nodes[v_rc].adj = realloc(g0->nodes[v_rc].adj, sizeof(gint_t));
		g0->edges[g0->edges[e].rc_id].count = g0->edges[e_rc0].count + g0->edges[e_rc1].count;
		++cnt;
	}
	__VERBOSE("Number of collapsed bubble: %ld\n", cnt);
}

static void print_debug_info(struct asm_graph_t *g)
{
	fprintf(stderr, "Edge 51446 info:\n");
	gint_t u, v;
	u = g->edges[51446].source;
	v = g->edges[51446].target;
	fprintf(stderr, "Source: %ld; Target: %ld\n", u, v);
	fprintf(stderr, "src_rc: %ld; tgt_rc: %ld\n", g->nodes[u].rc_id,
		g->nodes[v].rc_id);
	fprintf(stderr, "deg[u] = %ld\n", g->nodes[u].deg);
	fprintf(stderr, "deg[v] = %ld\n", g->nodes[v].deg);
	fprintf(stderr, "deg[u_rc] = %ld\n", g->nodes[g->nodes[u].rc_id].deg);
	fprintf(stderr, "deg[v_rc] = %ld\n", g->nodes[g->nodes[v].rc_id].deg);
}

void remove_bubble_and_loop(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g0);
	__VERBOSE("1 genome walk coverage: %lf\n", uni_cov);
	unroll_simple_loop(g0, uni_cov);
	// remove_self_loop(g0);
	remove_bubble_simple(g0, uni_cov);
	asm_condense2(g0, g);
	// print_debug_info(g);
}

static int dfs_dead_end(struct asm_graph_t *g, gint_t u,
			gint_t len, gint_t max_len, int ksize)
{
	if (len >= max_len)
		return 0;
	gint_t j, e;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		int k = dfs_dead_end(g, g->edges[e].target,
			len + g->edges[e].seq_len - ksize, max_len, ksize);
		if (k == 0)
			return 0;
	}
	return 1;
}

struct edge_cov_t {
	double cov;
	gint_t e;
};

static inline void ec_sort(struct edge_cov_t *b, struct edge_cov_t *e)
{
	struct edge_cov_t *i, *j, t;
	for (i = b + 1; i < e; ++i) {
		if (i->cov < (i - 1)->cov) {
			t = *i;
			for (j = i; j > b && t.cov < (j - 1)->cov; --j)
				*j = *(j - 1);
			*j = t;
		}
	}
}

void remove_tips_topology(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	/*             o
	 *             ^      o
	 *       o<--+ |      ^
	 *            \|     /         TIP
	 *        o<---o--->o--->o
	 *             ^
	 *            /
	 * o-------->o-------->o------>o------>o
	 */
	gint_t u, e;
	uint32_t *flag = calloc((g0->n_e + 31) >> 5, sizeof(uint32_t));
	struct edge_cov_t *tmp = NULL;
	gint_t tmp_len = 0;
	gint_t count_rm = 0;
	for (u = 0; u < g0->n_v; ++u) {
		if (u % 1000 == 0)
			fprintf(stderr, "\rRemoving dead-end node %ld", u);
		gint_t c, deg, cnt;
		deg = g0->nodes[u].deg;
		if (deg > tmp_len) {
			tmp_len = deg;
			tmp = realloc(tmp, tmp_len * sizeof(struct edge_cov_t));
		}
		cnt = 0;
		for (c = 0; c < g0->nodes[u].deg; ++c) {
			gint_t e_id;
			e_id = g0->nodes[u].adj[c];
			int ret = dfs_dead_end(g0, g0->edges[e_id].target,
					g0->ksize, MAX_TIPS_LEN, g0->ksize);
			if (ret) {
				tmp[cnt].e = e_id;
				tmp[cnt].cov = g0->edges[e_id].count * 1.0 /
					(g0->edges[e_id].seq_len - g0->ksize);
				++cnt;
			}
		}
		if (!cnt)
			continue;
		ec_sort(tmp, tmp + cnt);
		for (c = 0; c < cnt && deg > 1; ++c) {
			gint_t e_id, e_rc;
			e_id = tmp[c].e;
			e_rc = g0->edges[e_id].rc_id;
			flag[e_id >> 5] |= (uint32_t)1 << (e_id & 31);
			flag[e_rc >> 5] |= (uint32_t)1 << (e_rc & 31);
			--deg;
			++count_rm;
		}
	}
	fprintf(stderr, "\nNumber of removed edges: %ld\n", count_rm);
	free(tmp);
	for (e = 0; e < g0->n_e; ++e) {
		if (!((flag[e >> 5] >> (e & 31)) & 1))
			continue;
		u = g0->edges[e].source;
		gint_t c = find_adj_idx(g0->nodes[u].adj, g0->nodes[u].deg, e);
		assert(c != -1);
		g0->nodes[u].adj[c] = -1;
	}
	for (u = 0; u < g0->n_v; ++u) {
		gint_t c, deg;
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
	free(flag);

	asm_condense(g0, g);
}

void remove_tips(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	/*
	 *                      o
	 *                      ^
	 *                     /
	 *              cov~1 /
	 *                   /
	 *       cov~10     /   cov~10
	 * o-------------->o------------->o---------->o
	 */
	gint_t u;
	for (u = 0; u < g0->n_v; ++u) {
		double max_cov = 0.0, cov;
		gint_t deg = g0->nodes[u].deg, c;
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
			if (cov / max_cov < TIPS_RATIO_THRESHOLD) {
				gint_t v, v_rc;
				v = g0->edges[e_id].target;
				v_rc = g0->nodes[v].rc_id;
				e_rc = g0->edges[e_id].rc_id;
				gint_t j = find_adj_idx(g0->nodes[v_rc].adj,
						g0->nodes[v_rc].deg, e_rc);
				assert(j != -1);
				/* disconnect edge */
				g0->nodes[u].adj[c] = -1;
				g0->nodes[v_rc].adj[j] = -1;
			}
		}
	}

	asm_condense(g0, g);
}

int test_split(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t u, v, u_rc, v_rc, e_rc;
	u = g->edges[e].source;
	v = g->edges[e].target;
	e_rc = g->edges[e].rc_id;
	u_rc = g->nodes[u].rc_id;
	v_rc = g->nodes[v].rc_id;
	if (g->nodes[u_rc].deg < 1 || g->nodes[u].deg != 1 || g->nodes[v].deg != 1)
		return 0;
	// double cov, sum_fcov;
	// int e_cov, sum_cov;
	// cov = __get_edge_cov(g->edges + e, g->ksize);
	// e_cov = (int)(cov / uni_cov + 0.499999999);
	// sum_cov = 0;
	// sum_fcov = 0.0;
	double e_cov, sum_cov, cov;
	e_cov = __get_edge_cov(g->edges + e, g->ksize);
	// fprintf(stderr, "consider edge %ld: cov = %.3lf\n", e, e_cov);
	sum_cov = 0.0;
	gint_t j;
	uint32_t e_len, max_len;
	e_len = get_edge_len(g->edges + e);
	max_len = 0;
	for (j = 0; j < g->nodes[u_rc].deg; ++j) {
		gint_t n, n_rc;
		n_rc = g->nodes[u_rc].adj[j];
		n = g->edges[n_rc].rc_id;
		cov = __get_edge_cov(g->edges + n, g->ksize);
		// fprintf(stderr, "\tsattelite: edge %ld ~ %.3lf\n", n, cov);
		sum_cov += cov;
		// cov = __get_edge_cov(g->edges + n, g->ksize);
		// sum_fcov += cov;
		// sum_cov += (int)(cov / uni_cov + 0.499999999);
		uint32_t len = get_edge_len(g->edges + n);
		max_len = __max(max_len, len);
	}
	// fprintf(stderr, "consider edge %ld: cov = %d, sum sattelite cov = %d\n",
	// 	e, e_cov, sum_cov);
	if (e_len > 500 || max_len > 500) {
		if (e_cov < sum_cov * 0.75)
			return -1;
	}
	for (j = 0; j < g->nodes[u_rc].deg; ++j) {
		gint_t n, n_rc;
		n_rc = g->nodes[u_rc].adj[j];
		n = g->nodes[n_rc].rc_id;
		if (n == e || n_rc == e)
			return -1;
	}
	g->nodes[v_rc].adj = realloc(g->nodes[v_rc].adj,
		(g->nodes[v_rc].deg + g->nodes[u_rc].deg) * sizeof(gint_t));
	for (j = 0; j < g->nodes[u_rc].deg; ++j) {
		gint_t n, n_rc;
		n_rc = g->nodes[u_rc].adj[j];
		n = g->edges[n_rc].rc_id;
		// __VERBOSE("Joining edge: %ld -> %ld\n", n, e);
		// fprintf(stderr, "concat edge: (%ld %ld) -> (%ld %ld); (%ld %ld) -> (%ld %ld)\n",
		// 	g->edges[n].source, g->edges[n].target,
		// 	g->edges[e].source, g->edges[e].target,
		// 	g->edges[e_rc].source, g->edges[e_rc].target,
		// 	g->edges[n_rc].source, g->edges[n_rc].target);
		cov = __get_edge_cov(g->edges + n, g->ksize);
		uint64_t split_count = (uint64_t)(cov / sum_cov * g->edges[e].count);

		asm_append_edge_seq(g->edges + n, g->edges + e, g->ksize);
		g->edges[n].count += split_count;
		g->edges[n].target = v;

		asm_clean_edge_seq(g->edges + n_rc);
		asm_clone_reverse(g->edges + n_rc, g->edges + n);
		g->edges[n_rc].source = v_rc;
		g->edges[n_rc].target = g->nodes[g->edges[n].source].rc_id;
		g->nodes[v_rc].adj[g->nodes[v_rc].deg++] = n_rc;
		g->edges[n].rc_id = n_rc;
		g->edges[n_rc].rc_id = n;
		// if (get_edge_len(g->edges + n) >= 5000 && __get_edge_cov(g->edges + n, g->ksize) < 100.0)
		// 	__VERBOSE("old_cov = %.9lf; new_cov = %.9lf\n", cov, __get_edge_cov(g->edges + n, g->ksize));
	}
	asm_remove_edge(g, e);
	asm_remove_edge(g, e_rc);
	free(g->nodes[u_rc].adj);
	g->nodes[u_rc].adj = NULL;
	g->nodes[u_rc].deg = 0;
	return 1;
}

int graph_expanding(struct asm_graph_t *g, double uni_cov)
{
	int cnt, cnt_fp, step;
	cnt = cnt_fp = step =  0;
	gint_t e;
	while (1) {
		int local_cnt = 0;
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			int ret = test_split(g, e, uni_cov);
			if (ret == -1)
				++cnt_fp;
			else if (ret == 1)
				++local_cnt;
		}
		cnt += local_cnt;
		if (local_cnt == 0)
			break;
	}
	__VERBOSE("Number of expanding edges: %d\n", cnt);
	__VERBOSE("Number of edges cannot expand: %d\n", cnt_fp);
	return cnt;
}

int test_bubble2(struct asm_graph_t *g, gint_t u)
{
	if (g->nodes[u].deg < 2)
		return 0;
	gint_t j, k, e, ke, v, deg;
	deg = 0;
	int count_resolve = 0;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		if (e == -1)
			continue;
		v = g->edges[e].target;
		float cov, cur_cov = 0.0;
		uint32_t max_len = 0;
		uint32_t base_len = get_edge_len(g->edges + e);
		gint_t idx = -1;
		int cnt = 0;
		for (k = 0; k < g->nodes[u].deg; ++k) {
			e = g->nodes[u].adj[k];
			if (e == -1 || v != g->edges[e].target)
				continue;
			uint32_t len = get_edge_len(g->edges + e);
			if (!((len < 500 && base_len < 500) ||
				(base_len - LEN_VAR <= len && base_len + LEN_VAR >= len)))
				continue;
			cov = __get_edge_cov(g->edges + e, g->ksize);
			if (cov > cur_cov) {
				cur_cov = cov;
				idx = k;
			}
			max_len = __max(max_len, len);
			++cnt;
		}
		if (cnt < 2 || v == g->nodes[u].rc_id) {
			e = g->nodes[u].adj[j];
			g->nodes[u].adj[j] = -1;
			g->nodes[u].adj[deg++] = e;
			continue;
		}
		assert(v != u);
		ke = g->nodes[u].adj[idx];
		uint64_t kmer_count = 0;
		for (k = 0; k < g->nodes[u].deg; ++k) {
			e = g->nodes[u].adj[k];
			if (e == -1 || v != g->edges[e].target)
				continue;
			uint32_t len = get_edge_len(g->edges + e);
			if (!((len < 500 && base_len < 500) ||
				(base_len - LEN_VAR <= len && base_len + LEN_VAR >= len)))
				continue;
			kmer_count += g->edges[e].count;
			if (e != ke) {
				g->edges[e].source = g->edges[e].target = -1;
				asm_remove_edge(g, g->edges[e].rc_id);
				g->nodes[u].adj[k] = -1;
			}
		}
		g->edges[ke].count = kmer_count;
		g->edges[g->edges[ke].rc_id].count = kmer_count;
		g->nodes[u].adj[idx] = -1;
		g->nodes[u].adj[deg++] = ke;
		++count_resolve;
	}
	g->nodes[u].adj = realloc(g->nodes[u].adj, deg * sizeof(gint_t));
	g->nodes[u].deg = deg;
	return count_resolve;
}

int resolve_bubble2(struct asm_graph_t *g)
{
	gint_t v;
	int cnt = 0;
	for (v = 0; v < g->n_v; ++v) {
		int ret = test_bubble2(g, v);
		cnt += ret;
	}
	__VERBOSE("Number of collapse bubble: %d\n", cnt);
	return cnt;
}

void resolve_chain(struct asm_graph_t *g0, struct asm_graph_t *g1)
{
	int step = 0;
	while (1) {
		__VERBOSE("Iteration [%d]\n", step++);
		int cnt_expand, cnt_collapse, cnt_loop;
		double uni_cov = get_genome_coverage(g0);
		__VERBOSE("Genome coverage: %.9lf\n", uni_cov);
		cnt_loop = unroll_simple_loop(g0, uni_cov);
		cnt_collapse = resolve_bubble2(g0);
		cnt_expand = graph_expanding(g0, uni_cov);
		asm_condense(g0, g1);
		test_asm_graph(g1);
		if (cnt_expand == 0 && cnt_collapse == 0 && cnt_loop == 0)
			break;
		*g0 = *g1;
	}
}

static inline int cb_add(khash_t(gint) *h, gint_t k)
{
	int ret;
	kh_put(gint, h, k, &ret);
	return ret;
}

void find_region(struct asm_graph_t *g, gint_t se, uint32_t min_contig_len,
		uint32_t max_edge_count, khash_t(gint) *set_v, khash_t(gint) *set_e)
{
	gint_t *q = malloc(max_edge_count * sizeof(gint_t));
	int l, r;
	l = r = 0;
	uint32_t edge_count = 1;
	q[0] = se;
	cb_add(set_e, se);
	cb_add(set_e, g->edges[se].rc_id);
	cb_add(set_v, g->edges[se].target);
	while (l <= r) {
		gint_t e, u, u_rc, v, j;
		e = q[l++];
		v = g->edges[e].target;
		for (j = 0; j < g->nodes[v].deg; ++j) {
			gint_t ne = g->nodes[v].adj[j], ne_rc;
			ne_rc = g->edges[ne].rc_id;
			uint32_t len = get_edge_len(g->edges + ne);
			if (len < min_contig_len) {
				if (cb_add(set_e, ne) == 1) {
					if (edge_count == max_edge_count)
						goto clean_up;
					q[++r] = ne;
				}
				if (cb_add(set_e, ne_rc) == 1) {
					if (edge_count == max_edge_count)
						goto clean_up;
					q[++r] = ne_rc;
				}
				cb_add(set_v, g->edges[ne].target);
				cb_add(set_v, g->edges[ne_rc].target);
				cb_add(set_v, g->nodes[g->edges[ne].target].rc_id);
				cb_add(set_v, g->nodes[g->edges[ne_rc].target].rc_id);
			} else {
				cb_add(set_e, ne);
				cb_add(set_e, ne_rc);
				// cb_add(set_v, g->edges[ne].target);
				// cb_add(set_v, g->edges[ne_rc].target);
			}
		}
	}
clean_up:
	free(q);
}

void kh_merge_set(khash_t(gint) *dst, khash_t(gint) *src)
{
	khiter_t k;
	int ret;
	for (k = kh_begin(src); k != kh_end(src); ++k) {
		if (!kh_exist(src, k))
			continue;
		kh_put(gint, dst, kh_key(src, k), &ret);
	}
}

void print_set_info(khash_t(gint) *set_e)
{
	khiter_t k;
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (kh_exist(set_e, k))
			__VERBOSE("%ld\n", kh_key(set_e, k));
	}
}

static int is_hang_edge(struct asm_graph_t *g, khash_t(gint) *set_v,
				khash_t(gint) *set_e, gint_t e)
{
	gint_t v = g->edges[e].target;
	if (kh_get(gint, set_v, v) == kh_end(set_v))
		return 1;
	else
		return 0;
}

void detect_leg(struct asm_graph_t *g, khash_t(gint) *set_v,
			khash_t(gint) *set_e, khash_t(gint) *set_leg)
{
	khiter_t k;
	gint_t e;
	int ret;
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (!kh_exist(set_e, k))
			continue;
		e = kh_key(set_e, k);
		if (is_hang_edge(g, set_v, set_e, e))
			kh_put(gint, set_leg, e, &ret);
	}
}

int dfs_check_path(struct asm_graph_t *g, khash_t(gint) *set_e,
					khash_t(gint) *vis, gint_t u, gint_t t)
{
	gint_t j, e, v;
	int ret;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		if (kh_get(gint, set_e, e) == kh_end(set_e))
			continue;
		v = g->edges[e].target;
		if (v == t)
			return 1;
		if (kh_get(gint, vis, v) == kh_end(vis)) {
			kh_put(gint, vis, v, &ret);
			ret = dfs_check_path(g, set_e, vis, v, t);
			if (ret)
				return 1;
		}
	}
	return 0;
}

int check_path(struct asm_graph_t *g, khash_t(gint) *set_e, gint_t s, gint_t e)
{
	if (s == e)
		return 1;
	int ret;
	khash_t(gint) *vis;
	vis = kh_init(gint);
	kh_put(gint, vis, s, &ret);
	ret = dfs_check_path(g, set_e, vis, s, e);
	kh_destroy(gint, vis);
	return ret;
}

void collapse_1_1_complex(struct asm_graph_t *g, khash_t(gint) *set_e,
				khash_t(gint) *set_leg, double uni_cov)
{
	khiter_t k;
	gint_t legs[10];
	int n_leg, cov1, cov2;
	n_leg = 0;
	for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
		if (kh_exist(set_leg, k)) {
			legs[n_leg++] = kh_key(set_leg, k);
			kh_del(gint, set_e, kh_key(set_leg, k));
			kh_del(gint, set_e, g->edges[kh_key(set_leg, k)].rc_id);
		}
	}
	assert(n_leg == 2);
	gint_t e1, e2;
	e1 = legs[0];
	e2 = legs[1];
	cov1 = (int)(__get_edge_cov(g->edges + e1, g->ksize) / uni_cov + 0.499999);
	cov2 = (int)(__get_edge_cov(g->edges + e2, g->ksize) / uni_cov + 0.499999);
	if (cov1 != cov2) {
		__VERBOSE("WARNING: not collapse 1-1 jungle\n");
		return;
	}
	if (!check_path(g, set_e, g->nodes[g->edges[legs[0]].source].rc_id,
			g->edges[legs[1]].source)) {
		__VERBOSE("WARNING: not collapse 1-1 jungle, no path\n");
		return;
	}
	uint32_t gap_size = 0;
	uint64_t gap_count = 0;
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (!kh_exist(set_e, k))
			continue;
		gint_t e = kh_key(set_e, k);
		gint_t len = get_edge_len(g->edges + e);
		int cov = (int)(__get_edge_cov(g->edges + e, g->ksize) / uni_cov + 0.499999);
		gap_size += cov * (len - g->ksize);
		gap_count += g->edges[e].count;
		/* Remove edges , e.i isolate the nodes */
		__VERBOSE("Removing edge %ld\n", e);
		asm_remove_edge(g, e);
	}
	gap_size /= 2;
	gap_count /= 2;
	__VERBOSE("Joining edge %ld_%ld <-> %ld_%ld, gap size = %u\n",
		g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id, gap_size);
	asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
				e2, g->edges[e2].rc_id, gap_size, gap_count);
}

void detect_complex_1_1(struct asm_graph_t *g, uint32_t min_contig_size,
		uint32_t max_edge_count)
{
	double uni_cov = get_genome_coverage(g);
	__VERBOSE("Genome coverage: %.9lf\n", uni_cov);
	khash_t(gint) *visited, *set_e, *set_v, *set_leg;
	visited = kh_init(gint);
	set_e = kh_init(gint);
	set_v = kh_init(gint);
	set_leg = kh_init(gint);
	gint_t e;
	for (e = 0; e < g->n_e; ++e) {
		uint32_t len = get_edge_len(g->edges + e);
		if (kh_get(gint, visited, g->edges[e].target) != kh_end(visited) ||
			len < min_contig_size || g->edges[e].source == -1)
			continue;
		find_region(g, e, min_contig_size, max_edge_count, set_v, set_e);
		if (kh_size(set_e) < max_edge_count) {
			kh_merge_set(visited, set_v);
			detect_leg(g, set_v, set_e, set_leg);
			// __VERBOSE("Edge: %ld; len: %u; Edge count in region: %u; Number of leg: %u\n",
			// 	e, len, kh_size(set_e), kh_size(set_leg));
			if (kh_size(set_leg) == 2)
				collapse_1_1_complex(g, set_e, set_leg, uni_cov);
			// if (kh_size(set_e) < 20)
			// 	print_set_info(set_e);
		}
		kh_clear(gint, set_leg);
		kh_clear(gint, set_e);
		kh_clear(gint, set_v);
	}
}

void collapse_1_1_jungle(struct asm_graph_t *g0, struct asm_graph_t *g1)
{
	detect_complex_1_1(g0, 3000, 1000);
	asm_condense(g0, g1);
}

gint_t bc_find_best_pair(struct asm_graph_t *g, gint_t se,
						gint_t *adj, gint_t deg)
{
	gint_t j, e, se_rc, ret_e;
	se_rc = g->edges[se].rc_id;
	double ret_ratio = -1;
	ret_e = -1;
	for (j = 0; j < deg; ++j) {
		e = adj[j];
		if (e == se || e == se_rc)
			continue;
		double ratio = get_barcode_ratio(g, se, e);
		if (ratio > ret_ratio) {
			ret_ratio = ratio;
			ret_e = e;
		}
	}
	if (__positive_ratio(ret_ratio))
		return ret_e;
	return -1;
}

static inline void asm_add_node_adj(struct asm_graph_t *g, gint_t u, gint_t e)
{
	g->nodes[u].adj = realloc(g->nodes[u].adj, (g->nodes[u].deg + 1) * sizeof(gint_t));
	g->nodes[u].adj[g->nodes[u].deg++] = e;
}

void check_n_m_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t u, v, v_rc, u_rc, e1, e2, e_rc;
	e_rc = g->edges[e].rc_id;
	v = g->edges[e].target;
	v_rc = g->nodes[v].rc_id;
	u = g->edges[e].source;
	u_rc = g->nodes[u].rc_id;
	if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
		(g->nodes[u_rc].deg < 2 && g->nodes[v].deg < 2))
		return;
	int ecov, cov1, cov2;
	uint64_t e_uni_count;
	ecov = (int)(__get_edge_cov(g->edges + e, g->ksize) / uni_cov + 0.499999);
	if (ecov == 0)
		return;
	e_uni_count = g->edges[e].count / ecov;
	__VERBOSE("%ld-%ld bridge: %ld_%ld ~ %d\n", g->nodes[u_rc].deg,
		g->nodes[v].deg, e, g->edges[e].rc_id, ecov);
	int i, k;
	// for (i = 0; i < 2; ++i) {
	// 	e1 = g->nodes[u_rc].adj[i];
	// 	cov = (int)(__get_edge_cov(g->edges + e1, g->ksize) / uni_cov + 0.499999);
	// 	__VERBOSE("Leg left %i: %ld_%ld ~ %d\n",
	// 		i, e1, g->edges[e1].rc_id, cov);
	// }
	// for (i = 0; i < 2; ++i) {
	// 	e2 = g->nodes[v].adj[i];
	// 	cov = (int)(__get_edge_cov(g->edges + e2, g->ksize) / uni_cov + 0.499999);
	// 	__VERBOSE("Leg right %i: %ld_%ld ~ %d\n",
	// 		i, e2, g->edges[e2].rc_id, cov);
	// }

	int resolve;
	do {
		resolve = 0;
		for (i = 0; i < g->nodes[u_rc].deg; ++i) {
			e1 = g->nodes[u_rc].adj[i];
			cov1 = (int)(__get_edge_cov(g->edges + e1, g->ksize) / uni_cov + 0.4999999);
			if (cov1 != 1)
				continue;
			e2 = bc_find_best_pair(g, e1, g->nodes[v].adj, g->nodes[v].deg);
			if (e2 == -1)
				continue;
			// gint_t et1 = bc_find_best_pair(g, e2, g->nodes[u_rc].adj, g->nodes[u_rc].deg);
			// __VERBOSE("e1 = %ld; e2 = %ld; et1 = %ld\n", e1, e2, et1);
			// if (et1 != e1)
			// 	continue;
			cov2 = (int)(__get_edge_cov(g->edges + e2, g->ksize) / uni_cov + 0.4999999);
			if (cov2 == 1) {
				/* join edge e1.rc -> e -> e2 */
				__VERBOSE("Joining edge %ld -> %ld -> %ld\n", e1, e, e2);
				asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
					e2, g->edges[e2].rc_id, e_uni_count * cov1);
				ecov -= cov1;
			} else if (cov2 > 1) {
				g->edges = realloc(g->edges, (g->n_e + 2) * sizeof(struct asm_edge_t));
				g->n_e += 2;
				asm_clone_edge2(g, g->n_e - 2, e2);
				asm_clone_edge2(g, g->n_e - 1, g->edges[e2].rc_id);
				g->edges[g->n_e - 2].rc_id = g->n_e - 1;
				g->edges[g->n_e - 1].rc_id = g->n_e - 2;
				g->edges[g->n_e - 2].count = g->edges[g->n_e - 1].count = g->edges[e2].count / cov2;
				g->edges[e2].count = g->edges[g->edges[e2].rc_id].count = g->edges[e2].count / cov2 * (cov2 - 1);
				asm_add_node_adj(g, g->edges[g->n_e - 2].source, g->n_e - 2);
				asm_add_node_adj(g, g->edges[g->n_e - 1].source, g->n_e - 1);
				asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
					g->n_e - 2, g->n_e - 1, e_uni_count * cov1);
				ecov -= cov1;
			} else {
				/* cov2 == 0 */
				continue;
			}
			// if (cov1 != cov2 || cov1 != 1)
			// 	continue;
			// g->edges[e].count = g->edges[e_rc].count = g->edges[e].count * (ecov - cov1) / ecov;
			resolve = 1;
			break;
		}
	} while (resolve);

	/* special case when resolved all but 1 pair left */
	if (g->nodes[u_rc].deg == 1 && g->nodes[v].deg == 1) {
		e1 = g->nodes[u_rc].adj[0];
		e2 = g->nodes[v].adj[0];
		cov1 = (int)(__get_edge_cov(g->edges + e1, g->ksize) / uni_cov + 0.4999999);
		cov2 = (int)(__get_edge_cov(g->edges + e2, g->ksize) / uni_cov + 0.4999999);
		double ratio = get_barcode_ratio(g, e1, e2);
		/* FIXME: may need more complicated resolve here */
		__VERBOSE("Leftover path: %ld -> %ld -> %ld, ratio: %.6lf\n",
				e1, e, e2, ratio);
		if (__positive_ratio(ratio) ||
			(ratio < 0 && ecov >= cov1 && cov1 == cov2)) {
			__VERBOSE("Joining edge %ld -> %ld -> %ld\n", e1, e, e2);
			asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
				e2, g->edges[e2].rc_id, e_uni_count * cov1);
		}
		/* destroy the bridge */
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
	} else if (g->nodes[u_rc].deg == 0 && g->nodes[v].deg == 0) {
		/* all pairs are resolved */
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
	} else {
		g->edges[e].count = g->edges[e_rc].count = e_uni_count * ecov;
	}
}

void check_2_2_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t u, v, v_rc, u_rc, e1, e2, e_rc;
	e_rc = g->edges[e].rc_id;
	v = g->edges[e].target;
	v_rc = g->nodes[v].rc_id;
	u = g->edges[e].source;
	u_rc = g->nodes[u].rc_id;
	if (g->nodes[v].deg != 2 || g->nodes[u_rc].deg != 2 ||
		g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1)
		return;
	int ecov, cov1, cov2;
	uint64_t e_uni_count;
	ecov = (int)(__get_edge_cov(g->edges + e, g->ksize) / uni_cov + 0.499999);
	e_uni_count = g->edges[e].count / ecov;
	__VERBOSE("2-2 bridge: %ld_%ld ~ %d\n", e, g->edges[e].rc_id, ecov);
	int i, k;
	// for (i = 0; i < 2; ++i) {
	// 	e1 = g->nodes[u_rc].adj[i];
	// 	cov = (int)(__get_edge_cov(g->edges + e1, g->ksize) / uni_cov + 0.499999);
	// 	__VERBOSE("Leg left %i: %ld_%ld ~ %d\n",
	// 		i, e1, g->edges[e1].rc_id, cov);
	// }
	// for (i = 0; i < 2; ++i) {
	// 	e2 = g->nodes[v].adj[i];
	// 	cov = (int)(__get_edge_cov(g->edges + e2, g->ksize) / uni_cov + 0.499999);
	// 	__VERBOSE("Leg right %i: %ld_%ld ~ %d\n",
	// 		i, e2, g->edges[e2].rc_id, cov);
	// }

	int resolve;
	do {
		resolve = 0;
		for (i = 0; i < g->nodes[u_rc].deg; ++i) {
			e1 = g->nodes[u_rc].adj[i];
			e2 = bc_find_best_pair(g, e1, g->nodes[v].adj, g->nodes[v].deg);
			if (e2 == -1)
				continue;
			gint_t et1 = bc_find_best_pair(g, e2, g->nodes[u_rc].adj, g->nodes[u_rc].deg);
			__VERBOSE("e1 = %ld; e2 = %ld; et1 = %ld\n", e1, e2, et1);
			if (et1 != e1)
				continue;
			/* join edge e1.rc -> e -> e2 */
			cov1 = (int)(__get_edge_cov(g->edges + e1, g->ksize) / uni_cov + 0.4999999);
			cov2 = (int)(__get_edge_cov(g->edges + e2, g->ksize) / uni_cov + 0.4999999);
			if (cov1 != cov2 || cov1 != 1)
				continue;
			__VERBOSE("Joining edge %ld -> %ld -> %ld\n", e1, e, e2);
			asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
				e2, g->edges[e2].rc_id, e_uni_count * cov1);
			resolve = 1;
			break;
			// g->edges[e].count = g->edges[e_rc].count = g->edges[e].count * (ecov - cov1) / ecov;
			// ecov -= cov1;
		}
	} while (resolve);

	/* special case when resolved all but 1 pair left */
	if (g->nodes[u_rc].deg == 1 && g->nodes[v].deg == 1) {
		e1 = g->nodes[u_rc].adj[0];
		e2 = g->nodes[v].adj[0];
		cov1 = (int)(__get_edge_cov(g->edges + e1, g->ksize) / uni_cov + 0.4999999);
		cov2 = (int)(__get_edge_cov(g->edges + e2, g->ksize) / uni_cov + 0.4999999);
		double ratio = get_barcode_ratio(g, e1, e2);
		/* FIXME: may need more complicated resolve here */
		__VERBOSE("Leftover path: %ld -> %ld -> %ld, ratio: %.6lf\n",
				e1, e, e2, ratio);
		if (__positive_ratio(ratio) ||
			(ratio < 0 && ecov >= cov1 && cov1 == cov2)) {
			__VERBOSE("Joining edge %ld -> %ld -> %ld\n", e1, e, e2);
			asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
				e2, g->edges[e2].rc_id, e_uni_count * cov1);
		}
		/* destroy the bridge */
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
	} else if (g->nodes[u_rc].deg == 0 && g->nodes[v].deg == 0) {
		/* all pairs are resolved */
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
	}
}

void collapse_4_leg_complex(struct asm_graph_t *g, khash_t(gint) *set_e,
				khash_t(gint) *set_leg, double uni_cov)
{
	khiter_t k;
	gint_t legs[10];
	int n_leg, cov1, cov2;
	n_leg = 0;
	for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
		if (kh_exist(set_leg, k)) {
			legs[n_leg++] = kh_key(set_leg, k);
			kh_del(gint, set_e, kh_key(set_leg, k));
			kh_del(gint, set_e, g->edges[kh_key(set_leg, k)].rc_id);
		}
	}
	assert(n_leg == 4);
	/* only one bridge, use simpler version */
	if (kh_size(set_e) == 2) {
		__VERBOSE("2-1-2 bridge\n");
		gint_t e = -1;
		for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
			if (kh_exist(set_e, k)) {
				e = kh_key(set_e, k);
				break;
			}
		}
		check_2_2_bridge(g, e, uni_cov);
		return;
	}

	// gint_t e1, e2;
// 	e1 = legs[0];
// 	e2 = legs[1];
// 	cov1 = (int)(__get_edge_cov(g->edges + e1, g->ksize) / uni_cov + 0.499999);
// 	cov2 = (int)(__get_edge_cov(g->edges + e2, g->ksize) / uni_cov + 0.499999);
// 	if (cov1 != cov2) {
// 		__VERBOSE("WARNING: not collapse 1-1 jungle\n");
// 		return;
// 	}
// 	if (!check_path(g, set_e, g->edges[legs[0]].source,
// 			g->nodes[g->edges[legs[1]].source].rc_id)) {
// 		__VERBOSE("WARNING: not collapse 1-1 jungle, no path\n");
// 		return;
// 	}
	uint32_t gap_size = 0;
	uint64_t gap_count = 0;
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (!kh_exist(set_e, k))
			continue;
		gint_t e = kh_key(set_e, k);
		gint_t len = get_edge_len(g->edges + e);
		int cov = (int)(__get_edge_cov(g->edges + e, g->ksize) / uni_cov + 0.499999);
		gap_size += cov * (len - g->ksize);
		gap_count += g->edges[e].count;
		/* Remove edges , e.i isolate the nodes */
		// __VERBOSE("Removing edge %ld\n", e);
		// asm_remove_edge(g, e);
	}
	gap_size /= 2;
	gap_count /= 2;
	int resolve;
	do {
		resolve = 0;
		int i;
		gint_t e1, e2, et1;
		for (i = 0; i < n_leg; ++i) {
			e1 = legs[i];
			if (g->edges[e1].source == -1)
				__ERROR("Noob!");
			e2 = bc_find_best_pair(g, e1, legs, n_leg);
			if (g->edges[e2].source == -1)
				__ERROR("Noob!");
			et1 = bc_find_best_pair(g, e2, legs, n_leg);
			if (e1 != et1) {
				__VERBOSE("WARNING: not best pair %ld - %ld - %ld\n", e1, e2, et1);
				continue;
			}
			if (!check_path(g, set_e,
					g->nodes[g->edges[e1].source].rc_id,
					g->edges[e2].source)) {
				__VERBOSE("WARNING: not join edge %ld - %ld, no path\n", e1, e2);
				continue;
			}
			__VERBOSE("Joining edge with gap %ld -> %ld\n", e1, e2);
			asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1, e2,
				g->edges[e2].rc_id, gap_size / 2, gap_count / 2);
			/* remove leg */
			gint_t j;
			j = find_adj_idx(legs, n_leg, e1);
			if (j != -1)
				legs[j] = legs[--n_leg];
			j = find_adj_idx(legs, n_leg, e2);
			if (j != -1)
				legs[j] = legs[--n_leg];
			// asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
			// 	e2, g->edges[e2].rc_id, e_uni_count * cov1);
			resolve = 1;
			break;
		}
	} while (resolve);
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (!kh_exist(set_e, k))
			continue;
		gint_t e = kh_key(set_e, k);
		/* Remove edges , e.i isolate the nodes */
		__VERBOSE("Removing edge %ld\n", e);
		asm_remove_edge(g, e);
	}
}

void collapse_2_2_jungle(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	__VERBOSE("Genome coverage: %.9lf\n", uni_cov);
	khash_t(gint) *visited, *set_e, *set_v, *set_leg;
	visited = kh_init(gint);
	set_e = kh_init(gint);
	set_v = kh_init(gint);
	set_leg = kh_init(gint);
	gint_t e;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		uint32_t len = get_edge_len(g->edges + e);
		if (kh_get(gint, visited, g->edges[e].target) != kh_end(visited) ||
			len < MIN_SCAFFOLD_LEN || g->edges[e].source == -1)
			continue;
		find_region(g, e, MIN_SCAFFOLD_LEN, MAX_EDGE_COUNT, set_v, set_e);
		if (kh_size(set_e) < MAX_EDGE_COUNT) {
			kh_merge_set(visited, set_v);
			detect_leg(g, set_v, set_e, set_leg);
			// __VERBOSE("Edge: %ld; len: %u; Edge count in region: %u; Number of leg: %u\n",
			// 	e, len, kh_size(set_e), kh_size(set_leg));
			if (kh_size(set_leg) == 4) {
				__VERBOSE("Edge: %ld, len = %u; Edge count in region: %u; Number of leg: %u\n",
					e, len, kh_size(set_e), kh_size(set_leg));
				collapse_4_leg_complex(g, set_e, set_leg, uni_cov);
			}
				// collapse_1_1_complex(g, set_e, set_leg, uni_cov);
			// if (kh_size(set_e) < 20)
			// 	print_set_info(set_e);
		}
		kh_clear(gint, set_leg);
		kh_clear(gint, set_e);
		kh_clear(gint, set_v);
	}
}

void collapse_n_m_bridge(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t e;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		check_n_m_bridge(g, e, uni_cov);
	}
}

void collapse_n_m_jungle(struct asm_graph_t *g0, struct asm_graph_t *g1)
{
	collapse_2_2_jungle(g0);
	collapse_n_m_bridge(g0);
	asm_condense(g0, g1);
}

