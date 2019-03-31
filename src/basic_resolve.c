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

/* return 0: not at all
 * return 1: loop
 * return 2: false loop, just go through
 */
static inline int check_simple_loop(struct asm_graph_t *g, gint_t e)
{
	if (g->edges[e].source == -1) /* edge is removed */
		return 0;
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
	float cov_left, cov_right, cov_ahead, cov_return;
	cov_left = __get_edge_cov(g->edges + e_left, g->ksize);
	cov_right = __get_edge_cov(g->edges + e_right, g->ksize);
	cov_ahead = __get_edge_cov(g->edges + e, g->ksize);
	cov_return = __get_edge_cov(g->edges + e_return, g->ksize);

	if (__int_cov_ratio(cov_left, cov_right) != 1)
		return 0;
	if (__int_cov_ratio(cov_ahead, cov_return) < 2) {
		if (__int_cov_ratio(cov_left, cov_ahead) == 1) {
			asm_remove_edge(g, e_return);
			asm_remove_edge(g, e_return_rc);
			/* remove the return edge */
			return 2;
		} else {
			/* still resolve */
			asm_append_edge_seq(g->edges + e_return,
						g->edges + e, g->ksize);
			asm_append_edge_seq(g->edges + e,
						g->edges + e_return, g->ksize);
			g->edges[e].count += g->edges[e_return].count;
			asm_append_edge_seq(g->edges + e_return_rc,
						g->edges + e_rc, g->ksize);
			asm_append_edge_seq(g->edges + e_rc,
					g->edges + e_return_rc, g->ksize);
			asm_remove_edge(g, e_return);
			asm_remove_edge(g, e_return_rc);
			return 1;
		}
	} else {
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
		asm_remove_edge(g, e_return);
		asm_remove_edge(g, e_return_rc);
		return 1;
	}
}

void unroll_simple_loop(struct asm_graph_t *g)
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
		ret = check_simple_loop(g, e);
		if (ret == 1)
			++cnt_loop;
		else if (ret == 2)
			++cnt_removed;
	}
	__VERBOSE("Number of unroll loops: %ld\n", cnt_loop);
	__VERBOSE("Number of false loops: %ld\n", cnt_removed);
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
	unroll_simple_loop(g0);
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
	double cov, sum_fcov;
	int e_cov, sum_cov;
	cov = __get_edge_cov(g->edges + e, g->ksize);
	e_cov = (int)(cov / uni_cov + 0.499999999);
	sum_cov = 0;
	sum_fcov = 0.0;
	gint_t j;
	uint32_t e_len, max_len;
	e_len = get_edge_len(g->edges + e);
	max_len = 0;
	for (j = 0; j < g->nodes[u_rc].deg; ++j) {
		gint_t n, n_rc;
		n_rc = g->nodes[u_rc].adj[j];
		n = g->edges[n_rc].rc_id;
		cov = __get_edge_cov(g->edges + n, g->ksize);
		sum_fcov += cov;
		sum_cov += (int)(cov / uni_cov + 0.499999999);
		uint32_t len = get_edge_len(g->edges + n);
		max_len = __max(max_len, len);
	}
	fprintf(stderr, "consider edge %ld: cov = %d, sum sattelite cov = %d\n",
		e, e_cov, sum_cov);
	// fprintf(stderr, "consider edge %ld: cov = %.3f, sum sattelite cov = %.3f\n",
		// e, e_cov, sum_cov);
	if (e_len > 500 || max_len > 500) {
		if (e_cov < sum_cov)
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
		// fprintf(stderr, "concat edge: (%ld %ld) -> (%ld %ld); (%ld %ld) -> (%ld %ld)\n",
		// 	g->edges[n].source, g->edges[n].target,
		// 	g->edges[e].source, g->edges[e].target,
		// 	g->edges[e_rc].source, g->edges[e_rc].target,
		// 	g->edges[n_rc].source, g->edges[n_rc].target);
		cov = __get_edge_cov(g->edges + n, g->ksize);
		uint64_t split_count = (uint64_t)(cov / sum_fcov * g->edges[e].count);

		asm_append_edge_seq(g->edges + n, g->edges + e, g->ksize);
		g->edges[n].count += split_count;
		g->edges[n].target = v;

		struct asm_edge_t tmp;
		asm_clone_edge(&tmp, g->edges + e_rc);
		tmp.count = g->edges[n_rc].count + split_count;
		asm_append_edge_seq(&tmp, g->edges + n_rc, g->ksize);
		tmp.source = v_rc;
		tmp.target = g->edges[n_rc].target;
		asm_clean_edge_seq(g->edges + n_rc);

		g->edges[n_rc] = tmp;
		g->nodes[v_rc].adj[g->nodes[v_rc].deg++] = n_rc;
		g->edges[n].rc_id = n_rc;
		g->edges[n_rc].rc_id = n;
		assert(g->edges[n].source == g->nodes[g->edges[n_rc].target].rc_id);
		assert(g->edges[n].target == g->nodes[g->edges[n_rc].source].rc_id);
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
		gint_t idx = -1;
		int cnt = 0;
		for (k = 0; k < g->nodes[u].deg; ++k) {
			e = g->nodes[u].adj[k];
			if (e == -1 || v != g->edges[e].target)
				continue;
			cov = __get_edge_cov(g->edges + e, g->ksize);
			if (cov > cur_cov) {
				cur_cov = cov;
				idx = k;
			}
			uint32_t len = get_edge_len(g->edges + e);
			max_len = __max(max_len, len);
			++cnt;
		}
		if (cnt < 2 || v == g->nodes[u].rc_id || max_len > 500) {
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
		int cnt_expand, cnt_collapse;
		double uni_cov = get_genome_coverage(g0);
		__VERBOSE("Genome coverage: %.9lf\n", uni_cov);
		cnt_expand = graph_expanding(g0, uni_cov);
		cnt_collapse = resolve_bubble2(g0);
		asm_condense(g0, g1);
		test_asm_graph(g1);
		if (cnt_expand == 0 && cnt_collapse == 0)
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
		// u = g->edges[e].source;
		// u_rc = g->nodes[u].rc_id;
		// for (j = 0; j < g->nodes[u_rc].deg; ++j) {
		// 	gint_t ne = g->nodes[u_rc].adj[j];
		// 	uint32_t len = get_edge_len(g->edges + ne);
		// 	if (len < min_contig_len) {
		// 		if (cb_add(set_e, ne) == 1) {
		// 			if (edge_count == max_edge_count)
		// 				goto clean_up;
		// 			q[++r] = ne;
		// 		}
		// 		cb_add(set_v, g->edges[ne].target);
		// 	} else {
		// 		cb_add(set_e, ne);
		// 		cb_add(set_v, g->edges[ne].target);
		// 	}
		// }
	}
clean_up:
	free(q);
}

// struct dfs_info_t {
// 	gint_t num;
// 	gint_t low;
// 	gint_t parent;
// };

// static void rm_edge(struct asm_graph_t *g, gint_t e, khash_t(gint) *set_e)
// {
// 	gint_t e_rc = g->nodes[e].rc_id;
// 	if (e > e_rc)
// 		kh_del(gint, set_e, e_rc);
// 	else
// 		kh_del(gint, set_e, e);
// }

// static struct dfs_info_t *add_node(struct asm_graph_t *g, gint_t u,
// 							khash_t(dfs) *nodes)
// {
// 	khiter_t k;
// 	int ret;
// 	gint_t u_rc = g->nodes[u].rc_id;
// 	if (u > u_rc)
// 		k = kh_put(dfs, nodes, u_rc, &ret);
// 	else
// 		k = kh_put(dfs, nodes, u, &ret);
// 	if (ret == 1)
// 		kh_value(nodes, k).parent = -1;
// 	return &(kh_value(nodes, k));
// }

// void dfs_bridge(struct asm_graph_t *g, gint_t u, gint_t *cnt,
// 		khash_t(dfs) *nodes, khash_t(gint) *set_e, khash_t(gint) *leg)
// {
// 	struct dfs_info_t *it_u, *it_v;
// 	it_u = add_node(g, u, nodes);
// 	it_u->num = cnt;
// 	it_u->low = cnt;
// 	++*cnt;
// 	u_rc = u;
// 	for (j = 0; j < g->nodes[u].deg; ++j) {
// 		e = g->nodes[u].adj[j];
// 		/* edges not present on graph */
// 		if (kh_get(gint, set_e, e) == kh_end(set_e))
// 			continue;
// 		v = g->nodes[e].target;
// 		it_v = add_node(g, v, nodes);
// 		if (it_v->parent == -1) {
// 			it_v->parent = u;
// 			dfs_bridge(g, v, cnt, nodes, set_e, leg);
// 			it_u->low = __min(it_u->low, it_v->low);
// 		} else {
// 			it_u->low  __min(it_u->low, it_v->num);
// 		}
// 		if (it_v->low > it_u->num)
// 			kh_put(gint, leg, e, &ret);
// 		rm_edge(g, e, set_e);
// 	}

// 	for (j = 0; j < g->nodes[u_rc].deg; ++j) {
// 		e = g->nodes[u_rc].adj[j];
// 		if (kh_get(gint, set_e, e) == kh_end(set_e))
// 			continue;
// 		v = g->nodes[e].target;
// 		it_v = add_node(g, v, nodes);
// 		if (it_v->parent == -1) {
// 			it_v->parent = u;
// 			dfs_bridge(g, v, cnt, nodes, set_e, leg);
// 			it_u->low = __min(it_u->low, it_v->low);
// 		} else {
// 			it_u->low = __min(it_u->low, it_v->num);
// 		}
// 		if (it_v->low > it_u->num)
// 			kh_put(gint, leg, e, &ret);
// 		rm_edge(g, e, set_e);
// 	}
// }

// void calibrate_leg(struct asm_graph_t *g, gint_t se, khash_t(gint) *set_v,
// 	khash_t(gint) *set_e, khash_t(gint) *leg, , int min_contig_len)
// {
// 	gint_t cnt;
// }

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

void detect_complex(struct asm_graph_t *g, uint32_t min_contig_size, uint32_t max_edge_count)
{
	khash_t(gint) *visited, *set_e, *set_v;
	visited = kh_init(gint);
	set_e = kh_init(gint);
	set_v = kh_init(gint);
	gint_t e;
	for (e = 0; e < g->n_e; ++e) {
		uint32_t len = get_edge_len(g->edges + e);
		if (kh_get(gint, visited, g->edges[e].target) != kh_end(visited) || len < min_contig_size)
			continue;
		find_region(g, e, min_contig_size, max_edge_count, set_v, set_e);
		if (kh_size(set_e) < max_edge_count) {
			kh_merge_set(visited, set_v);
			__VERBOSE("Edge: %ld; len: %u; Edge count in region: %u\n",
				e, len, kh_size(set_e));
			if (kh_size(set_e) < 20)
				print_set_info(set_e);
		}
		kh_clear(gint, set_e);
		kh_clear(gint, set_v);
	}
}
