#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "io_utils.h"
#include "khash.h"
#include "resolve.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

#define LEN_VAR				20
#define MAX_JOIN_LEN			2000
#define __get_edge_cov_int(g, e, uni_cov) (int)((g)->edges[e].count * 1.0 /    \
	((g)->edges[e].seq_len - ((g)->edges[e].n_holes + 1) * (g)->ksize) /   \
	(uni_cov) + 0.499999999)

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

static inline void asm_add_node_adj(struct asm_graph_t *g, gint_t u, gint_t e)
{
	g->nodes[u].adj = realloc(g->nodes[u].adj, (g->nodes[u].deg + 1) * sizeof(gint_t));
	g->nodes[u].adj[g->nodes[u].deg++] = e;
}

static inline void asm_remove_node_adj(struct asm_graph_t *g, gint_t u, gint_t e)
{
	gint_t j = find_adj_idx(g->nodes[u].adj, g->nodes[u].deg, e);
	if (j == -1)
		return;
	g->nodes[u].adj[j] = g->nodes[u].adj[--g->nodes[u].deg];
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
					if (g0->nodes[v].deg != 1) {
						fprintf(stderr, "Node %ld, deg = %ld\n",
							v, g0->nodes[v].deg);
						assert(0 && "Middle node degree is not equal to 1");
					}
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
 * return 1: self - loop
 * return 2: self - loop - reverse
 * return 3: double - loop
 * return -1: false loop, just go through
 */
int check_simple_loop(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t e_rc, u, v, u_rc, v_rc, e1, e2, e_rc1, e_rc2, e_return, e_return_rc;
	int cov, rep, rep_e, rep_e_return;
	double fcov, fcov1, fcov2, fcov_mean;
	e_rc = g->edges[e].rc_id;
	u = g->edges[e].source;
	v = g->edges[e].target;
	u_rc = g->nodes[u].rc_id;
	v_rc = g->nodes[v].rc_id;
	cov = __get_edge_cov_int(g, e, uni_cov);
	fcov = __get_edge_cov(g->edges + e, g->ksize);
	if (u == v) { /* self loop */
		if (cov == 0 && get_edge_len(g->edges + e) < MIN_NOTICE_LEN) {
			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
			return -1;
		}
		if (g->nodes[u_rc].deg > 2 || g->nodes[u].deg > 2)
			return 0;
		// asm_duplicate_edge_seq(g, e, cov);
		// asm_duplicate_edge_seq(g, e_rc, cov);
		/* split the node */
		v = asm_create_node(g);
		v_rc = g->nodes[v].rc_id;
		/* set source-target of edge e */
		g->edges[e].target = v;
		asm_remove_node_adj(g, u_rc, e_rc);
		g->edges[e_rc].source = v_rc;
		asm_add_node_adj(g, v_rc, e_rc);
		/* move edge from node u to node v */
		g->nodes[v].adj = malloc(g->nodes[u].deg * sizeof(gint_t));
		memcpy(g->nodes[v].adj, g->nodes[u].adj, g->nodes[u].deg * sizeof(gint_t));
		g->nodes[v].deg = g->nodes[u].deg;
		asm_remove_node_adj(g, v, e);
		/* node u has only edge e left */
		g->nodes[u].adj = realloc(g->nodes[u].adj, sizeof(gint_t));
		g->nodes[u].adj[0] = e;
		g->nodes[u].deg = 1;
		gint_t j;
		for (j = 0; j < g->nodes[v].deg; ++j) {
			gint_t e_t = g->nodes[v].adj[j];
			// fprintf(stderr, "e_t[%ld_%ld] %ld->%ld\n",
			// 	e_t, g->edges[e_t].rc_id, g->edges[e_t].source,
			// 	g->edges[e_t].target);
			g->edges[e_t].source = v;
			g->edges[g->edges[e_t].rc_id].target = v_rc;
		}
		return 1;
	} else if (u == v_rc) { /* self loop reverse */
		if (cov == 0) {
			__VERBOSE("Remove edges %ld[len=%u]\n", e, get_edge_len(g->edges + e));
			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
			return -1;
		}
		// if (g->nodes[u].deg > 2)
		// 	return 0;
		// if (g->nodes[v].deg == 1) {
		// 	e_rc1 = g->nodes[v].adj[0];
		// 	e1 = g->edges[e_rc1].rc_id;
		// 	asm_join_edge_loop_reverse(g, e1, e, e_rc, e_rc1);
		// 	asm_remove_edge(g, e);
		// 	asm_remove_edge(g, e_rc);
		// 	return 2;
		// } else if (g->nodes[v].deg == 2) {
		// 	e_rc1 = g->nodes[v].adj[0];
		// 	e1 = g->edges[e_rc1].rc_id;
		// 	e2 = g->nodes[v].adj[1];
		// 	e_rc2 = g->edges[e2].rc_id;
		// 	asm_join_edge3(g, e1, e_rc1, e, e_rc, e2, e_rc2,
		// 				g->edges[e].count);
		// 	asm_remove_edge(g, e);
		// 	asm_remove_edge(g, e_rc);
		// 	return 2;
		// }
		return 0;
	} else {
		if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
			g->nodes[u_rc].deg > 2 || g->nodes[v].deg > 2)
			return 0;
		e_return = -1;
		gint_t j;
		for (j = 0; j < g->nodes[v].deg; ++j) {
			if (g->edges[g->nodes[v].adj[j]].target == u) {
				e_return = g->nodes[v].adj[j];
				break;
			}
		}
		if (e_return == -1)
			return 0;
		e_return_rc = g->edges[e_return].rc_id;
		fcov1 = fcov2 = -1;
		for (j = 0; j < g->nodes[v].deg; ++j) {
			if (g->nodes[v].adj[j] != e_return)
				fcov1 = __get_edge_cov(g->edges + g->nodes[v].adj[j], g->ksize);
		}
		for (j = 0; j < g->nodes[u_rc].deg; ++j) {
			if (g->nodes[u_rc].adj[j] != e_return_rc)
				fcov2 = __get_edge_cov(g->edges + g->nodes[u_rc].adj[j], g->ksize);
		}
		if (fcov1 > 0 && fcov2 > 0) {
			if ((int)(fcov1 / fcov2 + 0.499999999) > 2)
				return 0;
			fcov_mean = (fcov1 + fcov2) / 2;
		} else if (fcov1 > 0) {
			fcov_mean = fcov1;
		} else if (fcov2 > 0) {
			fcov_mean = fcov2;
		} else {
			fcov_mean = uni_cov;
		}
		rep_e = __get_edge_cov_int(g, e, fcov_mean) - 1;
		rep_e_return = __get_edge_cov_int(g, e_return, fcov_mean);
		if (rep_e_return == 0) {
			asm_remove_edge(g, e_return);
			asm_remove_edge(g, e_return_rc);
			return -1;
		}
		if (g->edges[e].seq_len > g->edges[e_return].seq_len)
			rep = rep_e - 1;
		else
			rep = rep_e_return;
		if (rep > 0) {
			asm_duplicate_edge_seq2(g, e, e_return, rep);
			asm_duplicate_edge_seq2(g, e_rc, e_return_rc, rep);
		}
		asm_remove_edge(g, e_return);
		asm_remove_edge(g, e_return_rc);
		return 3;
	}
	return 0;
}

int unroll_simple_loop(struct asm_graph_t *g, double uni_cov)
{
	gint_t e, cnt_self, cnt_self_reverse, cnt_double, cnt_false;
	int ret;
	cnt_self = cnt_self_reverse = cnt_double = cnt_false = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		ret = check_simple_loop(g, e, uni_cov);
		cnt_self += (ret == 1);
		cnt_self_reverse += (ret == 2);
		cnt_double += (ret == 3);
		cnt_false += (ret == -1);
	}
	__VERBOSE("Number of unroll: self loop (%ld), self loop reverse (%ld), double loop (%ld), false loop (%ld)\n",
		cnt_self, cnt_self_reverse, cnt_double, cnt_false);
	return cnt_self + cnt_self_reverse + cnt_double + cnt_false;
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
			len + get_edge_len(g->edges + e) - ksize, max_len, ksize);
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
					get_edge_len(g0->edges + e_id) - g0->ksize, MAX_TIPS_LEN, g0->ksize);
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
			if (cov / max_cov < TIPS_RATIO_THRESHOLD && max_cov < 600.0) {
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
	gint_t u, v, u_rc, v_rc, e_rc, ec, ec_rc, j;
	u = g->edges[e].source;
	v = g->edges[e].target;
	e_rc = g->edges[e].rc_id;
	u_rc = g->nodes[u].rc_id;
	v_rc = g->nodes[v].rc_id;
	if (g->nodes[u_rc].deg < 1 || g->nodes[u].deg != 1 || g->nodes[v].deg > 1)
		return 0;
	double e_cov, cov, sum_fcov = 0;
	struct cov_range_t e_rcov, rcov;
	int sum_min_cov = 0;
	e_cov = __get_edge_cov(g->edges + e, g->ksize) / uni_cov;
	e_rcov = convert_cov_range(e_cov);
	// fprintf(stderr, "consider edge %ld: cov = %.3lf\n", e, e_cov);
	uint32_t e_len, max_len;
	e_len = get_edge_len(g->edges + e);
	if (e_len > 2000)
		return -1;
	max_len = 0;
	for (j = 0; j < g->nodes[u_rc].deg; ++j) {
		gint_t n, n_rc;
		n_rc = g->nodes[u_rc].adj[j];
		n = g->edges[n_rc].rc_id;
		cov = __get_edge_cov(g->edges + n, g->ksize) / uni_cov;
		sum_fcov += cov;
		rcov = convert_cov_range(cov);
		sum_min_cov += rcov.lo;
		// fprintf(stderr, "\tsattelite: edge %ld ~ %.3lf\n", n, cov);
		uint32_t len = get_edge_len(g->edges + n);
		max_len = __max(max_len, len);
	}
	if (e_len > MIN_NOTICE_LEN || max_len > MIN_NOTICE_LEN) {
		// if (e_rcov.hi < sum_min_cov || !__diff_accept(sum_fcov, e_cov))
		if (e_rcov.hi < sum_min_cov || e_cov + 0.5 < sum_fcov)
			return -1;
	}
	for (j = 0; j < g->nodes[u_rc].deg; ++j) {
		gint_t n, n_rc;
		n_rc = g->nodes[u_rc].adj[j];
		n = g->edges[n_rc].rc_id;
		if (n == e || n_rc == e)
			return -1;
	}
	// g->nodes[v_rc].adj = realloc(g->nodes[v_rc].adj,
	// 	(g->nodes[v_rc].deg + g->nodes[u_rc].deg) * sizeof(gint_t));
	while (g->nodes[u_rc].deg) {
		gint_t n, n_rc;
		n_rc = g->nodes[u_rc].adj[0];
		n = g->edges[n_rc].rc_id;
		ec = asm_clone_edge3(g, e);
		ec_rc = g->edges[ec].rc_id;
		cov = __get_edge_cov(g->edges + n, g->ksize) / uni_cov;
		g->edges[ec].count = g->edges[ec_rc].count =
				(uint64_t)(cov / sum_fcov * g->edges[e].count);

		asm_join_edge(g, n, n_rc, ec, ec_rc);
		// asm_append_edge_seq(g->edges + n, g->edges + e, g->ksize);
		// g->edges[n].count += split_count;
		// g->edges[n].target = v;

		// asm_clean_edge_seq(g->edges + n_rc);
		// asm_clone_reverse(g->edges + n_rc, g->edges + n);
		// g->edges[n_rc].source = v_rc;
		// g->edges[n_rc].target = g->nodes[g->edges[n].source].rc_id;
		// g->nodes[v_rc].adj[g->nodes[v_rc].deg++] = n_rc;
		// g->edges[n].rc_id = n_rc;
		// g->edges[n_rc].rc_id = n;
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
	// __VERBOSE("Number of edges cannot expand: %d\n", cnt_fp);
	return cnt;
}

int test_bubble2(struct asm_graph_t *g, gint_t u)
{
	if (g->nodes[u].deg < 2)
		return 0;
	gint_t j, k, e, ke, v, deg, best_e;
	int resolve, count_resolve = 0;
	do {
		resolve = 0;
		for (j = 0; j < g->nodes[u].deg; ++j) {
			e = g->nodes[u].adj[j];
			if (v == g->nodes[u].rc_id)
				continue;
			v = g->edges[e].target;
			uint32_t max_len, base_len;
			base_len = get_edge_len(g->edges + e);
			gint_t idx = -1;
			uint64_t cur_count, sum_count;
			sum_count = 0;
			int cnt = 0;
			for (k = 0; k < g->nodes[u].deg; ++k) {
				e = g->nodes[u].adj[k];
				if (v != g->edges[e].target)
					continue;
				uint32_t len = get_edge_len(g->edges + e);
				if ((len >= MAX_JOIN_LEN || base_len >= MAX_JOIN_LEN) &&
					(base_len - LEN_VAR > len || base_len + LEN_VAR < len))
					continue;
				if (g->edges[e].count > cur_count) {
					cur_count = g->edges[e].count;
					best_e = e;
				}
				max_len = __max(max_len, len);
				sum_count += g->edges[e].count;
				++cnt;
			}
			if (cnt < 2)
				continue;
			g->edges[best_e].count = g->edges[g->edges[best_e].rc_id].count = sum_count;
			for (k = 0; k < g->nodes[u].deg; ++k) {
				e = g->nodes[u].adj[k];
				if (v != g->edges[e].target)
					continue;
				if (e != best_e) {
					asm_remove_edge(g, g->edges[e].rc_id);
					g->nodes[u].adj[k] = -1;
					g->edges[e].source = g->edges[e].target = -1;
				}
			}
			deg = 0;
			for (k = 0; k < g->nodes[u].deg; ++k) {
				if (g->nodes[u].adj[k] != -1)
					g->nodes[u].adj[deg++] = g->nodes[u].adj[k];
			}
			g->nodes[u].deg = deg;
			resolve = 1;
			break;
		}
		count_resolve += resolve;
	} while (resolve);
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

gint_t remove_low_cov_edge(struct asm_graph_t *g0, struct asm_graph_t *g1)
{
	double uni_cov, cov;
	struct cov_range_t rcov;
	gint_t e, e_rc;
	uni_cov = get_genome_coverage(g0);
	gint_t cnt = 0;
	for (e = 0; e < g0->n_e; ++e) {
		if (g0->edges[e].source == -1)
			continue;
		e_rc = g0->edges[e].rc_id;
		cov = __get_edge_cov(g0->edges + e, g0->ksize) / uni_cov;
		rcov = convert_cov_range(cov);
		if (rcov.hi == 0) {
			asm_remove_edge(g0, e);
			asm_remove_edge(g0, e_rc);
			++cnt;
		}
	}
	asm_condense(g0, g1);
	__VERBOSE("Number of low coverage edge removed: %ld\n", cnt);
	return cnt;
}

void resolve_chain(struct asm_graph_t *g0, struct asm_graph_t *g1)
{
	int step = 0;
	while (1) {
		while (1) {
			__VERBOSE("Iteration [%d]\n", step++);
			gint_t cnt_expand, cnt_collapse, cnt_loop;
			double uni_cov = get_genome_coverage(g0);
			__VERBOSE("Genome coverage: %.9lf\n", uni_cov);
			cnt_loop = unroll_simple_loop(g0, uni_cov);
			cnt_collapse = resolve_bubble2(g0);
			cnt_expand = graph_expanding(g0, uni_cov);
			asm_condense(g0, g1);
			test_asm_graph(g1);
			*g0 = *g1;
			if (cnt_expand == 0 && cnt_collapse == 0 && cnt_loop == 0)
				break;
		}
		gint_t cnt_trim = remove_low_cov_edge(g0, g1);
		if (cnt_trim == 0)
			break;
		*g0 = *g1;
	}
}

