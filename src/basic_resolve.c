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

void asm_lazy_condense(struct asm_graph_t *g)
{
	/* remove unused links */
	gint_t u, u_rc, deg_fw, deg_rv, e1, e2;
	for (u = 0; u < g->n_v; ++u) {
		gint_t c, deg;
		deg = g->nodes[u].deg;
		g->nodes[u].deg = 0;
		for (c = 0; c < deg; ++c) {
			gint_t e_id = g->nodes[u].adj[c];
			if (e_id != -1)
				g->nodes[u].adj[g->nodes[u].deg++] = e_id;
		}
		g->nodes[u].adj = realloc(g->nodes[u].adj,
					g->nodes[u].deg * sizeof(gint_t));
	}

	/* join non-branching path */
	for (u = 0; u < g->n_v; ++u) {
		u_rc = g->nodes[u].rc_id;
		deg_fw = g->nodes[u].deg;
		deg_rv = g->nodes[u_rc].deg;
		if (deg_fw == 1 && deg_rv == 1) {
			e1 = g->nodes[u].adj[0];
			e2 = g->nodes[u_rc].adj[0];
			if (e1 == e2 || e1 == g->edges[e2].rc_id)
				continue;
			asm_join_edge(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
		}
	}
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
			asm_clone_seq(edges + p, g0->edges + e);
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
					asm_append_seq(edges + p,
						g0->edges + e, g0->ksize);
					edges[p].count += g0->edges[e].count;
				} else {
					break;
				}
			} while (1);
			edges[p].source = x;
			edges[p].target = node_id[v];
			asm_clone_seq_reverse(edges + q, edges + p);
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
	g->aux_flag = 0;
	g->n_v = n_v;
	g->n_e = n_e;
	g->nodes = nodes;
	g->edges = edges;
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

void find_topo(struct asm_graph_t *g, gint_t *d, gint_t *deg, uint32_t max_len)
{
	gint_t u, u_rc, v, v_rc, e, l, r, j, ksize;
	gint_t *q;
	q = malloc(g->n_v * sizeof(gint_t));
	l = 0; r = -1;
	ksize = g->ksize;
	for (u = 0; u < g->n_v; ++u) {
		deg[u] = g->nodes[u].deg;
		if (g->nodes[u].deg == 0) {
			q[++r] = u;
			d[u] = 0;
		} else {
			d[u] = 0;
		}
	}
	while (l <= r) {
		u = q[l++];
		u_rc = g->nodes[u].rc_id;
		for (j = 0; j < g->nodes[u_rc].deg; ++j) {
			e = g->nodes[u_rc].adj[j];
			v_rc = g->edges[e].target;
			v = g->nodes[v_rc].rc_id;
			--deg[v];
			if (d[v] == -1 || d[u] + g->edges[e].seq_len - ksize > d[v])
				d[v] = g->edges[e].seq_len - ksize;
			if (d[v] > max_len)
				d[v] = max_len;
			if (deg[v] == 0) {
				q[++r] = v;
			}
		}
	}
	free(q);
}

gint_t remove_tips_topo(struct asm_graph_t *g)
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
	gint_t *d, *degs;
	d = malloc(g->n_v * sizeof(gint_t));
	degs = malloc(g->n_v * sizeof(gint_t));
	find_topo(g, d, degs, TIPS_LEN_THRES);

	gint_t u, u_rc, v, e, e_rc, j, cnt_removed;
	double cov, cov_fw, cov_rv, max_cov;
	uint32_t len_fw, len_rv, extend_left, extend_right;
	cnt_removed = 0;
	for (u = 0; u < g->n_v; ++u) {
		u_rc = g->nodes[u].rc_id;
		cov_fw = cov_rv = 0;
		len_fw = len_rv = 0;
		extend_left = extend_right = 0;
		for (j = 0; j < g->nodes[u].deg; ++j) {
			e = g->nodes[u].adj[j];
			cov = __get_edge_cov(g->edges + e, g->ksize);
			cov_fw = __max(cov_fw, cov);
			len_fw = __max(len_fw, g->edges[e].seq_len);
			v = g->edges[e].target;
			extend_left |= (degs[v] != 0 || d[v] == -1 || d[v] + g->edges[e].seq_len - g->ksize >= MIN_TIPS_LEG);
		}
		for (j = 0; j < g->nodes[u_rc].deg; ++j) {
			e = g->nodes[u_rc].adj[j];
			cov = __get_edge_cov(g->edges + e, g->ksize);
			cov_rv = __max(cov_rv, cov);
			len_rv = __max(len_rv, g->edges[e].seq_len);
			v = g->edges[e].target;
			extend_right |= (degs[v] != 0 || d[v] == -1 || d[v] + g->edges[e].seq_len - g->ksize >= MIN_TIPS_LEG);
		}
		max_cov = __max(cov_fw, cov_rv);
		for (j = 0; j < g->nodes[u].deg; ++j) {
			e = g->nodes[u].adj[j];
			v = g->edges[e].target;
			cov = __get_edge_cov(g->edges + e, g->ksize);
			if (degs[v] == 0 && d[v] != -1 && cov < max_cov &&
				((d[v] + g->edges[e].seq_len - g->ksize < TIPS_LEN_THRES && ((extend_left && extend_right && cov < 30) || (cov < cov_fw))) ||
				 (cov < TIPS_COV_THRES && cov < max_cov * TIPS_RATIO_THRES) ||
				 (len_fw >= MIN_TIPS_LEG && len_rv >= MIN_TIPS_LEG && cov < max_cov * TIPS_RATIO_THRES))) {
				e_rc = g->edges[e].rc_id;
				asm_remove_edge(g, e);
				asm_remove_edge(g, e_rc);
				++cnt_removed;
				--j;
			}
		}
	}
	free(d);
	free(degs);
	__VERBOSE("Number of tips remove using graph topology: %ld\n", cnt_removed);
	return cnt_removed;
}

gint_t remove_tips(struct asm_graph_t *g)
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
	gint_t u, u_rc, j, e, e_rc, v, cnt_removed;
	double cov_fw, cov_rv, cov, max_cov;
	uint32_t len_fw, len_rv, extend_left, extend_right;
	cnt_removed = 0;
	for (u = 0; u < g->n_v; ++u) {
		u_rc = g->nodes[u].rc_id;
		cov_fw = cov_rv = 0;
		len_fw = len_rv = 0;
		extend_left = extend_right = 0;
		for (j = 0; j < g->nodes[u].deg; ++j) {
			e = g->nodes[u].adj[j];
			cov = __get_edge_cov(g->edges + e, g->ksize);
			cov_fw = __max(cov_fw, cov);
			len_fw = __max(len_fw, g->edges[e].seq_len);
			v = g->edges[e].target;
			extend_left |= (g->nodes[v].deg != 0 || g->edges[e].seq_len >= MIN_TIPS_LEG);
		}
		for (j = 0; j < g->nodes[u_rc].deg; ++j) {
			e = g->nodes[u_rc].adj[j];
			cov = __get_edge_cov(g->edges + e, g->ksize);
			cov_rv = __max(cov_rv, cov);
			len_rv = __max(len_rv, g->edges[e].seq_len);
			v = g->edges[e].target;
			extend_right |= (g->nodes[v].deg != 0 || g->edges[e].seq_len >= MIN_TIPS_LEG);
		}
		max_cov = __max(cov_fw, cov_rv);
		for (j = 0; j < g->nodes[u].deg; ++j) {
			e = g->nodes[u].adj[j];
			v = g->edges[e].target;
			cov = __get_edge_cov(g->edges + e, g->ksize);
			if (g->nodes[v].deg == 0 && cov < max_cov &&
				((g->edges[e].seq_len < TIPS_LEN_THRES && extend_left && extend_right && cov < 30) ||
				 (cov < TIPS_COV_THRES && cov < max_cov * TIPS_RATIO_THRES) ||
				 (len_fw >= MIN_TIPS_LEG && len_rv >= MIN_TIPS_LEG && cov < max_cov * TIPS_RATIO_THRES))) {
				e_rc = g->edges[e].rc_id;
				asm_remove_edge(g, e);
				asm_remove_edge(g, e_rc);
				++cnt_removed;
				--j;
			}
		}
	}
	__VERBOSE("Number of trivial tips removed: %ld\n", cnt_removed);
	return cnt_removed;
}

static inline double get_max_out_cov(struct asm_graph_t *g, gint_t u)
{
	double cur_cov, cov;
	gint_t k, ep;
	cur_cov = 0.0;
	for (k = 0; k < g->nodes[u].deg; ++k) {
		ep = g->nodes[u].adj[k];
		if (g->edges[ep].source == -1)
			continue;
		cov = __get_edge_cov(g->edges + ep, g->ksize);
		cur_cov = __max(cur_cov, cov);
	}
	return cur_cov;
}

gint_t remove_chimeric(struct asm_graph_t *g)
{
	gint_t e, e_rc, u, u_rc, v, v_rc, cnt_removed;
	double cov, cov_fw, cov_rv;
	cnt_removed = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		e_rc = g->edges[e].rc_id;
		u = g->edges[e].source;
		u_rc = g->nodes[u].rc_id;
		v = g->edges[e].target;
		v_rc = g->nodes[v].rc_id;
		cov = __get_edge_cov(g->edges + e, g->ksize);
		cov_fw = get_max_out_cov(g, u);
		cov_fw = __min(cov_fw, get_max_out_cov(g, u_rc));
		cov_rv = get_max_out_cov(g, v);
		cov_rv = __min(cov_rv, get_max_out_cov(g, v_rc));
		if ((cov < CHIMERIC_RATIO_THRES * cov_fw ||
			cov < CHIMERIC_RATIO_THRES * cov_rv) &&
			cov < CHIMERIC_COV_THRES) {
			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
			++cnt_removed;
		}
	}
	__VERBOSE("Number of chimeric edge removed: %ld\n", cnt_removed);
	return cnt_removed;
}

int check_simple_loop(struct asm_graph_t *g, gint_t e)
{
	gint_t e_rc, u, v, u_rc, v_rc, e1, e2, e_rc1, e_rc2, e_return, e_return_rc;
	int rep, rep_e, rep_e_return;
	e_rc = g->edges[e].rc_id;
	u = g->edges[e].source;
	v = g->edges[e].target;
	u_rc = g->nodes[u].rc_id;
	v_rc = g->nodes[v].rc_id;
	if (u == v) { /* self loop */
		if (g->edges[e].seq_len < MIN_NOTICE_LEN) {
			/* split the node */
		}
	}
}

/* return 0: not at all
 * return 1: self - loop
 * return 2: self - loop - reverse
 * return 3: double - loop
 * return -1: false loop, just go through
 */
// int check_simple_loop(struct asm_graph_t *g, gint_t e, double uni_cov)
// {
// 	gint_t e_rc, u, v, u_rc, v_rc, e1, e2, e_rc1, e_rc2, e_return, e_return_rc;
// 	int cov, rep, rep_e, rep_e_return;
// 	double fcov, fcov1, fcov2, fcov_mean;
// 	e_rc = g->edges[e].rc_id;
// 	u = g->edges[e].source;
// 	v = g->edges[e].target;
// 	u_rc = g->nodes[u].rc_id;
// 	v_rc = g->nodes[v].rc_id;
// 	cov = __get_edge_cov_int(g, e, uni_cov);
// 	fcov = __get_edge_cov(g->edges + e, g->ksize);
// 	if (u == v) { /* self loop */
// 		if (cov == 0 && get_edge_len(g->edges + e) < MIN_NOTICE_LEN) {
// 			asm_remove_edge(g, e);
// 			asm_remove_edge(g, e_rc);
// 			return -1;
// 		}
// 		if (g->nodes[u_rc].deg > 2 || g->nodes[u].deg > 2)
// 			return 0;
// 		// asm_duplicate_edge_seq(g, e, cov);
// 		// asm_duplicate_edge_seq(g, e_rc, cov);
// 		/* split the node */
// 		v = asm_create_node(g);
// 		v_rc = g->nodes[v].rc_id;
// 		/* set source-target of edge e */
// 		g->edges[e].target = v;
// 		asm_remove_node_adj(g, u_rc, e_rc);
// 		g->edges[e_rc].source = v_rc;
// 		asm_add_node_adj(g, v_rc, e_rc);
// 		/* move edge from node u to node v */
// 		g->nodes[v].adj = malloc(g->nodes[u].deg * sizeof(gint_t));
// 		memcpy(g->nodes[v].adj, g->nodes[u].adj, g->nodes[u].deg * sizeof(gint_t));
// 		g->nodes[v].deg = g->nodes[u].deg;
// 		asm_remove_node_adj(g, v, e);
// 		/* node u has only edge e left */
// 		g->nodes[u].adj = realloc(g->nodes[u].adj, sizeof(gint_t));
// 		g->nodes[u].adj[0] = e;
// 		g->nodes[u].deg = 1;
// 		gint_t j;
// 		for (j = 0; j < g->nodes[v].deg; ++j) {
// 			gint_t e_t = g->nodes[v].adj[j];
// 			// fprintf(stderr, "e_t[%ld_%ld] %ld->%ld\n",
// 			// 	e_t, g->edges[e_t].rc_id, g->edges[e_t].source,
// 			// 	g->edges[e_t].target);
// 			g->edges[e_t].source = v;
// 			g->edges[g->edges[e_t].rc_id].target = v_rc;
// 		}
// 		return 1;
// 	} else if (u == v_rc) { /* self loop reverse */
// 		if (cov == 0) {
// 			__VERBOSE("Remove edges %ld[len=%u]\n", e, get_edge_len(g->edges + e));
// 			asm_remove_edge(g, e);
// 			asm_remove_edge(g, e_rc);
// 			return -1;
// 		}
// 		// if (g->nodes[u].deg > 2)
// 		// 	return 0;
// 		// if (g->nodes[v].deg == 1) {
// 		// 	e_rc1 = g->nodes[v].adj[0];
// 		// 	e1 = g->edges[e_rc1].rc_id;
// 		// 	asm_join_edge_loop_reverse(g, e1, e, e_rc, e_rc1);
// 		// 	asm_remove_edge(g, e);
// 		// 	asm_remove_edge(g, e_rc);
// 		// 	return 2;
// 		// } else if (g->nodes[v].deg == 2) {
// 		// 	e_rc1 = g->nodes[v].adj[0];
// 		// 	e1 = g->edges[e_rc1].rc_id;
// 		// 	e2 = g->nodes[v].adj[1];
// 		// 	e_rc2 = g->edges[e2].rc_id;
// 		// 	asm_join_edge3(g, e1, e_rc1, e, e_rc, e2, e_rc2,
// 		// 				g->edges[e].count);
// 		// 	asm_remove_edge(g, e);
// 		// 	asm_remove_edge(g, e_rc);
// 		// 	return 2;
// 		// }
// 		return 0;
// 	} else {
// 		if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
// 			g->nodes[u_rc].deg > 2 || g->nodes[v].deg > 2)
// 			return 0;
// 		e_return = -1;
// 		gint_t j;
// 		for (j = 0; j < g->nodes[v].deg; ++j) {
// 			if (g->edges[g->nodes[v].adj[j]].target == u) {
// 				e_return = g->nodes[v].adj[j];
// 				break;
// 			}
// 		}
// 		if (e_return == -1)
// 			return 0;
// 		e_return_rc = g->edges[e_return].rc_id;
// 		if (g->edges[e].seq_len >= MIN_CONTIG_READPAIR ||
// 			g->edges[e_return].seq_len >= MIN_CONTIG_READPAIR)
// 			return 0;
// 		double fcov_e, fcov_e_return;
// 		fcov_e = __get_edge_cov(g->edges + e, g->ksize) / uni_cov;
// 		fcov_e_return = __get_edge_cov(g->edges + e_return, g->ksize) / uni_cov;
// 		struct cov_range_t rcov_e, rcov_e_return;
// 		rcov_e = convert_cov_range(fcov_e);
// 		rcov_e_return = convert_cov_range(fcov_e_return);
// 		int rep = __min(rcov_e.lo - 1, rcov_e_return.lo);
// 		if (rep <= 0)
// 			rep = 1;
// 		asm_unroll_loop_forward(g, e, e_return, rep);
// 		asm_unroll_loop_forward(g, e_rc, e_return_rc, rep);
// 		asm_remove_edge(g, e_return);
// 		asm_remove_edge(g, e_return_rc);
// 		return 3;
// 	}
// 	return 0;
// }

gint_t unroll_simple_loop(struct asm_graph_t *g)
{
	gint_t e, cnt_self, cnt_self_rv, cnt_double, cnt_false, ret;
	cnt_self = cnt_self_rv = cnt_double = cnt_false = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		ret = check_simple_loop(g, e);
		cnt_self += (ret == 1);
		cnt_self_rv += (ret == 2);
		cnt_double += (ret == 3);
		cnt_false += (ret == -1);
	}
	__VERBOSE("Number of unroll: self loop (%ld), self loop reverse (%ld), double loop (%ld), false loop (%ld)\n",
		cnt_self, cnt_self_rv, cnt_double, cnt_false);
	return cnt_self + cnt_self_rv + cnt_double + cnt_false;
}

static gint_t bubble_keep_longest(struct asm_graph_t *g, gint_t *b, gint_t n)
{
	gint_t e_kept, e, e_rc, i;
	uint32_t max_len = 0;
	uint64_t sum_cnt = 0;
	e_kept = -1;
	for (i = 0; i < n; ++i) {
		e = b[i];
		if (g->edges[e].seq_len > max_len) {
			max_len = g->edges[e].seq_len;
			e_kept = e;
		}
		sum_cnt += g->edges[e].count;
	}
	for (i = 0; i < n; ++i) {
		e = b[i];
		if (e != e_kept) {
			e_rc = g->edges[e].rc_id;
			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
		}
	}
	g->edges[e_kept].count = sum_cnt;
	e_rc = g->edges[e_kept].rc_id;
	g->edges[e_rc].count = sum_cnt;
	return n - 1;
}

static gint_t check_simple_bubble(struct asm_graph_t *g, gint_t se)
{
	gint_t u, v, e, ret, n, j;
	gint_t *branch;
	u = g->edges[se].source;
	v = g->edges[se].target;
	if (u == g->nodes[v].rc_id) /* self loop reverse */
		return 0;
	branch = alloca(g->nodes[u].deg * sizeof(gint_t));
	n = 0;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		if (g->edges[e].target == v && g->edges[e].seq_len < MIN_NOTICE_LEN)
			branch[n++] = e;
	}
	if (n < 2)
		return 0;
	bubble_keep_longest(g, branch, n);
	return n;
}

gint_t resolve_simple_bubble(struct asm_graph_t *g)
{
	gint_t e, cnt_collapsed;
	cnt_collapsed = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		cnt_collapsed += check_simple_bubble(g, e);
	}
	__VERBOSE("Number of collapsed bubble: %ld\n", cnt_collapsed);
	return cnt_collapsed;
}

void resolve_graph_operation(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	gint_t cnt_tips, cnt_tips_complex, cnt_chimeric, cnt_loop, cnt_collapse;
	int iter = 0;
	do {
		__VERBOSE("Iteration [%d]\n", ++iter);
		cnt_tips = cnt_tips_complex = cnt_chimeric = 0;

		cnt_tips = remove_tips(g0);
		asm_condense(g0, g);
		asm_graph_destroy(g0);
		*g0 = *g;

		cnt_tips_complex = remove_tips_topo(g0);
		asm_condense(g0, g);
		asm_graph_destroy(g0);
		*g0 = *g;

		cnt_chimeric = remove_chimeric(g0);
		asm_condense(g0, g);
		asm_graph_destroy(g0);
		*g0 = *g;

		do {
			cnt_loop = cnt_collapse = 0;

			// cnt_loop = unroll_simple_loop(g0);
			cnt_collapse = resolve_simple_bubble(g0);
			asm_lazy_condense(g0);
		} while (cnt_loop + cnt_collapse);

		asm_condense(g0, g);
		asm_graph_destroy(g0);
		*g0 = *g;
	} while (cnt_tips + cnt_tips_complex + cnt_chimeric);
}

