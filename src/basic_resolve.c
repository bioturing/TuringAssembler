#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "io_utils.h"
#include "khash.h"
#include "resolve.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"
#include "graph_search.h"
#include "kmer_hash.h"
#include "process.h"
#define KSIZE_CHECK (g->ksize + 6)
#define LEN_VAR				20
#define MAX_JOIN_LEN			2000
#define __get_edge_cov_int(g, e, uni_cov) (int)((g)->edges[e].count * 1.0 /    \
	((g)->edges[e].seq_len - ((g)->edges[e].n_holes + 1) * (g)->ksize) /   \
	(uni_cov) + 0.499999999)
#define JUNGLE_RADIUS 10
#define MIN_NOTICE_BRIDGE 4000
#define MAX_DUMP_EDGE_LEN 200
#define MIN_COVERAGE_RATIO 0.3
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

void asm_condense_barcode(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	if (!(g0->aux_flag & ASM_HAVE_BARCODE))
		__ERROR("Graph must have barcode\n");
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
		int is_single_loop = 0;
		if (deg_fw == 1 && deg_rv == 1){
			int fw_e = g0->nodes[u].adj[0];
			int rv_e = g0->edges[g0->nodes[g0->nodes[u].rc_id].adj[0]].rc_id;
			if (fw_e == rv_e)
				is_single_loop = 1;
		}
		if (!is_single_loop && ((deg_fw == 1 && deg_rv == 1)
				|| deg_fw + deg_rv == 0 || is_dead_end(g0, u)))
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
			edges[p].barcodes = calloc(3, sizeof(struct barcode_hash_t));
			for (int i = 0; i < 3; ++i){
				barcode_hash_clone(edges[p].barcodes + i,
						g0->edges[e].barcodes + i);
			}
			int n_mid = 0;
			int *mid_nodes = calloc(g0->n_v, sizeof(int));
			do {
				v = g0->edges[e].target;
				mid_nodes[n_mid++] = v;
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
					asm_append_barcode_edge(edges + p,
							g0->edges + e);
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
			edges[q].barcodes = calloc(3, sizeof(struct barcode_hash_t));
			for (int i = 0; i < 3; ++i){
				barcode_hash_clone(edges[q].barcodes + i,
						g0->edges[e_rc].barcodes + i);
			}
			for (int i = n_mid - 2; i >= 0; --i){
				int v = g0->nodes[mid_nodes[i]].rc_id;
				int e_rc = g0->nodes[v].adj[0];
				asm_append_barcode_edge(edges + q, g0->edges + e_rc);
			}
			free(mid_nodes);
			gint_t j = find_adj_idx(g0->nodes[v_rc].adj,
						g0->nodes[v_rc].deg, e_rc);
			assert(j >= 0);
			g0->nodes[v_rc].adj[j] = -1; y_rc = node_id[v_rc];
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
	g->aux_flag = g0->aux_flag;
	g->n_v = n_v;
	g->n_e = n_e;
	g->nodes = nodes;
	g->edges = edges;
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
		int is_single_loop = 0;
		if (deg_fw == 1 && deg_rv == 1){
			int fw_e = g0->nodes[u].adj[0];
			int rv_e = g0->edges[g0->nodes[g0->nodes[u].rc_id].adj[0]].rc_id;
			if (fw_e == rv_e)
				is_single_loop = 1;
		}
		if (!is_single_loop && ((deg_fw == 1 && deg_rv == 1)
				|| deg_fw + deg_rv == 0 || is_dead_end(g0, u)))
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
						log_debug("Middle node degree is not equal to 1");
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
	log_debug("Number of tips remove using graph topology: %ld", cnt_removed);
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
	log_debug("Number of trivial tips removed: %ld", cnt_removed);
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
			g->edges[e].seq_len < CHIMERIC_LEN_THRES) {
			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
			++cnt_removed;
		}
	}
	log_debug("Number of chimeric edge removed: %ld", cnt_removed);
	return cnt_removed;
}

int check_simple_loop(struct asm_graph_t *g, gint_t e)
{
	if (g->edges[e].seq_len >= MIN_NOTICE_LEN)
		return 0;
	gint_t e_rc, u, v, u_rc, v_rc, e1, e2, e_rc1, e_rc2, e_return, e_return_rc, j;
	int rep, rep_e, rep_e_return, n_edges;
	double cov, sum_cov;
	e_rc = g->edges[e].rc_id;
	u = g->edges[e].source;
	v = g->edges[e].target;
	u_rc = g->nodes[u].rc_id;
	v_rc = g->nodes[v].rc_id;
	cov = __get_edge_cov(g->edges + e, g->ksize);
	if (u == v) { /* self loop */
		sum_cov = 0;
		n_edges = 0;
		e1 = e2 = -1;
		for (j = 0; j < g->nodes[u_rc].deg; ++j) {
			if (g->nodes[u_rc].adj[j] != e) {
				e1 = g->nodes[u_rc].adj[j];
				sum_cov += __get_edge_cov(g->edges + e1, g->ksize);
				++n_edges;
			}
		}
		for (j = 0; j < g->nodes[u].deg; ++j) {
			if (g->nodes[u].adj[j] != e) {
				e2 = g->nodes[u].adj[j];
				sum_cov += __get_edge_cov(g->edges + e2, g->ksize);
				++n_edges;
			}
		}
		if (e1 == -1 && e2 == -1)
			return 0;
		if (cov < sum_cov / n_edges * 0.5) {
			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
			return -1;
		}
		if (g->nodes[u_rc].deg > 2 || g->nodes[u].deg > 2)
			return 0;
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
			g->edges[e_t].source = v;
			g->edges[g->edges[e_t].rc_id].target = v_rc;
		}
		return 1;
	} else if (u == v_rc) { /* self loop reverse */
		sum_cov = 0;
		n_edges = 0;
		for (j = 0; j < g->nodes[u_rc].deg; ++j) {
			e1 = g->nodes[u_rc].adj[j];
			sum_cov += __get_edge_cov(g->edges + e1, g->ksize);
			++n_edges;
		}
		for (j = 0; j < g->nodes[u].deg; ++j) {
			if (g->nodes[u].adj[j] != e && g->nodes[u].adj[j] != e_rc) {
				e2 = g->nodes[u].adj[j];
				sum_cov += __get_edge_cov(g->edges + e2, g->ksize);
				++n_edges;
			}
		}
		if (cov < sum_cov / n_edges * 0.5) {
			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
			return -1;
		}
	} else {
		if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
			g->nodes[u_rc].deg > 2 || g->nodes[v].deg > 2)
			return 0;
		e1 = e2 = e_return = e_return_rc = -1;
		for (j = 0; j < g->nodes[v].deg; ++j) {
			if (g->edges[g->nodes[v].adj[j]].target == u) {
				e_return = g->nodes[v].adj[j];
			} else {
				e2 = g->nodes[v].adj[j];
			}
		}
		for (j = 0; j < g->nodes[u_rc].deg; ++j) {
			if (g->edges[g->nodes[u_rc].adj[j]].target == v_rc)
				e_return_rc = g->nodes[u_rc].adj[j];
			else
				e1 = g->nodes[u_rc].adj[j];
		}
		if (e_return == -1 || e_return_rc == -1)
			return 0;
		if (g->edges[e_return].seq_len >= MIN_NOTICE_LEN)
			return 0;
		if (e1 == -1 && e2 == -1)
			return 0;
		double fcov_e, fcov_e_return, mean_cov;
		struct cov_range_t rcov_e, rcov_e_return;
		if (e1 == -1)
			mean_cov = __get_edge_cov(g->edges + e2, g->ksize);
		else if (e2 == -1)
			mean_cov = __get_edge_cov(g->edges + e1, g->ksize);
		else
			mean_cov = (__get_edge_cov(g->edges + e1, g->ksize) +
				__get_edge_cov(g->edges + e2, g->ksize)) / 2;
		fcov_e = __get_edge_cov(g->edges + e, g->ksize) / mean_cov;
		fcov_e_return = __get_edge_cov(g->edges + e_return, g->ksize) / mean_cov;
		rcov_e = convert_cov_range(fcov_e);
		rcov_e_return = convert_cov_range(fcov_e_return);
		rep = __min(rcov_e.lo - 1, rcov_e_return.lo);
		if (rep <= 0)
			rep = 1;
		asm_unroll_loop_forward(g, e, e_return, rep);
		asm_unroll_loop_forward(g, e_rc, e_return_rc, rep);
		asm_remove_edge(g, e_return);
		asm_remove_edge(g, e_return_rc);
		return 3;
	}
	return 0;
}

/* return 0: not at all
 * return 1: self - loop
 * return 2: self - loop - reverse
 * return 3: double - loop
 * return -1: false loop, just go through
 */
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
	log_debug("Number of unroll: self loop (%ld), self loop reverse (%ld), double loop (%ld), false loop (%ld)",
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

static int bubble_check_align_edge(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	uint32_t m, n, i, j, c1, c2;
	m = g->edges[e1].seq_len;
	n = g->edges[e2].seq_len;
	int score, ret;
	int *A = malloc((m + 1) * (n + 1) * sizeof(int));
	A[0] = 0;
	for (i = 1; i <= m; ++i)
		A[i * (n + 1)] = -i * 3;
	for (j = 1; j <= n; ++j)
		A[j] = -j * 3;
	for (i = 1; i <= m; ++i) {
		for (j = 1; j <= n; ++j) {
			c1 = __binseq_get(g->edges[e1].seq, i - 1);
			c2 = __binseq_get(g->edges[e2].seq, j - 1);
			score = (c1 == c2 && c1 < 4) ? 1 : -1;
			A[i * (n + 1) + j] = __max(A[i * (n + 1) + j - 1], A[(i - 1) * (n + 1) + j]) - 3;
			A[i * (n + 1) + j] = __max(A[i * (n + 1) + j], A[(i - 1) * (n + 1) + j - 1] + score);
		}
	}
	ret = (A[(m + 1) * (n + 1) - 1] * 100 > 50 * (int)__max(m, n) && __max(m, n) - A[(m + 1) * (n + 1) - 1] < MIN_NOTICE_LEN * 2);
	free(A);
	return ret;
}

static gint_t check_align_bubble(struct asm_graph_t *g, gint_t se)
{
	gint_t u, v, e, ret, n, j;
	gint_t *branch;
	u = g->edges[se].source;
	v = g->edges[se].target;
	if (u == g->nodes[v].rc_id) /* self loop reverse */
		return 0;
	if (g->edges[se].seq_len >= 1000)
		return 0;
	branch = alloca(g->nodes[u].deg * sizeof(gint_t));
	n = 1;
	branch[0] = se;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		if (g->edges[e].seq_len < 1000 && g->edges[e].target == v &&
			e != se && bubble_check_align_edge(g, se, e))
			branch[n++] = e;
	}
	if (n < 2)
		return 0;
	bubble_keep_longest(g, branch, n);
	return n;
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
	log_debug("Number of collapsed bubble: %ld", cnt_collapsed);
	return cnt_collapsed;
}

gint_t resolve_align_bubble(struct asm_graph_t *g)
{
	gint_t e, cnt_collapsed;
	cnt_collapsed = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		cnt_collapsed += check_align_bubble(g, e);
	}
	log_debug("Number of collapsed aligned bubble: %ld", cnt_collapsed);
	return cnt_collapsed;
}

void resolve_local_graph_operation(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	gint_t cnt_tips, cnt_tips_complex, cnt_chimeric, cnt_loop, cnt_collapse;
	int iter = 0;
	do {
		log_debug("Iteration [%d]", ++iter);
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

			cnt_loop = unroll_simple_loop(g0);
			cnt_collapse = resolve_simple_bubble(g0);
			cnt_collapse += resolve_align_bubble(g0);
			cnt_loop += resolve_loop(g0);
			asm_lazy_condense(g0);
		} while (cnt_loop + cnt_collapse);

		asm_condense(g0, g);
		asm_graph_destroy(g0);
		*g0 = *g;

	} while (cnt_tips + cnt_tips_complex + cnt_chimeric);
}

void resolve_graph_operation(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	gint_t cnt_tips, cnt_tips_complex, cnt_chimeric, cnt_loop, cnt_collapse;
	int iter = 0;
	do {
		log_debug("Iteration [%d]", ++iter);
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

			cnt_loop = unroll_simple_loop(g0);
			cnt_collapse = resolve_simple_bubble(g0);
			cnt_collapse += resolve_align_bubble(g0);
			cnt_loop += resolve_loop(g0);
			asm_lazy_condense(g0);
		} while (cnt_loop + cnt_collapse);

		asm_condense(g0, g);
		asm_graph_destroy(g0);
		*g0 = *g;

	} while (cnt_tips + cnt_tips_complex + cnt_chimeric);
}

int check_loop(struct asm_graph_t *g, int i_e2)
{
	/*
 	 *
 	 * -------a'<------b'------
 	 *        |\      /|   
 	 *        | \    / |   
 	 *        |  \  /  |  
 	 *        |   \/   |     
 	 *        |   /\   |    
 	 *        |  /  \  |       
 	 *        | /    \ |          
 	 *        |/      \|          
 	 *        |v      v|          
 	 * ------>a ------>b ----->
 	 *    e1      e2     e3
 	 * check deg 
 	 * check dep of e2 > e4 leng e4 < 200
 	 * leng e1 > 1000 leng e3 > 1000
 	 */
	struct asm_edge_t *e2 = &g->edges[i_e2];
	int i_a = e2->source, i_b = e2->target;
	struct asm_node_t *a = &g->nodes[i_a], *b = &g->nodes[i_b];
	int i_a_rc = a->rc_id, i_b_rc = b->rc_id;
	struct asm_node_t *a_rc = &g->nodes[i_a_rc], *b_rc = &g->nodes[i_b_rc];
	if (a->deg != 1)
		return 0;
	if (b->deg != 1)
		return 0;
	if (a_rc->deg != 2) 
		return 0;
	if (b_rc->deg != 2)
		return 0;
	int b1 = 0;
	int i_e4 = 0, i_e1 = 0, i_e3 = 0;
	struct asm_edge_t *e1 = NULL, *e3 = NULL, *e4 = NULL;
	for (int i = 0; i < 2; i++) {
		struct asm_edge_t *e = &g->edges[a_rc->adj[i]];
		if (e->target != i_b)  {
			i_e1 = e->rc_id;
		} else {
			b1 = 1;
		}
	}
	if (b1 == 0) 
		return 0;
	i_e3 = b->adj[0];
	for (int i = 0; i < 2; i++) {
		struct asm_edge_t *e = &g->edges[b_rc->adj[i]];
		if (e->target == i_a) {
			i_e4 = b_rc->adj[i];
		} else {
			if (e->target != i_a_rc)
				return 0;
		}
	}
	e1 = &g->edges[i_e1];
	e3 = &g->edges[i_e3];
	e4 = &g->edges[i_e4];
//	if (e1->seq_len < 1000)
//		return 0;
//	if (e3->seq_len < 1000)
//		return 0;
	float cov_e2 = __get_edge_cov(e2, g->ksize);
	float cov_e4 = __get_edge_cov(e4, g->ksize);
	log_debug("cov e2 %f e4 %f e4len %d", cov_e2, cov_e4, e4->seq_len);
	if (cov_e2 < cov_e4)
		return 0;
	if (e4->seq_len > 200)
		return 0;
	log_debug("check cov ok");
	asm_remove_edge(g, i_e4);
	int i_e4_rc = g->edges[i_e4].rc_id;
	asm_remove_edge(g, i_e4_rc);
	return 1;
}

int resolve_loop(struct asm_graph_t *g0)
{
	int count = 0;
	for (int i_e2 = 0; i_e2 < g0->n_e; i_e2++) if (g0->edges[i_e2].source != -1) {
		count += check_loop(g0, i_e2);
	}
	log_debug("remove %d loop", count);
	return count;
}

int asm_resolve_dump_loop_ite(struct asm_graph_t *g)
{
	int ite = 0;
	int res = 0;
	do{
		int resolved = asm_resolve_dump_loop(g);
		if (!resolved)
			break;
		res += resolved;
		++ite;
		log_debug("%d-th iteration: %d loop(s) resolved", ite, resolved);
	} while(1);
	log_info("%d dump loop(s) resolved after %d iterations", res, ite);
	return res;
}

int asm_resolve_dump_loop(struct asm_graph_t *g)
{
	int res = 0;
	int tmp_n_e = g->n_e;
	for (int e = 0; e < tmp_n_e; ++e){
		int rc = g->edges[e].rc_id;
		if (e > rc)
			continue;
		int tg = g->edges[e].target;
		if (tg == -1)
			continue;
		int sr = g->nodes[g->edges[e].source].rc_id;
		if (g->nodes[tg].deg == 2 && g->nodes[sr].deg == 2){
			int loop_e = -1;
			for (int i = 0; loop_e == -1 && i < 2; ++i){
				for (int j = 0; loop_e == -1 && j < 2; ++j){
					if (g->nodes[tg].adj[i] ==
						g->edges[g->nodes[sr].adj[j]].rc_id)
						loop_e = g->nodes[tg].adj[i];
				}
			}
			if (loop_e == -1)
				continue;
			int e1 = g->edges[g->nodes[sr].adj[0]].rc_id != loop_e ?
				g->edges[g->nodes[sr].adj[0]].rc_id :
				g->edges[g->nodes[sr].adj[1]].rc_id;
			int e2 = g->nodes[tg].adj[0] != loop_e ?
				g->nodes[tg].adj[0] : g->nodes[tg].adj[1];
			if (e1 == e2 || e == loop_e)
				continue;
			log_debug("Dump loop detected, e1: %d, e: %d, loop e: %d, e2: %d",
					e1, e, loop_e, e2);
			asm_append_seq(g->edges + loop_e, g->edges + e,
					g->ksize);
			asm_append_barcode_readpair(g, loop_e, e);
			asm_append_seq(g->edges + e, g->edges + loop_e,
					g->ksize);
			asm_append_barcode_readpair(g, e, loop_e);
			g->edges[e].count += g->edges[e].count + g->edges[loop_e].count;
			int loop_e_rc = g->edges[loop_e].rc_id;
			int e_rc = g->edges[e].rc_id;
			asm_append_seq(g->edges + loop_e_rc, g->edges + e_rc,
					g->ksize);
			asm_append_barcode_readpair(g, loop_e_rc, e_rc);
			asm_append_seq(g->edges + e_rc, g->edges + loop_e_rc,
					g->ksize);
			asm_append_barcode_readpair(g, e_rc, loop_e_rc);
			g->edges[e_rc].count = g->edges[e].count;
			asm_remove_edge(g, loop_e);
			asm_remove_edge(g, loop_e_rc);
			++res;
		}
	}
	test_asm_graph(g);
	return res;
}

int asm_resolve_dump_branch(struct asm_graph_t *g)
{
	int res = 0;
	for (int e = 0; e < g->n_e; ++e){
		int rc = g->edges[e].rc_id;
		if (e > rc)
			continue;
		int tg = g->edges[e].target;
		if (tg == -1)
			continue;
		if (g->nodes[tg].deg != 2)
			continue;
		int next_edge[2] = {-1, -2};
		int mid_edge[2];
		for (int i = 0; i < 2; ++i){
			int mid_e = g->nodes[tg].adj[i];
			mid_edge[i] = mid_e;
			int mid_tg = g->edges[mid_e].target;
			if (g->nodes[mid_tg].deg != 1)
				break;
			next_edge[i] = g->nodes[mid_tg].adj[0];
		}
		if (next_edge[0] != next_edge[1] || next_edge[0] == e)
			continue;
		__VERBOSE("Dump branch detected, e1: %d, e2: %d, branches: %d %d\n",
				e, next_edge[0], mid_edge[0], mid_edge[1]);
		int trash_e = __get_edge_cov(g->edges + mid_edge[0], g->ksize)
			< __get_edge_cov(g->edges + mid_edge[1], g->ksize) ?
			mid_edge[0] : mid_edge[1];
		asm_remove_edge(g, trash_e);
		trash_e = g->edges[trash_e].rc_id;
		asm_remove_edge(g, trash_e);
		++res;
	}
	struct asm_graph_t g1;
	asm_condense(g, &g1);
	asm_graph_destroy(g);
	*g = g1;
	return res;
}

int asm_resolve_dump_jungle_ite(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	int ite = 0;
	int res = 0;
	do{
		int resolved = asm_resolve_dump_jungle(opt, g);
		if (!resolved)
			break;
		res += resolved;
		++ite;
		log_debug("%d-th iteration: %d jungle(s) resolved", ite, resolved);
		/*char graph[1024]; // Save graph for debugging
		sprintf(graph, "level_pro_%d_ite", ite);
		save_graph_info(opt->out_dir, g, graph);
		graph[0] = '\0';
		sprintf(graph, "%s/graph_k_%d_level_pro_%d_ite.bin", opt->out_dir,
				g->ksize, ite);
		save_asm_graph(g, graph);*/
	} while(1);
	log_info("%d dump jungle(s) resolved after %d iterations",
			res, ite);
	return res;
}

int detect_dump_jungle(struct asm_graph_t *g, int e, int **dump_edges, int *n_dump)
{
	*dump_edges = NULL;
	*n_dump = 0;
	//			STAGE 1 			//
	int *nearby;
	int n_nb;
	struct graph_info_t ginfo;
	graph_info_init(g, &ginfo, e, e);
	get_nearby_edges(g, e, &ginfo, JUNGLE_RADIUS, &nearby, &n_nb);
	graph_info_destroy(&ginfo);

	int e1 = e;
	int e2 = -1;
	for (int i = 0; i < n_nb; ++i){
		if (nearby[i] == e || nearby[i] == g->edges[e].rc_id)
			continue;
		if (g->edges[nearby[i]].seq_len >= MIN_NOTICE_BRIDGE){
			e2 = nearby[i];
			break;
		}
	}
	if (e2 == -1)
		goto free_stage_1;

	//			STAGE 2 			//
	graph_info_init(g, &ginfo, e, e);
	mark_edge_trash(&ginfo, e1);
	mark_edge_trash(&ginfo, g->edges[e1].rc_id);
	mark_edge_trash(&ginfo, e2);
	mark_edge_trash(&ginfo, g->edges[e2].rc_id);
	int *nb1, n_nb1;
	get_nearby_edges(g, e1, &ginfo, JUNGLE_RADIUS, &nb1, &n_nb1);
	graph_info_destroy(&ginfo);

	graph_info_init(g, &ginfo, e, e);
	mark_edge_trash(&ginfo, e1);
	mark_edge_trash(&ginfo, g->edges[e1].rc_id);
	mark_edge_trash(&ginfo, e2);
	mark_edge_trash(&ginfo, g->edges[e2].rc_id);
	int *nb2, n_nb2;
	get_nearby_edges(g, g->edges[e2].rc_id, &ginfo, JUNGLE_RADIUS, &nb2, &n_nb2);
	graph_info_destroy(&ginfo);

	for (int i = 0; i < n_nb1; ++i){
		int e = nb1[i];
		if (e == e1 || e == g->edges[e1].rc_id
				|| e == e2 || e == g->edges[e2].rc_id)
			continue;
		if (g->edges[e].seq_len >= MAX_DUMP_EDGE_LEN){
			e2 = -1;
			break;
		}
	}
	if (e2 == -1)
		goto free_stage_2;

	for (int i = 0; i < n_nb2; ++i){
		int e = g->edges[nb2[i]].rc_id;
		if (e == e1 || e == g->edges[e1].rc_id
				|| e == e2 || e == g->edges[e2].rc_id)
			continue;
		if (g->edges[e].seq_len >= MAX_DUMP_EDGE_LEN){
			e2 = -1;
			break;
		}
	}
	if (e2 == -1)
		goto free_stage_2;


	//			STAGE 3 			//
	int *mark = calloc(g->n_e, sizeof(int));
	for (int i = 0; i < n_nb1; ++i)
		mark[nb1[i]] = 1;
	for (int i = 0; i < n_nb1 && e2 != -1; ++i){
		int e = nb1[i];
		int tg = g->edges[e].target;
		for (int j = 0; j < g->nodes[tg].deg; ++j){
			int next_e = g->nodes[tg].adj[j];
			if (next_e == e1 || next_e == g->edges[e1].rc_id
				|| next_e == e2 || next_e == g->edges[e2].rc_id)
				continue;
			if (!mark[next_e]){
				e2 = -1;
				break;
			}
		}
	}
	if (e2 == -1)
		goto free_stage_3;

	memset(mark, 0, g->n_e * sizeof(int));
	for (int i = 0; i < n_nb2; ++i)
		mark[nb2[i]] = 1;
	for (int i = 0; i < n_nb2 && e2 != -1; ++i){
		int e = nb2[i];
		int tg = g->edges[e].target;
		for (int j = 0; j < g->nodes[tg].deg; ++j){
			int next_e = g->nodes[tg].adj[j];
			if (next_e == e1 || next_e == g->edges[e1].rc_id
				|| next_e == e2 || next_e == g->edges[e2].rc_id)
				continue;
			if (!mark[next_e]){
				e2 = -1;
				break;
			}
		}
	}
	if (e2 == -1)
		goto free_stage_3;
	log_debug("Dump jungle detected at %d and %d", e1, e2);
	*dump_edges = calloc(n_nb1 + n_nb2 - 2, sizeof(int));
	for (int i = 1; i < n_nb1; ++i)
		(*dump_edges)[(*n_dump)++] = nb1[i];
	for (int i = 1; i < n_nb2; ++i)
		(*dump_edges)[(*n_dump)++] = g->edges[nb2[i]].rc_id;

free_stage_3:
	free(mark);
free_stage_2:
	free(nb1);
	free(nb2);
free_stage_1:
	free(nearby);
	return e2;
}

int asm_resolve_dump_jungle(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	int res = 0;
	struct read_path_t read_sorted_path;
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
		log_error("Reads must be sorted");
	}
	khash_t(bcpos) *dict = kh_init(bcpos);
	construct_read_index(&read_sorted_path, dict);
	int tmp = g->n_e;
	int m_e = g->n_e;
	for (int e1 = 0; e1 < tmp; ++e1){
		if (g->edges[e1].target == -1)
			continue;
		if (g->edges[e1].seq_len < MIN_NOTICE_BRIDGE)
			continue;
		//		STAGE 1 detect dump jungle and find path  		//
		int *dump_edges;
		int n_dump;
		int e2 = detect_dump_jungle(g, e1, &dump_edges, &n_dump);
		if (e2 == -1)
			continue;
		struct path_info_t pinfo;
		path_info_init(&pinfo);
		struct edge_map_info_t emap1;
		struct edge_map_info_t emap2;
		emap1.lc_e = e1;
		emap2.lc_e = e2;

		log_debug("Get local reads");
		struct read_path_t local_read_path;
		get_union_barcode_reads(opt, g, e1, e2, dict, &read_sorted_path,
				&local_read_path);
		khash_t(kmer_int) *kmer_count = get_kmer_hash(local_read_path.R1_path,
				local_read_path.R2_path, KSIZE_CHECK);
		log_debug("Finding paths between %d and %d", e1, e2);
		get_all_paths_kmer_check(g, &emap1, &emap2, &pinfo, KSIZE_CHECK,
				kmer_count);
		int longest_path = 0;
		for (int i = 0; i < pinfo.n_paths; ++i){
			if (pinfo.path_lens[i] > pinfo.path_lens[longest_path])
				longest_path = i;
		}
		int *path = pinfo.paths[longest_path];
		int len = pinfo.path_lens[longest_path];
		if (len <= 2){
			log_debug("No reliable paths found, continue");
			goto ignore_stage_1;
		}

		//		STAGE 2 condence the jungle 		//
		if (g->n_e == m_e){
			g->edges = (struct asm_edge_t *) realloc(g->edges,
					sizeof(struct asm_edge_t) * (m_e << 1));
			memset(g->edges + m_e, 0, m_e * sizeof(struct asm_edge_t));
			m_e <<= 1;
		}
		int p = g->n_e;
		int q = g->n_e + 1;
		asm_clone_edge(g, p, path[0]);
		for (int i = 1; i < len; ++i){
			int e1 = p;
			int e2 = path[i];
			asm_append_seq(g->edges + e1, g->edges + e2, g->ksize);
			asm_append_barcode_readpair(g, e1, e2);
			g->edges[e1].count += g->edges[e2].count;
		}

		asm_clone_edge(g, q, g->edges[path[len - 1]].rc_id);
		for (int i = len - 2; i >= 0; --i){
			int e1 = q;
			int e2 = g->edges[path[i]].rc_id;
			asm_append_seq(g->edges + e1, g->edges + e2, g->ksize);
			asm_append_barcode_readpair(g, e1, e2);
			g->edges[e1].count += g->edges[e2].count;
		}

		int u = g->edges[path[0]].source;
		int v = g->edges[path[len - 1]].target;
		g->edges[p].source = u;
		g->edges[p].target = v;
		g->edges[p].rc_id = q;
		g->nodes[u].adj = realloc(g->nodes[u].adj, (g->nodes[u].deg + 1)
				* sizeof(struct asm_edge_t));
		g->nodes[u].adj[g->nodes[u].deg] = p;
		++g->nodes[u].deg;
		int u_rc = g->nodes[u].rc_id;
		int v_rc = g->nodes[v].rc_id;
		g->edges[q].source = v_rc;
		g->edges[q].target = u_rc;
		g->edges[q].rc_id = p;
		g->nodes[v_rc].adj = realloc(g->nodes[v_rc].adj, (g->nodes[v_rc].deg + 1)
				* sizeof(struct asm_edge_t));
		g->nodes[v_rc].adj[g->nodes[v_rc].deg] = q;
		++g->nodes[v_rc].deg;
		g->n_e += 2;

		for (int i = 0; i < n_dump; ++i){
			int e = dump_edges[i];
			int e_rc = g->edges[e].rc_id;
			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
		}
		asm_remove_edge(g, e1);
		asm_remove_edge(g, g->edges[e1].rc_id);
		asm_remove_edge(g, e2);
		asm_remove_edge(g, g->edges[e2].rc_id);
		++res;
ignore_stage_1:
		path_info_destroy(&pinfo);
		kh_destroy(kmer_int, kmer_count);
		destroy_read_path(&local_read_path);
		free(dump_edges);
	}
	kh_destroy(bcpos, dict);
	test_asm_graph(g);
	return res;
}

