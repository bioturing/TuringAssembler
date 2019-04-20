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

#define __positive_ratio(r)		((r) + EPS >= 0.02)
#define MAX_EDGE_COUNT			10000
#define __get_edge_cov_int(g, e, uni_cov) (int)((g)->edges[e].count * 1.0 /    \
	((g)->edges[e].seq_len - ((g)->edges[e].n_holes + 1) * (g)->ksize) /   \
	(uni_cov) + 0.499999999)
#define __is_hang_edge(g, set_v, e) (kh_get(gint, set_v, (g)->edges[e].target) == kh_end(set_v))

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
	uint32_t l, r;
	l = r = 0;
	uint32_t edge_count = 1;
	q[0] = se;
	cb_add(set_e, se);
	cb_add(set_e, g->edges[se].rc_id);
	cb_add(set_v, g->edges[se].target);
	cb_add(set_v, g->nodes[g->edges[se].target].rc_id);
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
					if (r + 1 == max_edge_count)
						goto clean_up;
					q[++r] = ne;
				}
				if (cb_add(set_e, ne_rc) == 1) {
					if (r + 1 == max_edge_count)
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
			}
		}
		v = g->nodes[v].rc_id;
		for (j = 0; j < g->nodes[v].deg; ++j) {
			gint_t ne = g->nodes[v].adj[j], ne_rc;
			ne_rc = g->edges[ne].rc_id;
			uint32_t len = get_edge_len(g->edges + ne);
			if (len < min_contig_len) {
				if (cb_add(set_e, ne) == 1) {
					if (r + 1 == max_edge_count)
						goto clean_up;
					q[++r] = ne;
				}
				if (cb_add(set_e, ne_rc) == 1) {
					if (r + 1 == max_edge_count)
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

void detect_leg(struct asm_graph_t *g, uint32_t min_contig_len, khash_t(gint) *set_v,
	khash_t(gint) *set_e, khash_t(gint) *set_leg, khash_t(gint) *set_self)
{
	khiter_t k;
	gint_t e;
	int ret;
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (!kh_exist(set_e, k))
			continue;
		e = kh_key(set_e, k);
		if (__is_hang_edge(g, set_v, e))
			kh_put(gint, set_leg, e, &ret);
	}
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (!kh_exist(set_e, k))
			continue;
		e = kh_key(set_e, k);
		if (kh_get(gint, set_leg, e) == kh_end(set_leg) &&
			kh_get(gint, set_leg, g->edges[e].rc_id) == kh_end(set_leg) &&
			get_edge_len(g->edges + e) >= min_contig_len)
			kh_put(gint, set_self, e, &ret);
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

gint_t dfs_get_distance(struct asm_graph_t *g, uint32_t min_contig_size,
		khash_t(gint) *set_e, khash_t(gint) *vis, gint_t u, gint_t t)
{
	gint_t j, e, v;
	gint_t ret;
	int sret;
	uint32_t len;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		len = get_edge_len(g->edges + e);
		if (kh_get(gint, set_e, e) == kh_end(set_e) ||
			len >= min_contig_size)
			continue;
		v = g->edges[e].target;
		if (v == t)
			return 0;
		if (kh_get(gint, vis, v) == kh_end(vis)) {
			kh_put(gint, vis, v, &sret);
			ret = dfs_get_distance(g, min_contig_size, set_e,
								vis, v, t);
			if (ret != -1)
				return ret + (len - g->ksize);
		}
	}
	return -1;
}

gint_t get_distance(struct asm_graph_t *g, uint32_t min_contig_size,
		khash_t(gint) *set_e, gint_t s, gint_t e)
{
	if (s == e)
		return 0;
	int sret;
	gint_t ret;
	khash_t(gint) *vis;
	vis = kh_init(gint);
	kh_put(gint, vis, s, &sret);
	ret = dfs_get_distance(g, min_contig_size, set_e, vis, s, e);
	kh_destroy(gint, vis);
	return ret;
}

gint_t bc_find_pair_check_path(struct asm_graph_t *g, khash_t(gint) *set_e,
					gint_t se, gint_t *adj, gint_t deg)
{
	gint_t j, e, se_rc, ret_e;
	se_rc = g->edges[se].rc_id;
	double ret_ratio = -1;
	ret_e = -1;
	for (j = 0; j < deg; ++j) {
		e = adj[j];
		if (e == se || e == se_rc)
			continue;
		uint32_t dist = get_distance(g, MIN_CONTIG_BARCODE, set_e,
			g->nodes[g->edges[e].source].rc_id, g->edges[se].source);
		if (dist == -1)
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

gint_t bc_find_pair_check_path_strict(struct asm_graph_t *g, khash_t(gint) *set_e,
					gint_t se, gint_t *adj, gint_t deg)
{
	gint_t j, e, se_rc, ret_e, second_e;
	se_rc = g->edges[se].rc_id;
	double ret_ratio = -1, second_best = -1;
	ret_e = second_e = -1;
	for (j = 0; j < deg; ++j) {
		e = adj[j];
		if (e == se || e == se_rc)
			continue;
		uint32_t dist = get_distance(g, MIN_CONTIG_BARCODE, set_e,
			g->nodes[g->edges[e].source].rc_id, g->edges[se].source);
		if (dist == -1)
			continue;
		double ratio = get_barcode_ratio(g, se, e);
		if (ratio > ret_ratio) {
			second_best = ret_ratio;
			second_e = ret_e;
			ret_ratio = ratio;
			ret_e = e;
		} else if (ratio > second_best) {
			second_best = ratio;
			second_e = e;
		}
	}
	if (!__positive_ratio(ret_ratio))
		return -1;
	if (second_best != -1 && !__strictly_greater(ret_ratio, second_best)) {
		__VERBOSE("Best = %ld(~%.6lf); second best = %ld(~%.6lf)\n",
			ret_e, ret_ratio, second_e, second_best);
		return -2;
	}
	return ret_e;
}

gint_t bc_find_pair_strict(struct asm_graph_t *g, gint_t se, gint_t *adj,
								gint_t deg)
{
	gint_t j, e, se_rc, ret_e, second_e;
	se_rc = g->edges[se].rc_id;
	double ret_ratio = -1, second_best = -1;
	ret_e = -1;
	for (j = 0; j < deg; ++j) {
		e = adj[j];
		if (e == se || e == se_rc)
			continue;
		double ratio = get_barcode_ratio(g, se, e);
		if (ratio > ret_ratio) {
			second_best = ret_ratio;
			second_e = ret_e;
			ret_ratio = ratio;
			ret_e = e;
		} else if (ratio > second_best) {
			second_best = ratio;
			second_e = e;
		}
	}
	if (!__positive_ratio(ret_ratio))
		return -1;
	if (second_best != -1 && !__strictly_greater(ret_ratio, second_best)) {
		__VERBOSE("Best = %ld(~%.6lf); second best = %ld(~%.6lf)\n",
			ret_e, ret_ratio, second_e, second_best);
		return -2;
	}
	return ret_e;
}

gint_t bc_find_pair(struct asm_graph_t *g, gint_t se,
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

double callibrate_uni_cov(struct asm_graph_t *g, gint_t *legs, gint_t n_leg,
								double uni_cov)
{
	double sum_cov = 0, cov, ratio, ret;
	gint_t cnt = 0, i, e;
	for (i = 0; i < n_leg; ++i) {
		e = legs[i];
		cov = __get_edge_cov(g->edges + e, g->ksize);
		ratio = cov / uni_cov;
		if (ratio > 0.75 && ratio < 1.25) {
			sum_cov += cov;
			++cnt;
		}
	}
	if (cnt)
		ret = sum_cov / cnt;
	else
		ret = uni_cov;
//	__VERBOSE("Global cov ~ %.6lf. Local cov ~ %.6lf\n", uni_cov, ret);
	return ret;
}

gint_t check_2_2_small_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t e_rc, v, v_rc, u, u_rc;
	int i, k, flag;
	e_rc = g->edges[e].rc_id;
	v = g->edges[e].target;
	v_rc = g->nodes[v].rc_id;
	u = g->edges[e].source;
	u_rc = g->nodes[u].rc_id;
	if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
		g->nodes[u_rc].deg != 2 || g->nodes[v].deg != 2)
		return 0;

	/* callibrate uni_coverage */
	gint_t *legs = alloca(4 * sizeof(gint_t));
	legs[0] = g->nodes[u_rc].adj[0];
	legs[1] = g->nodes[u_rc].adj[1];
	legs[2] = g->nodes[v].adj[0];
	legs[3] = g->nodes[v].adj[1];
	double *fcov = alloca(4 * sizeof(double));
	struct cov_range_t *rcov = alloca(4 * sizeof(struct cov_range_t));
	for (i = 0; i < 4; ++i) {
		fcov[i] = __get_edge_cov(g->edges + legs[i], g->ksize) / uni_cov;
		rcov[i] = convert_cov_range(fcov[i]);
	}
	for (i = 0; i < 2; ++i) {
		flag = 0;
		for (k = 2; k < 4; ++k) {
			if (__cov_range_intersect(rcov[i], rcov[k]))
				flag = 1;
		}
		if (!flag)
			return 0;
	}
	double *ratio = alloca(4 * sizeof(double));
	ratio[0] = get_barcode_ratio_small(g, legs[0], legs[2]);
	ratio[1] = get_barcode_ratio_small(g, legs[0], legs[3]);
	ratio[2] = get_barcode_ratio_small(g, legs[1], legs[2]);
	ratio[3] = get_barcode_ratio_small(g, legs[1], legs[3]);
	if (__strictly_greater(ratio[0], ratio[1]) && __strictly_greater(ratio[3], ratio[2])) {
		if (!__positive_ratio(ratio[0]) || !__positive_ratio(ratio[3]))
			return 0;
		if (!__cov_range_intersect(rcov[0], rcov[2]) ||
			!__cov_range_intersect(rcov[1], rcov[3]) ||
			!__diff_accept(fcov[0], fcov[2]) ||
			!__diff_accept(fcov[1], fcov[3]))
			return 0;
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 1;
	} else if (__strictly_greater(ratio[1], ratio[0]) && __strictly_greater(ratio[2], ratio[3])) {
		if (!__positive_ratio(ratio[1]) || !__positive_ratio(ratio[2]))
			return 0;
		if (!__cov_range_intersect(rcov[0], rcov[3]) ||
			!__cov_range_intersect(rcov[1], rcov[2]) ||
			!__diff_accept(fcov[0], fcov[3]) ||
			!__diff_accept(fcov[1], fcov[2]))
			return 0;
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 1;
	} else {
		return 0;
	}
}

gint_t check_2_2_strict_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t e_rc, v, v_rc, u, u_rc;
	int i, k, flag;
	e_rc = g->edges[e].rc_id;
	v = g->edges[e].target;
	v_rc = g->nodes[v].rc_id;
	u = g->edges[e].source;
	u_rc = g->nodes[u].rc_id;
	if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
		g->nodes[u_rc].deg != 2 || g->nodes[v].deg != 2)
		return 0;

	/* callibrate uni_coverage */
	gint_t *legs = alloca(4 * sizeof(gint_t));
	legs[0] = g->nodes[u_rc].adj[0];
	legs[1] = g->nodes[u_rc].adj[1];
	legs[2] = g->nodes[v].adj[0];
	legs[3] = g->nodes[v].adj[1];
	for (i = 0; i < 4; ++i) {
		if (get_edge_len(g->edges + legs[i]) < MIN_CONTIG_BARCODE)
			return 0;
	}
	double uni_cov_local = callibrate_uni_cov(g, legs, 4, uni_cov);
	double *fcov = alloca(4 * sizeof(double));
	struct cov_range_t *rcov = alloca(4 * sizeof(struct cov_range_t));
	for (i = 0; i < 4; ++i) {
		fcov[i] = __get_edge_cov(g->edges + legs[i], g->ksize) / uni_cov_local;
		rcov[i] = convert_cov_range(fcov[i]);
	}
	for (i = 0; i < 2; ++i) {
		flag = 0;
		for (k = 2; k < 4; ++k) {
			if (__cov_range_intersect(rcov[i], rcov[k]))
				flag = 1;
		}
		if (!flag)
			return 0;
	}
	double *ratio = alloca(4 * sizeof(double));
	ratio[0] = get_barcode_ratio(g, legs[0], legs[2]);
	ratio[1] = get_barcode_ratio(g, legs[0], legs[3]);
	ratio[2] = get_barcode_ratio(g, legs[1], legs[2]);
	ratio[3] = get_barcode_ratio(g, legs[1], legs[3]);
	if (__strictly_greater(ratio[0], ratio[1]) && __strictly_greater(ratio[3], ratio[2])) {
		if (!__positive_ratio(ratio[0]) || !__positive_ratio(ratio[3]))
			return 0;
		if (!__cov_range_intersect(rcov[0], rcov[2]) ||
			!__cov_range_intersect(rcov[1], rcov[3]) ||
			!__diff_accept(fcov[0], fcov[2]) ||
			!__diff_accept(fcov[1], fcov[3]))
			return 0;
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 1;
	} else if (__strictly_greater(ratio[1], ratio[0]) && __strictly_greater(ratio[2], ratio[3])) {
		if (!__positive_ratio(ratio[1]) || !__positive_ratio(ratio[2]))
			return 0;
		if (!__cov_range_intersect(rcov[0], rcov[3]) ||
			!__cov_range_intersect(rcov[1], rcov[2]) ||
			!__diff_accept(fcov[0], fcov[3]) ||
			!__diff_accept(fcov[1], fcov[2]))
			return 0;
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 1;
	} else {
		return 0;
	}
}

gint_t check_n_m_bridge_strict(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t u, v, v_rc, u_rc, e1, e2, et1, e_rc, i;
	e_rc = g->edges[e].rc_id;
	v = g->edges[e].target;
	v_rc = g->nodes[v].rc_id;
	u = g->edges[e].source;
	u_rc = g->nodes[u].rc_id;
	if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
		(g->nodes[u_rc].deg < 2 && g->nodes[v].deg < 2))
		return 0;

	/* Callibrate uni_coverage */
	gint_t n_leg = g->nodes[u_rc].deg + g->nodes[v].deg + 1;
	gint_t *legs = alloca(n_leg * sizeof(gint_t));
	n_leg = 0;
	for (i = 0; i < g->nodes[u_rc].deg; ++i)
		legs[n_leg++] = g->nodes[u_rc].adj[i];
	for (i = 0; i < g->nodes[v].deg; ++i)
		legs[n_leg++] = g->nodes[v].adj[i];
	legs[n_leg++] = e;
	double uni_cov_local, e_cov, e_uni_cov, sum_fcov, fcov1, fcov2;
	struct cov_range_t rcov1, rcov2;
	uni_cov_local = callibrate_uni_cov(g, legs, n_leg, uni_cov);
	sum_fcov = 0;
	for (i = 0; i < g->nodes[u_rc].deg; ++i) {
		e1 = g->nodes[u_rc].adj[i];
		sum_fcov += __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
	}
	for (i = 0; i < g->nodes[v].deg; ++i) {
		e1 = g->nodes[v].adj[i];
		sum_fcov += __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
	}
	e_cov = __get_edge_cov(g->edges + e, g->ksize) / uni_cov_local;
	e_uni_cov = e_cov / sum_fcov;

	gint_t resolve, ret = 0;
	do {
		resolve = 0;
		for (i = 0; i < g->nodes[u_rc].deg; ++i) {
			e1 = g->nodes[u_rc].adj[i];
			fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
			rcov1 = convert_cov_range(fcov1);
			if (rcov1.hi != 1)
				continue;
			e2 = bc_find_pair_strict(g, e1, g->nodes[v].adj, g->nodes[v].deg);
			if (e2 < 0)
				continue;
			fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
			rcov2 = convert_cov_range(fcov2);
			if (rcov2.hi != 1)
				continue;
			et1 = bc_find_pair_strict(g, e2, g->nodes[u_rc].adj, g->nodes[u_rc].deg);
			if (et1 != e1) {
				__VERBOSE("Not best pair: (%ld, %ld) <-> %ld\n", e1, et1, e2);
				continue;
			}
			asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
				e2, g->edges[e2].rc_id, g->edges[e].count * (fcov1 + fcov2) / sum_fcov);
			e_cov -= e_uni_cov * (fcov1 + fcov2);
			resolve = 1;
			break;
		}
		ret += resolve;
	} while (resolve);
	g->edges[e].count = g->edges[e_rc].count = e_cov / e_uni_cov / sum_fcov * g->edges[e].count;
	/* special case when resolved all but 1 pair left */
	if (g->nodes[u_rc].deg == 1 && g->nodes[v].deg == 1) {
		e1 = g->nodes[u_rc].adj[0];
		e2 = g->nodes[v].adj[0];
		fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
		fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
		rcov1 = convert_cov_range(fcov1);
		rcov2 = convert_cov_range(fcov2);
		double ratio = get_barcode_ratio(g, e1, e2);
		/* FIXME: may need more complicated resolve here */
		__VERBOSE("Leftover path: %ld(~%.6lf) -> %ld -> %ld(~%.6lf), ratio: %.6lf\n",
				e1, fcov1, e, e2, fcov2, ratio);
		if ((__positive_ratio(ratio) || ratio < 0) &&
			__cov_range_intersect(rcov1, rcov2) &&
			__diff_accept(fcov1, fcov2)) {
			if (e1 != g->edges[e2].rc_id)
				// asm_join_edge_loop(g, e1, e2, e, e_rc,
				// 		g->edges[e].count);
			// else
				asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
					e2, g->edges[e2].rc_id, g->edges[e].count);
			++ret;
		}
		if (ret) {
			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
		}
	} else if (g->nodes[u_rc].deg == 0 && g->nodes[v].deg == 0) {
		/* destroy the bridge */
		if (ret) {
			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
		}
	}
	return ret;
}

// gint_t check_n_m_node_strict(struct asm_graph_t *g, gint_t *legs, gint_t n_leg,
// 								double uni_cov)
// {
// 	double uni_cov_local = callibrate_uni_cov(g, legs, n_leg, uni_cov);
// 	gint_t i, e1, e2, et1, ret, resolve;
// 	double fcov1, fcov2;
// 	struct cov_range_t rcov1, rcov2;
// 	ret = 0;
// 	do {
// 		resolve = 0;
// 		for (i = 0; i < n_leg; ++i) {
// 			e1 = legs[i];
// 			fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
// 			rcov1 = convert_cov_range(fcov1);
// 			if (rcov1.hi != 1)
// 				continue;
// 			e2 = bc_find_consecutive_strict(g, e1, legs, n_leg);
// 			if (e2 == -1)
// 				continue;
// 			fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
// 			rcov2 = convert_cov_range(fcov2);
// 			if (rcov2.hi != 1)
// 				continue;
// 			et1 = bc_find_consecutive_strict(g, e2, legs, n_leg);
// 			if (et1 != e1) {
// 				__VERBOSE("Not best pair (%ld, %ld) -> %ld\n",
// 					e1, et1, e2);
// 				continue;
// 			}
// 			asm_join_edge(g, g->edges[e1].rc_id, e1,
// 							e2, g->edges[e2].rc_id);
// 			gint_t j;
// 			j = find_adj_idx(legs, n_leg, e1);
// 			if (j != -1)
// 				legs[j] = legs[--n_leg];
// 			j = find_adj_idx(legs, n_leg, e2);
// 			if (j != -1)
// 				legs[j] = legs[--n_leg];
// 			resolve = 1;
// 			break;
// 		}
// 		ret += resolve;
// 	} while (resolve);
// 	return ret;
// }

static gint_t get_array_legs(struct asm_graph_t *g, gint_t *legs,
				khash_t(gint) *set_e, khash_t(gint) *set_leg)
{
	khiter_t k;
	gint_t e, ret = 0;
	for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
		if (kh_exist(set_leg, k)) {
			e = kh_key(set_leg, k);
			legs[ret++] = e;
			kh_del(gint, set_e, kh_get(gint, set_e, e));
			kh_del(gint, set_e, kh_get(gint, set_e, g->edges[e].rc_id));
		}
	}
	return ret;
}

gint_t check_simple_jungle_strict(struct asm_graph_t *g, khash_t(gint) *set_e,
				khash_t(gint) *set_leg, double uni_cov)
{
	gint_t *legs = alloca(kh_size(set_leg) * sizeof(gint_t));
	gint_t n_leg = get_array_legs(g, legs, set_e, set_leg);
	gint_t ret, resolve, i;
	// if (kh_size(set_e) == 0) {
	// 	ret = check_n_m_node_strict(g, legs, n_leg, uni_cov);
	// 	return ret;
	// }
	double uni_cov_local = callibrate_uni_cov(g, legs, n_leg, uni_cov);
	ret = 0;
	do {
		resolve = 0;
		gint_t e1, e2, et1;
		double fcov1, fcov2;
		struct cov_range_t rcov1, rcov2;
		uint32_t gap_size;
		for (i = 0; i < n_leg; ++i) {
			e1 = legs[i];
			fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
			rcov1 = convert_cov_range(fcov1);
			if (rcov1.hi != 1)
				continue;
			e2 = bc_find_pair_check_path_strict(g, set_e, e1, legs, n_leg);
			if (e2 < 0)
				continue;
			fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
			rcov2 = convert_cov_range(fcov2);
			if (rcov2.hi != 1)
				continue;
			et1 = bc_find_pair_check_path(g, set_e, e2, legs, n_leg);
			if (et1 != e1) {
				__VERBOSE("Not best pair (%ld, %ld) <-> %ld\n",
					e1, et1, e2);
				continue;
			}
			gap_size = get_distance(g, MIN_CONTIG_BARCODE, set_e,
				g->nodes[g->edges[e1].source].rc_id,
				g->edges[e2].source);
			if (gap_size)
				asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
					e2, g->edges[e2].rc_id, gap_size);
			else
				asm_join_edge(g, g->edges[e1].rc_id, e1,
							e2, g->edges[e2].rc_id);
			/* remove leg */
			gint_t j;
			j = find_adj_idx(legs, n_leg, e1);
			if (j != -1)
				legs[j] = legs[--n_leg];
			j = find_adj_idx(legs, n_leg, e2);
			if (j != -1)
				legs[j] = legs[--n_leg];
			resolve = 1;
			break;
		}
		ret += resolve;
	} while (resolve);
	return ret;
}

gint_t check_complex_jungle_strict(struct asm_graph_t *g, khash_t(gint) *set_e,
	khash_t(gint) *set_leg, khash_t(gint) *set_self, double uni_cov)
{
	gint_t *legs = alloca((kh_size(set_leg) + kh_size(set_self)) * sizeof(gint_t));
	gint_t n_leg = get_array_legs2(g, legs, set_e, set_leg, set_self);
	gint_t ret, resolve;
	double uni_cov_local = callibrate_uni_cov(g, legs, n_leg, uni_cov);
	do {
		resolve = 0;
		gint_t e1, e2, i;
		double fcov1, fcov2;
		struct cov_range_t rcov1, rcov2;
		uint32_t gap_size;
		for (i = 0; i < n_leg; ++i) {
			e1 = legs[i];
			if (kh_get(gint, set_leg, e1) == kh_end(set_leg))
				continue;
			fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
			rcov1 = convert_cov_range(fcov1);
			if (rcov1.hi != 1)
				continue;
			e2 = bc_find_pair_check_path_strict(g, set_e, e1, legs, n_leg);
			if (e2 < 0)
				continue;
			fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
			rcov2 = convert_cov_range(fcov2);
			if (rcov2.hi != 1)
				continue;
			gap_size = get_distance(g, MIN_CONTIG_BARCODE, set_e,
				g->nodes[g->edges[e1].source].rc_id,
				g->edges[e2].source);
			if (gap_size)
				asm_join_edge_with_pat(g, g->edges[e1].rc_id, e1,
					e2, g->edges[e2].rc_id, gap_size);
			else
				asm_join_edge(g, g->edges[e1].rc_id, e1,
							e2, g->edges[e2].rc_id);
			/* remove legs */
			gint_t j;
			if (kh_get(gint, set_self, e2) != kh_end(set_self)) {
				j = find_adj_idx(e2);
				if (j != -1)
					legs[j] = legs[--n_leg];
				j = find_adj_idx(g->edges[e2].rc_id);
				if (j != -1)
					legs[j] = legs[--n_leg];
			} else {
				j = find_adj_idx(e2);
				if (j != -1)
					legs[j] = legs[--n_leg];
				j = find_adj_idx(e1);
				if (j != -1)
					legs[j] = legs[--n_leg];
			}
		}
		ret += resolve;
	} while (resolve);
}

gint_t collapse_2_2_bridge_strict(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t e;
	gint_t cnt, cnt_local, ret;
	cnt = 0;
	do {
		cnt_local = 0;
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			ret = check_simple_loop(g, e, uni_cov);
			if (ret == 0) {
				ret = check_2_2_strict_bridge(g, e, uni_cov);
				cnt_local += ret;
			} else {
				cnt_local += ret;
			}
		}
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of resolved 2-2 strict bridges: %ld\n", cnt);
	return cnt;
}

gint_t collapse_2_2_small_bridge(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t cnt, cnt_local, ret, e;
	cnt = 0;
	do {
		cnt_local = 0;
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			ret = check_simple_loop(g, e, uni_cov);
			if (ret == 0) {
				ret = check_2_2_small_bridge(g, e, uni_cov);
				cnt_local += ret;
			} else {
				cnt_local += ret;
			}
		}
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of resolved 2-2 small bridges: %ld\n", cnt);
	return cnt;
}

gint_t collapse_simple_jungle_strict(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	khash_t(gint) *visited, *set_e, *set_v, *set_leg, *set_self;
	visited = kh_init(gint);
	set_e = kh_init(gint);
	set_v = kh_init(gint);
	set_leg = kh_init(gint);
	set_self = kh_init(gint);
	gint_t e, ret = 0;
	uint32_t n_leg, n_self;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		uint32_t len = get_edge_len(g->edges + e);
		if (kh_get(gint, visited, g->edges[e].target) != kh_end(visited) ||
			len < MIN_CONTIG_BARCODE || g->edges[e].source == -1)
			continue;
		find_region(g, e, MIN_CONTIG_BARCODE, MAX_EDGE_COUNT, set_v, set_e);
		if (kh_size(set_e) < MAX_EDGE_COUNT) {
			kh_merge_set(visited, set_v);
			detect_leg(g, MIN_CONTIG_BARCODE, set_v, set_e, set_leg, set_self);
			n_leg = kh_size(set_leg);
			n_self = kh_size(set_self);
			if (n_leg > 1 && n_self == 0)
				ret += check_simple_jungle_strict(g, set_e, set_leg, uni_cov);
		}
		kh_clear(gint, set_leg);
		kh_clear(gint, set_e);
		kh_clear(gint, set_v);
		kh_clear(gint, set_self);
	}
	__VERBOSE("Number of pair-edge join through simple jungle: %ld\n", ret);
	return ret;
}

gint_t collapse_complex_jungle_strict(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	khash_t(gint) *visited, *set_e, *set_v, *set_leg, *set_self;
	visited = kh_init(gint);
	set_e = kh_init(gint);
	set_v = kh_init(gint);
	set_leg = kh_init(gint);
	set_self = kh_init(gint);
	gint_t e, ret = 0;
	uint32_t n_leg, n_self;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		uint32_t len = get_edge_len(g->edges + e);
		if (kh_get(gint, visited, g->edges[e].target) != kh_end(visited) ||
			len < MIN_CONTIG_BARCODE || g->edges[e].source == -1)
			continue;
		find_region(g, e, MIN_CONTIG_BARCODE, MAX_EDGE_COUNT, set_v, set_e);
		if (kh_size(set_e) < MAX_EDGE_COUNT) {
			kh_merge_set(visited, set_v);
			detect_leg(g, MIN_CONTIG_BARCODE, set_v, set_e, set_leg, set_self);
			n_leg = kh_size(set_leg);
			n_self = kh_size(set_self);
			if (n_leg + n_self > 1)
				ret += check_complex_jungle_strict(g, set_e, set_leg, set_self, uni_cov);
		}
		kh_clear(gint, set_leg);
		kh_clear(gint, set_e);
		kh_clear(gint, set_v);
		kh_clear(gint, set_self);
	}
	__VERBOSE("Number of pair-edge join through complex jungle: %ld\n", ret);
	return ret;
}

gint_t collapse_n_m_bridge_strict(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t e, cnt_local, cnt = 0;
	do {
		cnt_local = 0;
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			cnt_local += check_n_m_bridge_strict(g, e, uni_cov);
		}
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of strictly joined path through bridge: %ld\n", cnt);
	return cnt;
}

gint_t collapse_n_m_node_strict(struct asm_graph_t *g0)
{
	double uni_cov = get_genome_coverage(g);
	gint_t u, cnt_local, cnt = 0;
	do {
		cnt_local = 0;
		for (u = 0; u < g->n_v; ++u)
			cnt_local += check_n_m_node_strict(g, u, uni_cov);
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of strictly joined path through node: %ld\n", cnt);
	return cnt;
}

void resolve_n_m_simple(struct asm_graph_t *g0, struct asm_graph_t *g1)
{
	gint_t cnt = 0, cnt_local;
	do {
		cnt_local = 0;
		cnt_local += collapse_2_2_bridge_strict(g0);
		cnt_local += collapse_2_2_small_bridge(g0);
		cnt_local += collapse_n_m_node_strict(g0);
		cnt_local += collapse_n_m_bridge_strict(g0);
		cnt += cnt_local;
	} while (cnt_local);
	test_asm_graph(g0);
	collapse_simple_jungle_strict(g0);
	collapse_complex_jungle_strict(g0);
	// collapse_simple_jungle(g0);
	// collapse_n_m_bridge(g0);
	asm_condense(g0, g1);
}

void resolve_n_m_jungle(struct asm_graph_t *g0, struct asm_graph_t *g1)
{
	// collapse_n_m_jungle(g0);
	asm_condense(g0, g1);
}

// void check_n_m_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
// {
// 	gint_t u, v, v_rc, u_rc, e1, e2, et1, e_rc, i;
// 	e_rc = g->edges[e].rc_id;
// 	v = g->edges[e].target;
// 	v_rc = g->nodes[v].rc_id;
// 	u = g->edges[e].source;
// 	u_rc = g->nodes[u].rc_id;
// 	if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
// 		(g->nodes[u_rc].deg < 2 && g->nodes[v].deg < 2))
// 		return;

// 	/* Callibrate uni_coverage */
// 	gint_t n_leg = g->nodes[u_rc].deg + g->nodes[v].deg + 1;
// 	gint_t *legs = alloca(n_leg * sizeof(gint_t));
// 	n_leg = 0;
// 	for (i = 0; i < g->nodes[u_rc].deg; ++i)
// 		legs[n_leg++] = g->nodes[u_rc].adj[i];
// 	for (i = 0; i < g->nodes[v].deg; ++i)
// 		legs[n_leg++] = g->nodes[v].adj[i];
// 	legs[n_leg++] = e;
// 	double uni_cov_local = callibrate_uni_cov(g, legs, n_leg, uni_cov);
// 	int ecov, cov1, cov2;
// 	uint64_t e_uni_count;
// 	ecov = __get_edge_cov_int(g, e, uni_cov_local);
// 	if (ecov == 0)
// 		ecov = 1;
// 	e_uni_count = g->edges[e].count / ecov;
// 	__VERBOSE("%ld-%ld bridge: %ld_%ld ~ %d\n", g->nodes[u_rc].deg,
// 		g->nodes[v].deg, e, g->edges[e].rc_id, ecov);

// 	int resolve;
// 	do {
// 		resolve = 0;
// 		for (i = 0; i < g->nodes[u_rc].deg; ++i) {
// 			e1 = g->nodes[u_rc].adj[i];
// 			cov1 = __get_edge_cov_int(g, e1, uni_cov_local);
// 			if (cov1 != 1)
// 				continue;
// 			e2 = bc_find_best_pair(g, e1, g->nodes[v].adj, g->nodes[v].deg);
// 			if (e2 == -1)
// 				continue;
// 			cov2 = __get_edge_cov_int(g, e2, uni_cov_local);
// 			if (cov2 == 1) {
// 				/* join edge e1.rc -> e -> e2 */
// 				__VERBOSE("cov1 = %d; cov2 = %d\n", cov1, cov2);
// 				et1 = bc_find_best_pair(g, e2, g->nodes[u_rc].adj, g->nodes[u_rc].deg);
// 				if (et1 != e1) {
// 					__VERBOSE("Not best pair: %ld <-> %ld && %ld <-> %ld\n",
// 						e1, e2, et1, e2);
// 					continue;
// 				}
// 				asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
// 					e2, g->edges[e2].rc_id, e_uni_count * cov1);
// 				ecov -= cov1;
// 			} else if (cov2 > 1) {
// 				g->edges = realloc(g->edges, (g->n_e + 2) * sizeof(struct asm_edge_t));
// 				g->n_e += 2;
// 				asm_clone_edge2(g, g->n_e - 2, e2);
// 				asm_clone_edge2(g, g->n_e - 1, g->edges[e2].rc_id);
// 				g->edges[g->n_e - 2].rc_id = g->n_e - 1;
// 				g->edges[g->n_e - 1].rc_id = g->n_e - 2;
// 				g->edges[g->n_e - 2].count = g->edges[g->n_e - 1].count = g->edges[e2].count / cov2;
// 				g->edges[e2].count = g->edges[g->edges[e2].rc_id].count = g->edges[e2].count / cov2 * (cov2 - 1);
// 				asm_add_node_adj(g, g->edges[g->n_e - 2].source, g->n_e - 2);
// 				asm_add_node_adj(g, g->edges[g->n_e - 1].source, g->n_e - 1);
// 				asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
// 					g->n_e - 2, g->n_e - 1, e_uni_count * cov1);
// 				ecov -= cov1;
// 			} else {
// 				/* cov2 == 0 */
// 				continue;
// 			}
// 			resolve = 1;
// 			break;
// 		}
// 	} while (resolve);

// 	/* special case when resolved all but 1 pair left */
// 	if (g->nodes[u_rc].deg == 1 && g->nodes[v].deg == 1) {
// 		e1 = g->nodes[u_rc].adj[0];
// 		e2 = g->nodes[v].adj[0];
// 		cov1 = __get_edge_cov_int(g, e1, uni_cov_local);
// 		cov2 = __get_edge_cov_int(g, e2, uni_cov_local);
// 		double ratio = get_barcode_ratio(g, e1, e2);
// 		/* FIXME: may need more complicated resolve here */
// 		__VERBOSE("Leftover path: %ld -> %ld -> %ld, ratio: %.6lf\n",
// 				e1, e, e2, ratio);
// 		if (__positive_ratio(ratio) ||
// 			(ratio < 0 && ecov >= cov1 && cov1 == cov2 && cov1 > 0)) {
// 			if (e1 == g->edges[e2].rc_id)
// 				asm_join_edge_loop(g, e1, e2, e, e_rc,
// 					e_uni_count * cov1);
// 			else
// 				asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
// 					e2, g->edges[e2].rc_id, e_uni_count * cov1);
// 			ecov -= cov1;
// 		}
// 		if (ecov && get_edge_len(g->edges + e) >= MIN_CONTIG_BARCODE) {
// 			/* Keep the bridge but make disconnected */
// 			gint_t v1, v2;
// 			v1 = asm_create_node(g);
// 			v2 = asm_create_node(g);
// 			asm_remove_node_adj(g, g->edges[e].source, e);
// 			asm_remove_node_adj(g, g->edges[e_rc].source, e_rc);
// 			asm_add_node_adj(g, v1, e);
// 			asm_add_node_adj(g, v2, e_rc);
// 			g->edges[e].target = g->nodes[v2].rc_id;
// 			g->edges[e_rc].target = g->nodes[v1].rc_id;
// 		} else {
// 			/* destroy the bridge */
// 			asm_remove_edge(g, e);
// 			asm_remove_edge(g, e_rc);
// 		}
// 	} else if (g->nodes[u_rc].deg == 0 && g->nodes[v].deg == 0) {
// 		if (ecov && get_edge_len(g->edges + e) >= MIN_CONTIG_BARCODE) {
// 			/* Keep the bridge but make disconnected */
// 			gint_t v1, v2;
// 			v1 = asm_create_node(g);
// 			v2 = asm_create_node(g);
// 			asm_remove_node_adj(g, g->edges[e].source, e);
// 			asm_remove_node_adj(g, g->edges[e_rc].source, e_rc);
// 			asm_add_node_adj(g, v1, e);
// 			asm_add_node_adj(g, v2, e_rc);
// 			g->edges[e].target = g->nodes[v2].rc_id;
// 			g->edges[e_rc].target = g->nodes[v1].rc_id;
// 		} else {
// 			/* destroy the bridge */
// 			asm_remove_edge(g, e);
// 			asm_remove_edge(g, e_rc);
// 		}
// 	} else {
// 		g->edges[e].count = g->edges[e_rc].count = e_uni_count * ecov;
// 	}
// }

// void collapse_n_m_complex(struct asm_graph_t *g, uint32_t min_contig_size,
// 	khash_t(gint) *set_e, khash_t(gint) *set_leg, khash_t(gint) *set_self,
// 	double uni_cov)
// {
// 	khiter_t k;
// 	int n_leg, n_self;
// 	// n_leg = kh_size(set_leg);
// 	// n_self = kh_size(set_self);
// 	n_leg = kh_size(set_leg) + kh_size(set_self);
// 	gint_t *legs, *self_loops;
// 	legs = alloca(n_leg * sizeof(gint_t));
// 	// self_loops = alloca(n_self * sizeof(gint_t));
// 	n_leg = 0;
// 	for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
// 		if (kh_exist(set_leg, k)) {
// 			gint_t e = kh_key(set_leg, k);
// 			legs[n_leg++] = e;
// 			kh_del(gint, set_e, kh_get(gint, set_e, e));
// 			kh_del(gint, set_e, kh_get(gint, set_e, g->edges[e].rc_id));
// 		}
// 	}
// 	// n_self = 0;
// 	// for (k = kh_begin(set_self); k != kh_end(set_self); ++k) {
// 	// 	if (kh_exist(set_self, k))
// 	// 		self_loops[n_self++] = kh_key(set_self, k);
// 	// }
// 	for (k = kh_begin(set_self); k != kh_end(set_self); ++k) {
// 		if (kh_exist(set_self, k))
// 			legs[n_leg++] = kh_key(set_self, k);
// 	}

// 	double uni_cov_local = callibrate_uni_cov(g, legs, n_leg, uni_cov);
// 	int resolve;
// 	do {
// 		resolve = 0;
// 		int i, cov1, cov2;
// 		gint_t e1, e2, et1;
// 		for (i = 0; i < n_leg; ++i) {
// 			e1 = legs[i];
// 			cov1 = __get_edge_cov_int(g, e1, uni_cov_local);
// 			if (cov1 != 1)
// 				continue;
// 			// e2 = bc_find_2set_strict(g, min_contig_size, set_e, e1,
// 			// 	legs, n_leg, self_loops, n_self);
// 			e2 = bc_find_pair_check_strict(g, min_contig_size, set_e, e1, legs, n_leg);
// 			if (e2 == -1)
// 				continue;
// 			gint_t gap_size;
// 			gap_size = get_distance(g, min_contig_size, set_e,
// 					g->nodes[g->edges[e1].source].rc_id,
// 					g->edges[e2].source);
// 			assert(gap_size != -1);
// 			// if (gap_size == -1) {
// 			// 	__VERBOSE("WARNING: not join edge %ld-%ld, no strict path\n",
// 			// 		e1, e2);
// 			// 	continue;
// 			// }
// 			cov2 = __get_edge_cov_int(g, e2, uni_cov_local);
// 			if (cov2 == 1) {
// 				// et1 = bc_find_2set_strict(g, min_contig_size, set_e, e2,
// 				// 	legs, n_leg, self_loops, n_self);
// 				et1 = bc_find_pair_check_strict(g, min_contig_size, set_e, e2, legs, n_leg);
// 				if (e1 != et1) {
// 					__VERBOSE("Not best pair %ld <-> %ld && %ld <-> %ld\n",
// 						e1, e2, et1, e2);
// 					continue;
// 				}
// 				asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1, e2,
// 					g->edges[e2].rc_id, gap_size);
// 				/* remove leg */
// 				gint_t j;
// 				j = find_adj_idx(legs, n_leg, e1);
// 				if (j != -1)
// 					legs[j] = legs[--n_leg];
// 				j = find_adj_idx(legs, n_leg, e2);
// 				if (j != -1)
// 					legs[j] = legs[--n_leg];
// 			} else if (cov2 > 1) {
// 				g->edges = realloc(g->edges, (g->n_e + 2) * sizeof(struct asm_edge_t));
// 				g->n_e += 2;
// 				asm_clone_edge2(g, g->n_e - 2, e2);
// 				asm_clone_edge2(g, g->n_e - 1, g->edges[e2].rc_id);
// 				g->edges[g->n_e - 2].rc_id = g->n_e - 1;
// 				g->edges[g->n_e - 1].rc_id = g->n_e - 2;
// 				g->edges[g->n_e - 2].count = g->edges[g->n_e - 1].count = g->edges[e2].count / cov2;
// 				g->edges[e2].count = g->edges[g->edges[e2].rc_id].count = g->edges[e2].count / cov2 * (cov2 - 1);
// 				asm_add_node_adj(g, g->edges[g->n_e - 2].source, g->n_e - 2);
// 				asm_add_node_adj(g, g->edges[g->n_e - 1].source, g->n_e - 1);
// 				asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
// 						g->n_e - 2, g->n_e - 1, gap_size);
// 				gint_t j;
// 				j = find_adj_idx(legs, n_leg, e1);
// 				if (j != -1)
// 					legs[j] = legs[--n_leg];
// 			} else {
// 				continue;
// 			}
// 			resolve = 1;
// 			break;
// 		}
// 	} while (resolve);
// }

// void collapse_simple_jungle(struct asm_graph_t *g)
// {
// 	double uni_cov = get_genome_coverage(g);
// 	__VERBOSE("Genome coverage: %.9lf\n", uni_cov);
// 	khash_t(gint) *visited, *set_e, *set_v, *set_leg, *set_self;
// 	visited = kh_init(gint);
// 	set_e = kh_init(gint);
// 	set_v = kh_init(gint);
// 	set_leg = kh_init(gint);
// 	set_self = kh_init(gint);
// 	gint_t e;
// 	uint32_t n_leg, n_self;
// 	for (e = 0; e < g->n_e; ++e) {
// 		if (g->edges[e].source == -1)
// 			continue;
// 		uint32_t len = get_edge_len(g->edges + e);
// 		if (kh_get(gint, visited, g->edges[e].target) != kh_end(visited) ||
// 			len < MIN_CONTIG_BARCODE || g->edges[e].source == -1)
// 			continue;
// 		find_region(g, e, MIN_CONTIG_BARCODE, MAX_EDGE_COUNT, set_v, set_e);
// 		if (kh_size(set_e) < MAX_EDGE_COUNT) {
// 			kh_merge_set(visited, set_v);
// 			detect_leg(g, MIN_CONTIG_BARCODE, set_v, set_e, set_leg, set_self);
// 			// __VERBOSE("Edge: %ld; len: %u; Edge count in region: %u; Number of leg: %u\n",
// 			// 	e, len, kh_size(set_e), kh_size(set_leg));
// 			n_leg = kh_size(set_leg);
// 			n_self = kh_size(set_self);
// 			// if (n_leg == 4 && n_self == 0) {
// 			// 	collapse_4_leg_complex(g, set_e, set_leg, uni_cov);
// 			// } else if (n_leg > 1 && n_self == 0) {
// 			// 	collapse_small_complex(g, set_e, set_leg, uni_cov);
// 			// }
// 		}
// 		kh_clear(gint, set_leg);
// 		kh_clear(gint, set_e);
// 		kh_clear(gint, set_v);
// 		kh_clear(gint, set_self);
// 	}
// }

// void collapse_n_m_bridge(struct asm_graph_t *g)
// {
// 	double uni_cov = get_genome_coverage(g);
// 	gint_t e;
// 	for (e = 0; e < g->n_e; ++e) {
// 		if (g->edges[e].source == -1)
// 			continue;
// 		check_n_m_bridge(g, e, uni_cov);
// 	}
// }

// void collapse_n_m_jungle(struct asm_graph_t *g)
// {
// 	double uni_cov = get_genome_coverage(g);
// 	khash_t(gint) *visited, *set_e, *set_v, *set_leg, *set_self;
// 	visited = kh_init(gint);
// 	set_e = kh_init(gint);
// 	set_v = kh_init(gint);
// 	set_leg = kh_init(gint);
// 	set_self = kh_init(gint);
// 	gint_t e;
// 	uint32_t n_leg, n_self;
// 	for (e = 0; e < g->n_e; ++e) {
// 		if (g->edges[e].source == -1)
// 			continue;
// 		uint32_t len = get_edge_len(g->edges + e);
// 		if (kh_get(gint, visited, g->edges[e].target) != kh_end(visited) ||
// 			len < MIN_CONTIG_BARCODE || g->edges[e].source == -1)
// 			continue;
// 		find_region(g, e, MIN_CONTIG_BARCODE, MAX_EDGE_COUNT, set_v, set_e);
// 		if (kh_size(set_e) < MAX_EDGE_COUNT) {
// 			kh_merge_set(visited, set_v);
// 			detect_leg(g, MIN_CONTIG_BARCODE, set_v, set_e, set_leg, set_self);
// 			// __VERBOSE("Edge: %ld; len: %u; Edge count in region: %u; Number of leg: %u\n",
// 			// 	e, len, kh_size(set_e), kh_size(set_leg));
// 			n_leg = kh_size(set_leg);
// 			n_self = kh_size(set_self);
// 			if (n_leg + n_self > 1) {
// 				__VERBOSE("Edge: %ld; Edge count in region: %u; Number of leg: %u; Number of pseudo-self-loop: %u\n",
// 					e, kh_size(set_e), n_leg, n_self);
// 				collapse_n_m_complex(g, MIN_CONTIG_BARCODE, set_e,
// 					set_leg, set_self, uni_cov);
// 			}
// 		}
// 		kh_clear(gint, set_leg);
// 		kh_clear(gint, set_e);
// 		kh_clear(gint, set_v);
// 		kh_clear(gint, set_self);
// 	}
// }

// void check_n_m_node(struct asm_graph_t *g, gint_t *legs,
// 						gint_t n_leg, double uni_cov)
// {
// 	double uni_cov_local = callibrate_uni_cov(g, legs, n_leg, uni_cov);
// 	gint_t i, e1, e2, et1;
// 	int cov1, cov2;
// 	int resolve;
// 	do {
// 		resolve = 0;
// 		for (i = 0; i < n_leg; ++i) {
// 			e1 = legs[i];
// 			cov1 = __get_edge_cov_int(g, e1, uni_cov_local);
// 			if (cov1 != 1)
// 				continue;
// 			e2 = bc_find_best_pair(g, e1, legs, n_leg);
// 			if (e2 == -1 || g->edges[e1].source != g->nodes[g->edges[e2].source].rc_id)
// 				continue;
// 			/* join edge e1.rc -> e -> e2 */
// 			cov2 = __get_edge_cov_int(g, e2, uni_cov_local);
// 			if (cov2 == 1) {
// 				et1 = bc_find_best_pair(g, e2, legs, n_leg);
// 				if (et1 != e1) {
// 					__VERBOSE("Not best pair %ld <-> %ld && %ld <-> %ld\n",
// 						e1, e2, et1, e2);
// 					continue;
// 				}
// 				asm_join_edge(g, g->edges[e1].rc_id, e1,
// 							e2, g->edges[e2].rc_id);
// 				/* remove leg */
// 				gint_t j;
// 				j = find_adj_idx(legs, n_leg, e1);
// 				if (j != -1)
// 					legs[j] = legs[--n_leg];
// 				j = find_adj_idx(legs, n_leg, e2);
// 				if (j != -1)
// 					legs[j] = legs[--n_leg];
// 			} else if (cov2 > 1) {
// 				g->edges = realloc(g->edges, (g->n_e + 2) * sizeof(struct asm_edge_t));
// 				g->n_e += 2;
// 				asm_clone_edge2(g, g->n_e - 2, e2);
// 				asm_clone_edge2(g, g->n_e - 1, g->edges[e2].rc_id);
// 				g->edges[g->n_e - 2].rc_id = g->n_e - 1;
// 				g->edges[g->n_e - 1].rc_id = g->n_e - 2;
// 				g->edges[g->n_e - 2].count = g->edges[g->n_e - 1].count = g->edges[e2].count / cov2;
// 				g->edges[e2].count = g->edges[g->edges[e2].rc_id].count = g->edges[e2].count / cov2 * (cov2 - 1);
// 				asm_add_node_adj(g, g->edges[g->n_e - 2].source, g->n_e - 2);
// 				asm_add_node_adj(g, g->edges[g->n_e - 1].source, g->n_e - 1);
// 				asm_join_edge(g, g->edges[e1].rc_id, e1,
// 							g->n_e - 2, g->n_e - 1);
// 				gint_t j;
// 				j = find_adj_idx(legs, n_leg, e1);
// 				if (j != -1)
// 					legs[j] = legs[--n_leg];
// 			} else {
// 				continue;
// 			}
// 			resolve = 1;
// 			break;
// 		}
// 	} while (resolve);
// 	/* Final case when more than 2 edges left, but cannot connect by barcode */
// 	if (n_leg >= 2) {
// 		g->nodes = realloc(g->nodes, (g->n_v + 2) * sizeof(struct asm_node_t));
// 		g->nodes[g->n_v].rc_id = g->n_v + 1;
// 		g->nodes[g->n_v + 1].rc_id = g->n_v;
// 		g->nodes[g->n_v].adj = g->nodes[g->n_v + 1].adj = NULL;
// 		g->nodes[g->n_v].deg = g->nodes[g->n_v + 1].deg = 0;
// 		gint_t src = g->edges[legs[0]].source, e, e_rc;
// 		for (i = 0; i < n_leg; ++i) {
// 			e = legs[i];
// 			if (g->edges[e].source == src) {
// 				e_rc = g->edges[e].rc_id;
// 				gint_t j = find_adj_idx(g->nodes[src].adj,
// 						g->nodes[src].deg, e);
// 				g->nodes[src].adj[j] = g->nodes[src].adj[--g->nodes[src].deg];
// 				g->edges[e].source = g->n_v;
// 				asm_add_node_adj(g, g->n_v, e);
// 				g->edges[e_rc].target = g->n_v + 1;
// 			}
// 		}
// 		g->n_v += 2;
// 		free(g->nodes[src].adj);
// 		g->nodes[src].adj = NULL;
// 		g->nodes[src].deg = 0;
// 	}
// }

// void collapse_small_complex(struct asm_graph_t *g, khash_t(gint) *set_e,
// 				khash_t(gint) *set_leg, double uni_cov)
// {
// 	khiter_t k;
// 	int n_leg, cov1, cov2;
// 	n_leg = kh_size(set_leg);
// 	__VERBOSE("n_leg = %d\n", n_leg);
// 	gint_t *legs = alloca(n_leg * sizeof(gint_t));
// 	n_leg = 0;
// 	for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
// 		if (kh_exist(set_leg, k)) {
// 			gint_t e = kh_key(set_leg, k);
// 			legs[n_leg++] = e;
// 			kh_del(gint, set_e, kh_get(gint, set_e, e));
// 			kh_del(gint, set_e, kh_get(gint, set_e, g->edges[e].rc_id));
// 		}
// 	}
// 	uint32_t gap_size = 0;
// 	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
// 		if (!kh_exist(set_e, k))
// 			continue;
// 		gint_t e = kh_key(set_e, k);
// 		gint_t len = get_edge_len(g->edges + e);
// 		int cov = __get_edge_cov_int(g, e, uni_cov);
// 		gap_size += cov * (len - g->ksize);
// 	}
// 	double uni_cov_local = callibrate_uni_cov(g, legs, n_leg, uni_cov);
// 	int resolve;
// 	do {
// 		resolve = 0;
// 		int i;
// 		gint_t e1, e2, et1;
// 		for (i = 0; i < n_leg; ++i) {
// 			e1 = legs[i];
// 			cov1 = __get_edge_cov_int(g, e1, uni_cov_local);
// 			if (cov1 != 1)
// 				continue;
// 			e2 = bc_find_pair_check_path(g, set_e, e1, legs, n_leg);
// 			if (e2 == -1)
// 				continue;
// 			// if (!check_path(g, set_e,
// 			// 		g->nodes[g->edges[e1].source].rc_id,
// 			// 		g->edges[e2].source)) {
// 			// 	__VERBOSE("WARNING: not join edge %ld - %ld, no path\n", e1, e2);
// 			// 	continue;
// 			// }
// 			cov2 = __get_edge_cov_int(g, e2, uni_cov_local);
// 			if (cov2 == 1) {
// 				et1 = bc_find_pair_check_path(g, set_e, e2, legs, n_leg);
// 				if (e1 != et1) {
// 					__VERBOSE("Not best pair %ld <-> %ld && %ld <-> %ld\n",
// 						e1, e2, et1, e2);
// 					continue;
// 				}
// 				asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1, e2,
// 					g->edges[e2].rc_id, gap_size);
// 				/* remove leg */
// 				gint_t j;
// 				j = find_adj_idx(legs, n_leg, e1);
// 				if (j != -1)
// 					legs[j] = legs[--n_leg];
// 				j = find_adj_idx(legs, n_leg, e2);
// 				if (j != -1)
// 					legs[j] = legs[--n_leg];
// 			} else if (cov2 > 1) {
// 				g->edges = realloc(g->edges, (g->n_e + 2) * sizeof(struct asm_edge_t));
// 				g->n_e += 2;
// 				asm_clone_edge2(g, g->n_e - 2, e2);
// 				asm_clone_edge2(g, g->n_e - 1, g->edges[e2].rc_id);
// 				g->edges[g->n_e - 2].rc_id = g->n_e - 1;
// 				g->edges[g->n_e - 1].rc_id = g->n_e - 2;
// 				g->edges[g->n_e - 2].count = g->edges[g->n_e - 1].count = g->edges[e2].count / cov2;
// 				g->edges[e2].count = g->edges[g->edges[e2].rc_id].count = g->edges[e2].count / cov2 * (cov2 - 1);
// 				asm_add_node_adj(g, g->edges[g->n_e - 2].source, g->n_e - 2);
// 				asm_add_node_adj(g, g->edges[g->n_e - 1].source, g->n_e - 1);
// 				asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
// 						g->n_e - 2, g->n_e - 1, gap_size);
// 				gint_t j;
// 				j = find_adj_idx(legs, n_leg, e1);
// 				if (j != -1)
// 					legs[j] = legs[--n_leg];
// 			} else {
// 				continue;
// 			}
// 			resolve = 1;
// 			break;
// 		}
// 	} while (resolve);
// 	// for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
// 	// 	if (!kh_exist(set_e, k))
// 	// 		continue;
// 	// 	gint_t e = kh_key(set_e, k);
// 	// 	/* Remove edges , e.i isolate the nodes */
// 	// 	__VERBOSE("Removing edge %ld\n", e);
// 	// 	asm_remove_edge(g, e);
// 	// }
// }

