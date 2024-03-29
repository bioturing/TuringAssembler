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
			log_info("%ld", kh_key(set_e, k));
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
		gint_t dist = get_distance(g, MIN_CONTIG_BARCODE, set_e,
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
		assert(g->edges[e].source != -1);
		gint_t dist = get_distance(g, MIN_CONTIG_BARCODE, set_e,
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
		log_info("Best = %ld(~%.6lf); second best = %ld(~%.6lf)",
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
		log_info("Best = %ld(~%.6lf); second best = %ld(~%.6lf)",
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
//	log_info("Global cov ~ %.6lf. Local cov ~ %.6lf", uni_cov, ret);
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
			if (__check_coverage(fcov[i], fcov[k], rcov[i], rcov[k]))
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
		if (!__check_coverage(fcov[0], fcov[2], rcov[0], rcov[2]) ||
			!__check_coverage(fcov[1], fcov[3], rcov[1], rcov[3]))
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
		if (!__check_coverage(fcov[0], fcov[3], rcov[0], rcov[3]) ||
			!__check_coverage(fcov[1], fcov[2], rcov[1], rcov[2]))
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
			if (__check_coverage(fcov[i], fcov[k], rcov[i], rcov[k]))
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
		if (!__check_coverage(fcov[0], fcov[2], rcov[0], rcov[2]) ||
			!__check_coverage(fcov[1], fcov[3], rcov[1], rcov[3]))
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
		if (!__check_coverage(fcov[0], fcov[3], rcov[0], rcov[3]) ||
			!__check_coverage(fcov[1], fcov[2], rcov[1], rcov[2]))
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
			if (rcov1.hi < 1 || rcov1.lo > 1)
				continue;
			// if (rcov1.hi != 1)
			// 	continue;
			e2 = bc_find_pair_strict(g, e1, g->nodes[v].adj, g->nodes[v].deg);
			if (e2 < 0)
				continue;
			fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
			rcov2 = convert_cov_range(fcov2);
			if (!__check_coverage(fcov1, fcov2, rcov1, rcov2))
				continue;
			et1 = bc_find_pair_strict(g, e2, g->nodes[u_rc].adj, g->nodes[u_rc].deg);
			if (et1 != e1) {
				log_info("Not best pair: (%ld, %ld) <-> %ld", e1, et1, e2);
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
		log_info("Leftover path: %ld(~%.6lf) -> %ld -> %ld(~%.6lf), ratio: %.6lf",
				e1, fcov1, e, e2, fcov2, ratio);
		if ((__positive_ratio(ratio) || ratio < 0) &&
			__check_coverage(fcov1, fcov2, rcov1, rcov2)) {
			if (e1 != g->edges[e2].rc_id) {
				asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
					e2, g->edges[e2].rc_id, g->edges[e].count);
				++ret;
			}
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

gint_t check_n_m_node_strict(struct asm_graph_t *g, gint_t u, double uni_cov)
{
	gint_t u_rc;
	u_rc = g->nodes[u].rc_id;
	if (u == u_rc)
		return 0;
	if (g->nodes[u].deg == 0 || g->nodes[u_rc].deg == 0)
		return 0;
	gint_t *legs_1, *legs_2;
	gint_t n_leg1, n_leg2, i, e, resolve, ret;
	legs_1 = alloca(g->nodes[u].deg * sizeof(gint_t));
	legs_2 = alloca(g->nodes[u_rc].deg * sizeof(gint_t));
	n_leg1 = n_leg2 = 0;
	for (i = 0; i < g->nodes[u].deg; ++i) {
		e = g->nodes[u].adj[i];
		if (get_edge_len(g->edges + e) < MIN_CONTIG_BARCODE)
			return 0;
		if (g->edges[e].source == g->edges[e].target ||
			g->edges[e].source == g->nodes[g->edges[e].target].rc_id)
			return 0;
		legs_1[n_leg1++] = e;
	}
	for (i = 0; i < g->nodes[u_rc].deg; ++i) {
		e = g->nodes[u_rc].adj[i];
		if (get_edge_len(g->edges + e) < MIN_CONTIG_BARCODE)
			return 0;
		if (g->edges[e].source == g->edges[e].target ||
			g->edges[e].source == g->nodes[g->edges[e].target].rc_id)
			return 0;
		legs_2[n_leg2++] = e;
	}
	ret = 0;
	gint_t e1, e2, et1;
	double fcov1, fcov2;
	struct cov_range_t rcov1, rcov2;
	do {
		resolve = 0;
		for (i = 0; i < g->nodes[u_rc].deg; ++i) {
			e1 = g->nodes[u_rc].adj[i];
			fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov;
			rcov1 = convert_cov_range(fcov1);
			if (rcov1.hi < 1 || rcov1.lo > 1)
				continue;
			e2 = bc_find_pair_strict(g, e1, g->nodes[u].adj, g->nodes[u].deg);
			if (e2 < 0)
				continue;
			fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov;
			rcov2 = convert_cov_range(fcov2);
			if (!__check_coverage(fcov1, fcov2, rcov1, rcov2))
				continue;
			et1 = bc_find_pair_strict(g, e2, g->nodes[u_rc].adj, g->nodes[u_rc].deg);
			if (et1 != e1) {
				log_info("Not best pair: (%ld, %ld) <-> %ld", e1, et1, e2);
				continue;
			}
			asm_join_edge(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
			resolve = 1;
			break;
		}
		ret += resolve;
	} while (resolve);
	if (g->nodes[u_rc].deg == 1 && g->nodes[u].deg == 1) {
		e1 = g->nodes[u].adj[0];
		e2 = g->nodes[u_rc].adj[0];
		fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov;
		fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov;
		rcov1 = convert_cov_range(fcov1);
		rcov2 = convert_cov_range(fcov2);
		double ratio = get_barcode_ratio(g, e1, e2);
		log_info("Leftover path: %ld(~%.6lf) ->%ld(~%.6lf), ratio: %.6lf",
			e1, fcov1, e2, fcov2, ratio);
		if ((__positive_ratio(ratio) || ratio < 0) &&
			__check_coverage(fcov1, fcov2, rcov1, rcov2)) {
			if (e1 != g->edges[e2].rc_id) {
				asm_join_edge(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
			} else {
				asm_remove_node_adj(g, g->edges[e1].source, e1);
				gint_t n = asm_create_node(g);
				g->edges[e1].source = n;
				g->nodes[n].adj = malloc(sizeof(gint_t));
				g->nodes[n].adj[0] = e1;
				g->nodes[n].deg = 1;
				g->edges[g->edges[e1].rc_id].target = g->nodes[n].rc_id;
			}
		} else {
			asm_remove_node_adj(g, g->edges[e1].source, e1);
			gint_t n = asm_create_node(g);
			g->edges[e1].source = n;
			g->nodes[n].adj = malloc(sizeof(gint_t));
			g->nodes[n].adj[0] = e1;
			g->nodes[n].deg = 1;
			g->edges[g->edges[e1].rc_id].target = g->nodes[n].rc_id;
		}
	}
	return ret;
}

static gint_t get_array_legs2(struct asm_graph_t *g, gint_t *legs,
		khash_t(gint) *set_e, khash_t(gint) *set_leg, khash_t(gint) *set_self)
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
	for (k = kh_begin(set_self); k != kh_end(set_self); ++k) {
		if (kh_exist(set_self, k)) {
			e = kh_key(set_self, k);
			legs[ret++] = e;
		}
	}
	return ret;
}

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
			if (rcov1.hi < 1 || rcov1.lo > 1)
				continue;
			e2 = bc_find_pair_check_path_strict(g, set_e, e1, legs, n_leg);
			if (e2 < 0)
				continue;
			fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
			rcov2 = convert_cov_range(fcov2);
			if (!__check_coverage(fcov1, fcov2, rcov1, rcov2))
				continue;
			et1 = bc_find_pair_check_path(g, set_e, e2, legs, n_leg);
			if (et1 != e1) {
				log_info("Not best pair (%ld, %ld) <-> %ld",
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
	ret = 0;
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
			if (rcov1.hi < 1 || rcov1.lo > 1)
				continue;
			e2 = bc_find_pair_check_path_strict(g, set_e, e1, legs, n_leg);
			if (e2 < 0)
				continue;
			fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
			rcov2 = convert_cov_range(fcov2);
			if (!__check_coverage(fcov1, fcov2, rcov1, rcov2))
				continue;
			gap_size = get_distance(g, MIN_CONTIG_BARCODE, set_e,
				g->nodes[g->edges[e1].source].rc_id,
				g->edges[e2].source);
			log_info("Join %ld v %ld (gap size = %u)", e1, e2, gap_size);
			if (gap_size)
				asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
					e2, g->edges[e2].rc_id, gap_size);
			else
				asm_join_edge(g, g->edges[e1].rc_id, e1,
							e2, g->edges[e2].rc_id);
			/* remove legs */
			gint_t j;
			j = find_adj_idx(legs, n_leg, e2);
			if (j != -1)
				legs[j] = legs[--n_leg];
			j = find_adj_idx(legs, n_leg, e1);
			if (j != -1)
				legs[j] = legs[--n_leg];
			if (kh_get(gint, set_self, e2) != kh_end(set_self))
				kh_del(gint, set_self, kh_get(gint, set_self, g->edges[e2].rc_id));
			resolve = 1;
			break;
		}
		ret += resolve;
	} while (resolve);
	return ret;
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
				++cnt_local;
			}
		}
		cnt += cnt_local;
	} while (cnt_local);
	log_info("Number of resolved 2-2 strict bridges: %ld", cnt);
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
				++cnt_local;
			}
		}
		cnt += cnt_local;
	} while (cnt_local);
	log_info("Number of resolved 2-2 small bridges: %ld", cnt);
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
	log_info("Number of pair-edge join through simple jungle: %ld", ret);
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
	log_info("Number of pair-edge join through complex jungle: %ld", ret);
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
	log_info("Number of strictly joined path through bridge: %ld", cnt);
	return cnt;
}

gint_t collapse_n_m_node_strict(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t u, cnt_local, cnt = 0;
	do {
		cnt_local = 0;
		for (u = 0; u < g->n_v; ++u)
			cnt_local += check_n_m_node_strict(g, u, uni_cov);
		cnt += cnt_local;
	} while (cnt_local);
	log_info("Number of strictly joined path through node: %ld", cnt);
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

