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

#define __is_hang_edge(g, set_v, e) (kh_get(gint, set_v, (g)->edges[e].target) == kh_end(set_v))

static inline int cb_add(khash_t(gint) *h, gint_t k)
{
	int ret;
	kh_put(gint, h, k, &ret);
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

static void kh_merge_set(khash_t(gint) *dst, khash_t(gint) *src)
{
	khiter_t k;
	int ret;
	for (k = kh_begin(src); k != kh_end(src); ++k) {
		if (!kh_exist(src, k))
			continue;
		kh_put(gint, dst, kh_key(src, k), &ret);
	}
}

void find_region(struct asm_graph_t *g, gint_t se, uint32_t min_contig_len,
				uint32_t max_edge_count, double genome_cov,
				khash_t(gint) *set_v, khash_t(gint) *set_e)
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
			struct cov_range_t rcov = get_edge_cov_range(g, ne, genome_cov);
			if (len < min_contig_len || rcov.hi > 1) {
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
			struct cov_range_t rcov = get_edge_cov_range(g, ne, genome_cov);
			if (len < min_contig_len || rcov.hi > 1) {
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
	return ret;
}

/****************************** Detail resolve ********************************/

gint_t check_2_2_large_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t e_rc, v, v_rc, u, u_rc;
	int i, k, flag;
	e_rc = g->edges[e].rc_id;
	v = g->edges[e].target;
	v_rc = g->nodes[v].rc_id;
	u = g->edges[e].source;
	u_rc = g->nodes[u].rc_id;
	/* condition for 2-2 edge */
	if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
		g->nodes[u_rc].deg != 2 || g->nodes[v].deg != 2)
		return 0;

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
	/* calculate the coverage range for all legs */
	struct cov_range_t rcov[4];
	double fcov[4];
	for (i = 0; i < 4; ++i) {
		fcov[i] = __get_edge_cov(g->edges + legs[i], g->ksize) / uni_cov_local;
		rcov[i] = convert_cov_range(fcov[i]);
	}
	/* calculate the ratio of each pair */
	double ratio[4];
	ratio[0] = get_barcode_ratio(g, legs[0], legs[2]);
	ratio[1] = get_barcode_ratio(g, legs[0], legs[3]);
	ratio[2] = get_barcode_ratio(g, legs[1], legs[2]);
	ratio[3] = get_barcode_ratio(g, legs[1], legs[3]);
	if (__ratio_greater(ratio[0], ratio[1]) && __ratio_greater(ratio[3], ratio[2])) {
		if ((!__positive_ratio(ratio[0]) || !__positive_ratio(ratio[3]))
			&& (__ratio_greater(ratio[0], ratio[3]) || __ratio_greater(ratio[3], ratio[0]))) {
			__VERBOSE("Negative ratio (0-2, 1-3) %ld(%ld) %.3f %.3f\n",
				e, e_rc, ratio[0], ratio[3]);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[2], rcov[0], rcov[2]) ||
			!__check_coverage(fcov[1], fcov[3], rcov[1], rcov[3])) {
			__VERBOSE("Incompatible coverage range %ld(%ld)\n", e, e_rc);
			return 0;
		}
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);

		__VERBOSE("[Large bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld) (ratio: %.3f)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, ratio[0]);
		__VERBOSE("[Large Bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld) (ratio: %.3f)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, ratio[3]);

		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;
	} else if (__ratio_greater(ratio[1], ratio[0]) && __ratio_greater(ratio[2], ratio[3])) {
		if ((!__positive_ratio(ratio[1]) || !__positive_ratio(ratio[2]))
			&& (__ratio_greater(ratio[1], ratio[2]) || __ratio_greater(ratio[2], ratio[1]))) {
			__VERBOSE("Negative ratio (0-3, 1-2) %ld(%ld) %.3f %.3f\n",
				e, e_rc, ratio[1], ratio[2]);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[3], rcov[0], rcov[3]) ||
			!__check_coverage(fcov[1], fcov[2], rcov[1], rcov[2])) {
			__VERBOSE("Incompatible coverage range %ld(%ld)\n", e, e_rc);
			return 0;
		}
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);

		__VERBOSE("[Large bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld) (ratio: %.3f)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, ratio[1]);
		__VERBOSE("[Large bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld) (ratio: %.3f)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, ratio[2]);

		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;
	}
	return 0;
}

gint_t check_2_2_medium_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t e_rc, v, v_rc, u, u_rc;
	int i, k, flag;
	e_rc = g->edges[e].rc_id;
	v = g->edges[e].target;
	v_rc = g->nodes[v].rc_id;
	u = g->edges[e].source;
	u_rc = g->nodes[u].rc_id;
	/* condition for 2-2 edge */
	if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
		g->nodes[u_rc].deg != 2 || g->nodes[v].deg != 2)
		return 0;

	gint_t *legs = alloca(4 * sizeof(gint_t));
	legs[0] = g->nodes[u_rc].adj[0];
	legs[1] = g->nodes[u_rc].adj[1];
	legs[2] = g->nodes[v].adj[0];
	legs[3] = g->nodes[v].adj[1];
	for (i = 0; i < 4; ++i) {
		if (get_edge_len(g->edges + legs[i]) < MIN_CONTIG_READPAIR)
			return 0;
	}
	double uni_cov_local = callibrate_uni_cov(g, legs, 4, uni_cov);
	/* calculate the coverage range for all legs */
	struct cov_range_t rcov[4];
	double fcov[4];
	for (i = 0; i < 4; ++i) {
		fcov[i] = __get_edge_cov(g->edges + legs[i], g->ksize) / uni_cov_local;
		rcov[i] = convert_cov_range(fcov[i]);
	}
	if (g->edges[legs[0]].best_mate_contigs == legs[2]) {
		if (g->edges[legs[1]].best_mate_contigs != legs[3]) {
			// print_test_pair_end(g, legs[0]);
			// print_test_pair_end(g, legs[1]);
			// print_test_pair_end(g, legs[2]);
			// print_test_pair_end(g, legs[3]);
			// print_test_barcode_edge(g, legs[1], legs[3]);
			__VERBOSE("[Med Bridge] Not so confident (0-2, 1-3) %ld(%ld)\n",
				e, e_rc);
			return 0;
		}

		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);

		__VERBOSE("[Med bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);
		__VERBOSE("[Med Bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);

		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;

	} else if (g->edges[legs[0]].best_mate_contigs == legs[3]) {
		if (g->edges[legs[1]].best_mate_contigs != legs[2]) {
			// print_test_pair_end(g, legs[0]);
			// print_test_pair_end(g, legs[1]);
			// print_test_pair_end(g, legs[2]);
			// print_test_pair_end(g, legs[3]);
			// print_test_barcode_edge(g, legs[1], legs[2]);
			__VERBOSE("[Med Bridge] Not so confident (0-3, 1-2) %ld(%ld)\n",
				e, e_rc);
			return 0;
		}
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);

		__VERBOSE("[Med bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);
		__VERBOSE("[Med bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);

		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;
	}
	return 0;
}

static inline int check_pair_contigs(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	// if (get_edge_len(g->edges + e1) >= MIN_CONTIG_BARCODE &&
	// 	get_edge_len(g->edges + e2) >= MIN_CONTIG_BARCODE) {
	// 	double ratio = get_barcode_ratio(g, e1, e2);
	// 	return __positive_ratio(ratio);
	// }
	if (get_edge_len(g->edges + e1) >= MIN_CONTIG_READPAIR &&
		get_edge_len(g->edges + e2) >= MIN_CONTIG_READPAIR) {
		return g->edges[e1].best_mate_contigs == e2 &&
			g->edges[e2].best_mate_contigs == e1;
	}
	return 0;
}

gint_t check_2_2_covbase_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t e_rc, v, v_rc, u, u_rc;
	int i, k, flag;
	e_rc = g->edges[e].rc_id;
	v = g->edges[e].target;
	v_rc = g->nodes[v].rc_id;
	u = g->edges[e].source;
	u_rc = g->nodes[u].rc_id;
	/* condition for 2-2 edge */
	if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
		g->nodes[u_rc].deg != 2 || g->nodes[v].deg != 2)
		return 0;

	gint_t *legs = alloca(4 * sizeof(gint_t));
	legs[0] = g->nodes[u_rc].adj[0];
	legs[1] = g->nodes[u_rc].adj[1];
	legs[2] = g->nodes[v].adj[0];
	legs[3] = g->nodes[v].adj[1];
	/* does not need to check length 
	 * for (i = 0; i < 4; ++i) {
	 * 	if (get_edge_len(g->edges + legs[i]) < MIN_CONTIG_READPAIR)
	 * 		return 0;
	 * }
	 */
	double uni_cov_local = callibrate_uni_cov(g, legs, 4, uni_cov);
	/* calculate the coverage range for all legs */
	struct cov_range_t rcov[4], ecov;
	double fcov[4];
	ecov = get_edge_cov_range(g, e, uni_cov_local);
	for (i = 0; i < 4; ++i) {
		fcov[i] = __get_edge_cov(g->edges + legs[i], g->ksize) / uni_cov_local;
		rcov[i] = convert_cov_range(fcov[i]);
	}
	if (ecov.hi < __max(rcov[0].lo + rcov[1].lo, rcov[2].lo + rcov[3].lo)) {
		__VERBOSE("[Covbase Bridge] Coverage of bridge is too low %ld(%ld)\n",
			e, e_rc);
		return 0;
	}
	if (check_pair_contigs(g, legs[0], legs[2]) || check_pair_contigs(g, legs[1], legs[3])) {
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);

		__VERBOSE("[Covbase bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);
		__VERBOSE("[Covbase Bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);

		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;
	} else if (check_pair_contigs(g, legs[0], legs[3]) || check_pair_contigs(g, legs[1], legs[2])) {
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);

		__VERBOSE("[Covbase bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);
		__VERBOSE("[Covbase bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);

		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;
	}
	return 0;
}

gint_t join_1_1_small_jungle(struct asm_graph_t *g, khash_t(gint) *set_e,
					khash_t(gint) *set_leg, double uni_cov)
{
	gint_t *legs = alloca(kh_size(set_leg) * sizeof(gint_t));
	gint_t n_leg = get_array_legs(g, legs, set_e, set_leg);
	gint_t e1 = legs[0], e2 = legs[1];
	double ratio = get_barcode_ratio(g, e1, e2);
	if (!__positive_ratio(ratio)) {
		__VERBOSE("rejected ratio = %.3f\n", ratio);
		return 0;
	}
	__VERBOSE("ok ratio = %.3f\n", ratio);
	uint32_t gap_size = 0;
	khiter_t k;
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (!kh_exist(set_e, k))
			continue;
		gint_t e = kh_key(set_e, k);
		if (kh_get(gint, set_leg, e) == kh_end(set_leg) &&
			kh_get(gint, set_leg, g->edges[e].rc_id) == kh_end(set_leg))
			continue;
		gap_size += g->edges[e].seq_len - g->ksize;
		asm_remove_edge(g, e);
	}
	gap_size = __max(gap_size / 2, (uint32_t)g->ksize);
	asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id, gap_size);
	return 1;
}

/*************************** Iterate regions **********************************/
gint_t collapse_1_1_small_jungle(struct asm_graph_t *g)
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
		find_region(g, e, MIN_CONTIG_BARCODE, MAX_EDGE_COUNT, uni_cov, set_v, set_e);
		if (kh_size(set_e) < MAX_EDGE_COUNT) {
			kh_merge_set(visited, set_v);
			detect_leg(g, MIN_LONG_CONTIG, set_v, set_e, set_leg, set_self);
			n_leg = kh_size(set_leg);
			n_self = kh_size(set_self);
			if (n_leg == 2 && n_self == 0)
				ret += join_1_1_small_jungle(g, set_e, set_leg, uni_cov);
		}
		kh_clear(gint, set_leg);
		kh_clear(gint, set_e);
		kh_clear(gint, set_v);
		kh_clear(gint, set_self);
	}
	__VERBOSE("Number of resolved 1-1 small jungle: %ld\n", ret);
	return ret;
}

gint_t collapse_2_2_large_bridge(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t e, cnt, cnt_local, ret;
	cnt = 0;
	do {
		cnt_local = 0;
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			ret = check_simple_loop(g, e, uni_cov);
			if (ret == 0) {
				ret = check_2_2_large_bridge(g, e, uni_cov);
				cnt_local += ret;
			} else {
				++cnt_local;
			}
		}
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of joined 2-2 large bridge (using barcode): %ld\n", cnt);
	return cnt;
}

gint_t collapse_2_2_medium_bridge(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t e, cnt, cnt_local, ret;
	cnt = 0;
	do {
		cnt_local = 0;
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			ret = check_simple_loop(g, e, uni_cov);
			if (ret == 0) {
				ret = check_2_2_medium_bridge(g, e, uni_cov);
				cnt_local += ret;
			} else {
				++cnt_local;
			}
		}
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of joined 2-2 medium bridge (using read pair): %ld\n", cnt);
	return cnt;
}

gint_t collapse_2_2_covbase_bridge(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t e, cnt, cnt_local, ret;
	cnt = 0;
	do {
		cnt_local = 0;
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			ret = check_simple_loop(g, e, uni_cov);
			if (ret == 0) {
				ret = check_2_2_covbase_bridge(g, e, uni_cov);
				cnt_local += ret;
			} else {
				++cnt_local;
			}
		}
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of joined 2-2 bridge (using both barcode and read pair): %ld\n", cnt);
	return cnt;
}


/*************************** Process entry point ******************************/

void resolve_n_m_simple(struct asm_graph_t *g0, struct asm_graph_t *g1)
{
	gint_t cnt = 0, cnt_local;
	do {
		cnt_local = 0;
		cnt_local += collapse_1_1_small_jungle(g0);
		cnt_local += collapse_2_2_large_bridge(g0);
		cnt_local += collapse_2_2_medium_bridge(g0);
		cnt_local += collapse_2_2_covbase_bridge(g0);
		cnt += cnt_local;
	} while (cnt_local);
	asm_condense(g0, g1);
}
