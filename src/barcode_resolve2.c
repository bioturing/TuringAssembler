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

gint_t join_1_1_small_jungle(struct asm_graph_t *g, khash_t(gint) *set_e,
					khash_t(gint) *set_leg, double uni_cov)
{
	gint_t *legs = alloca(kh_size(set_leg) * sizeof(gint_t));
	gint_t n_leg = get_array_legs(g, legs, set_e, set_leg);
	gint_t e1 = legs[0], e2 = legs[1];
	double ratio = get_barcode_ratio(g, e1, e2);
	uint32_t gap_size = 0;
	int ret = 0;
	khiter_t k;
	for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
		if (!kh_exist(set_leg, k))
			continue;
		gint_t e = kh_key(set_leg, k);
		gap_size += g->edges[e].seq_len - g->ksize;
		asm_remove_edge(g, e);
	}
	gap_size = __max(gap_size / 2, (uint32_t)g->ksize);
	if (__positive_ratio(ratio)) {
		asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id, gap_size);
		ret = 1;
	}
	return ret;
}

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

void resolve_n_m_simple(struct asm_graph_t *g0, struct asm_graph_t *g1)
{
	gint_t cnt = 0, cnt_local;
	do {
		cnt_local = 0;
		cnt_local += collapse_1_1_small_jungle(g0);
		cnt += cnt_local;
	} while (cnt_local);
}
