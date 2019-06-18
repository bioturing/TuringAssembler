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

#define __positive_ratio(r)		((r) + EPS >= 0.1)
#define __positive_ratio_unique(r)	((r) + EPS >= 0.08)
#define MAX_EDGE_COUNT			10000

#define __is_hang_edge(g, set_e, e) (kh_get(gint, set_e, (g)->edges[e].rc_id) == kh_end(set_e))

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
	gint_t *q = malloc(max_edge_count * 2 * sizeof(gint_t));
	uint32_t l, r, len;
	gint_t u, e, j;
	l = r = 0;
	cb_add(set_e, se);
	q[0] = g->nodes[g->edges[se].source].rc_id;
	cb_add(set_v, g->nodes[g->edges[se].source].rc_id);
	while (l <= r) {
		u = q[l++];
		for (j = 0; j < g->nodes[u].deg; ++j) {
			e = g->nodes[u].adj[j];
			cb_add(set_e, e);
			len = get_edge_len(g->edges + e);
			struct cov_range_t rcov = get_edge_cov_range(g, e, genome_cov);
			if (len < min_contig_len || (len < MIN_CONTIG_BARCODE && rcov.hi > 1)) {
				if (cb_add(set_v, g->edges[e].target) == 1) {
					if (r + 1 == max_edge_count * 2)
						goto clean_up;
					q[++r] = g->edges[e].target;
				}
			}
		}
		if (g->nodes[u].deg) {
			if (cb_add(set_v, g->nodes[u].rc_id) == 1) {
				if (r + 1 == max_edge_count * 2)
					goto clean_up;
				q[++r] = g->nodes[u].rc_id;
			}
		}
	}
clean_up:
	free(q);
}

void detect_leg(struct asm_graph_t *g, uint32_t min_contig_len,
	uint32_t max_molecule_len, khash_t(gint) *set_v, khash_t(gint) *set_e,
				khash_t(gint) *set_leg, khash_t(gint) *set_self)
{
	khiter_t k;
	gint_t e;
	int ret;
	uint32_t len;
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (!kh_exist(set_e, k))
			continue;
		e = kh_key(set_e, k);
		if (__is_hang_edge(g, set_e, e))
			kh_put(gint, set_leg, e, &ret);
	}
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (!kh_exist(set_e, k))
			continue;
		e = kh_key(set_e, k);
		if (kh_get(gint, set_leg, e) == kh_end(set_leg)) {
			len = get_edge_len(g->edges + e);
			if (len >= max_molecule_len) {
				kh_put(gint, set_leg, e, &ret);
				kh_put(gint, set_leg, g->edges[e].rc_id, &ret);
			} else if (len >= min_contig_len) {
				kh_put(gint, set_self, e, &ret);
			}
		}
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

/*************************** Helper check fucntions ***************************/

static inline int check_large_pair_superior(struct asm_graph_t *g, gint_t e1, gint_t e2, gint_t e2a)
{
	double r1, r2, r1u, r2u;
	r1 = get_barcode_ratio(g, e1, e2);
	r2 = get_barcode_ratio(g, e1, e2a);
	if (__ratio_greater(r1, r2) && __positive_ratio(r1))
		return 1;
	r1u = get_barcode_ratio_unique(g, e1, e2);
	r2u = get_barcode_ratio_unique(g, e1, e2a);
	if (__ratio_greater(r1u, r2u) && __positive_ratio_unique(r1u))
		return 1;
	return g->edges[e1].best_mate_contigs == e2 &&
		g->edges[e2].best_mate_contigs == e1;
}

static inline int check_medium_pair_superior(struct asm_graph_t *g, gint_t e1, gint_t e2, gint_t e2a)
{
	if (g->edges[e1].seq_len >= MIN_CONTIG_BARCODE &&
		g->edges[e2].seq_len >= MIN_CONTIG_BARCODE &&
		g->edges[e2a].seq_len >= MIN_CONTIG_BARCODE)
		return check_large_pair_superior(g, e1, e2, e2a);
	if (g->edges[e1].best_mate_contigs == e2 &&
		g->edges[e2].best_mate_contigs == e1)
		return 1;
	if (g->edges[e1].seq_len >= MIN_CONTIG_BARCODE &&
		g->edges[e2].seq_len >= MIN_CONTIG_BARCODE) {
		/* Do some naughty stuff */
	}
	return 0;
}

static inline int check_medium_pair_positive(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	if (g->edges[e1].seq_len >= MIN_CONTIG_BARCODE &&
		g->edges[e2].seq_len >= MIN_CONTIG_BARCODE) {
		double r = get_barcode_ratio(g, e1, e2);
		if (__positive_ratio(r))
			return 1;
	}
	return g->edges[e1].best_mate_contigs == e2 &&
		g->edges[e2].best_mate_contigs == e1;
}

static inline int check_large_pair_greater(struct asm_graph_t *g, gint_t e1, gint_t e2, gint_t e2a)
{
	double r1, r2, r1u, r2u;
	r1 = get_barcode_ratio(g, e1, e2);
	r2 = get_barcode_ratio(g, e1, e2a);
	if (r1 > r2 && __positive_ratio(r1))
		return 1;
	r1u = get_barcode_ratio_unique(g, e1, e2);
	r2u = get_barcode_ratio_unique(g, e1, e2a);
	if (r1u > r2u && __positive_ratio_unique(r1u))
		return 1;
	return g->edges[e1].best_mate_contigs == e2 &&
		g->edges[e2].best_mate_contigs == e1;
}

static inline int check_medium_pair_greater(struct asm_graph_t *g, gint_t e1, gint_t e2, gint_t e2a)
{
	if (g->edges[e1].seq_len >= MIN_CONTIG_BARCODE &&
		g->edges[e2].seq_len >= MIN_CONTIG_BARCODE &&
		g->edges[e2a].seq_len >= MIN_CONTIG_BARCODE)
		return check_large_pair_greater(g, e1, e2, e2a);
	if (g->edges[e1].best_mate_contigs == e2 &&
		g->edges[e2].best_mate_contigs == e1)
		return 1;
	return 0;
}

static inline int remove_array_element(gint_t *a, gint_t n, gint_t x)
{
	gint_t j;
	for (j = 0; j < n; ++j)
		if (a[j] == x)
			break;
	if (j < n) {
		a[j] = a[--n];
		return 1;
	} else {
		return 0;
	}
}

static inline void dfs_push_edge(gint_t **seq, gint_t *mseq, gint_t *lseq, gint_t u)
{
	if (*mseq == *lseq) {
		*mseq <<= 1;
		*seq = realloc(*seq, (*mseq) * sizeof(gint_t));
	}
	(*seq)[(*lseq)++] = u;
}

gint_t dfs_get_dist(struct asm_graph_t *g, khash_t(gint) *set_e,
	khash_t(gint) *vis, gint_t **seq, gint_t *mseq, gint_t *lseq,
							gint_t u, gint_t t)
{
	gint_t j, e, v, ret;
	int hash_ret;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		if (kh_get(gint, set_e, e) == kh_end(set_e))
			continue;
		v = g->edges[e].target;
		dfs_push_edge(seq, mseq, lseq, e);
		if (v == t)
			return 0;
		if (kh_get(gint, vis, v) == kh_end(vis)) {
			kh_put(gint, vis, v, &hash_ret);
			ret = dfs_get_dist(g, set_e, vis, seq, mseq, lseq, v, t);
			if (ret != -1)
				return ret + (g->edges[e].seq_len - g->ksize);
		}
		--(*lseq);
	}
	return -1;
}

gint_t dfs_get_dist_simple(struct asm_graph_t *g, khash_t(gint) *set_e,
					khash_t(gint) *vis, gint_t u, gint_t t)
{
	gint_t j, e, v, ret;
	int hash_ret;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		if (kh_get(gint, set_e, e) == kh_end(set_e))
			continue;
		v = g->edges[e].target;
		if (v == t)
			return 0;
		if (kh_get(gint, vis, v) == kh_end(vis)) {
			kh_put(gint, vis, v, &hash_ret);
			ret = dfs_get_dist_simple(g, set_e, vis, v, t);
			if (ret != -1)
				return ret + (g->edges[e].seq_len - g->ksize);
		}
	}
	return -1;
}

gint_t get_dist(struct asm_graph_t *g, khash_t(gint) *set_e, gint_t **seq,
				gint_t *mseq, gint_t *lseq, gint_t s, gint_t t)
{
	*lseq = 0;
	if (s == t) {
		return 0;
	}
	int hash_ret;
	khash_t(gint) *vis;
	vis = kh_init(gint);
	kh_put(gint, vis, s, &hash_ret);
	gint_t ret = dfs_get_dist(g, set_e, vis, seq, mseq, lseq, s, t);
	kh_destroy(gint, vis);
	return ret;
}

gint_t get_dist_simple(struct asm_graph_t *g, khash_t(gint) *set_e, gint_t s, gint_t t)
{
	if (s == t)
		return 0;
	int hash_ret;
	khash_t(gint) *vis;
	vis = kh_init(gint);
	kh_put(gint, vis, s, &hash_ret);
	gint_t ret = dfs_get_dist_simple(g, set_e, vis, s, t);
	kh_destroy(gint, vis);
	return ret;
}

void join_edge_path(struct asm_graph_t *g, gint_t e1, gint_t e2, gint_t *seq,
							gint_t lseq, double cov)
{
	gint_t i, e_rc1, e_rc2, add_count, add_len;
	add_len = 0;
	e_rc1 = g->edges[e1].rc_id;
	e_rc2 = g->edges[e2].rc_id;
	for (i = 0; i < lseq; ++i) {
		asm_append_seq(g->edges + e1, g->edges + seq[i], g->ksize);
		asm_append_seq(g->edges + e_rc2,
			g->edges + g->edges[seq[lseq - i - 1]].rc_id, g->ksize);
		add_len += g->edges[seq[i]].seq_len;
	}
	asm_join_edge(g, e1, g->edges[e1].rc_id, e2, g->edges[e2].rc_id);
	add_count = add_len * cov;
	g->edges[e1].count += add_count;
	g->edges[e_rc2].count += add_count;
}

static gint_t bc_find_pair(struct asm_graph_t *g, gint_t se, gint_t *adj, gint_t n)
{
	gint_t ret_e, sec_e, j, e;
	ret_e = -1; sec_e = -1;
	for (j = 0; j < n; ++j) {
		e = adj[j];
		if (ret_e == -1 || check_medium_pair_greater(g, se, e, ret_e)) {
			sec_e = ret_e;
			ret_e = e;
		} else if (sec_e == -1 || check_medium_pair_greater(g, se, e, sec_e)) {
			sec_e = e;
		}
	}
	if (ret_e == -1 || !check_medium_pair_positive(g, se, ret_e))
		return -1;
	if (sec_e != -1 && !check_medium_pair_superior(g, se, ret_e, sec_e)) {
		return -2;
	}
	return ret_e;
}

static gint_t bc_find_pair_check_path(struct asm_graph_t *g, khash_t(gint) *set_e,
					gint_t se, khash_t(gint) *set_leg)
{
	gint_t ret_e, sec_e, e;
	khiter_t k;
	ret_e = sec_e = -1;
	for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
		if (!kh_exist(set_leg, k))
			continue;
		e = kh_key(set_leg, k);
		if (get_dist_simple(g, set_e,
				g->nodes[g->edges[se].source].rc_id,
				g->edges[e].source) != -1) {
			if (ret_e == -1 || check_medium_pair_greater(g, se, e, ret_e)) {
				sec_e = ret_e;
				ret_e = e;
			} else if (sec_e == -1 || check_medium_pair_greater(g, se, e, sec_e)) {
				sec_e = e;
			}
		}
	}
	if (ret_e == -1 || !check_medium_pair_positive(g, se, ret_e)) {
		return -1;
	}
	if (sec_e != -1 && !check_medium_pair_superior(g, se, ret_e, sec_e)) {
		return -2;
	}
	return ret_e;
}

static gint_t bc_find_alter_check_path(struct asm_graph_t *g, khash_t(gint) *set_e,
		gint_t se, gint_t de, khash_t(gint) *set_candidate)
{
	gint_t e, sec_e, ret_e;
	khiter_t k;
	ret_e = sec_e = de;
	for (k = kh_begin(set_candidate); k != kh_end(set_candidate); ++k) {
		if (!kh_exist(set_candidate, k))
			continue;
		e = kh_key(set_candidate, k);
		if (get_dist_simple(g, set_e,
					g->nodes[g->edges[se].source].rc_id,
					g->edges[e].source) != -1) {
			if (ret_e < 0 || check_medium_pair_greater(g, se, e, ret_e)) {
				sec_e = ret_e;
				ret_e = e;
			} else if (sec_e < 0 || check_medium_pair_greater(g, se, e, sec_e)) {
				sec_e = e;
			}
		}
	}
	if (ret_e < 0 || !check_medium_pair_positive(g, se, ret_e)) {
		return -1;
	}
	if (sec_e >= 0 && !check_medium_pair_superior(g, se, ret_e, sec_e)) {
		return -2;
	}
	return ret_e;
}

static inline int is_seq_rc(uint32_t *seq1, uint32_t l1,
						uint32_t *seq2, uint32_t l2)
{
	if (l1 != l2)
		return 0;
	uint32_t c1, c2, i, k;
	for (i = 0; i < l1; ++i) {
		k = l1 - i - 1;
		c1 = __binseq_get(seq1, i);
		c2 = __binseq_get(seq2, k);
		if (c1 != (c2 ^ 3))
			return 0;
	}
	return 1;
}

int check_edge(struct asm_graph_t *g, gint_t e)
{
	gint_t e_rc = g->edges[e].rc_id;
	return is_seq_rc(g->edges[e].seq, g->edges[e].seq_len, g->edges[e_rc].seq, g->edges[e_rc].seq_len);
}

/*************************** Detail resolve ***********************************/

gint_t check_2_2_large_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t e_rc, v, v_rc, u, u_rc;
	int i, k, flag, cnt;
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
		if (g->edges[legs[i]].seq_len < MIN_CONTIG_BARCODE)
			return 0;
	}
	double uni_cov_local = callibrate_uni_cov(g, legs, 4, uni_cov);
	/* calculate the coverage range for all legs */
	struct cov_range_t rcov[4], ecov;
	double fcov[4], ratio[4];
	for (i = 0; i < 4; ++i) {
		fcov[i] = __get_edge_cov(g->edges + legs[i], g->ksize) / uni_cov_local;
		rcov[i] = convert_cov_range(fcov[i]);
	}
	// printf("%.3f\n", get_barcode_ratio(g, legs[0], legs[2]));
	// printf("%.3f\n", get_barcode_ratio(g, legs[0], legs[3]));
	// printf("%.3f\n", get_barcode_ratio(g, legs[1], legs[2]));
	// printf("%.3f\n", get_barcode_ratio(g, legs[1], legs[3]));
	/* calculate the ratio of each pair */
	if ((cnt = (check_large_pair_superior(g, legs[0], legs[2], legs[3]) +
			check_large_pair_superior(g, legs[1], legs[3], legs[2]))) >= 1) {
		if (cnt == 2) {
			if (!__check_coverage(fcov[0], fcov[2], rcov[0], rcov[2]) ||
				!__check_coverage(fcov[1], fcov[3], rcov[1], rcov[3])) {
				__VERBOSE("[Large Bridge] Incompatible coverage range %ld(%ld)\n", e, e_rc);
				return 0;
			}
			__VERBOSE("[Large Bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
				g->edges[legs[0]].rc_id, legs[0], e, e_rc,
				legs[2], g->edges[legs[2]].rc_id);
			__VERBOSE("[Large Bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
				g->edges[legs[1]].rc_id, legs[1], e, e_rc,
				legs[3], g->edges[legs[3]].rc_id);
			// printf("%.3f %.3f\n",
			// 	get_barcode_ratio(g, legs[0], legs[2]),
			// 	get_barcode_ratio_unique(g, legs[0], legs[2]));
			// printf("%.3f %.3f\n",
			// 	get_barcode_ratio(g, legs[1], legs[3]),
			// 	get_barcode_ratio_unique(g, legs[1], legs[3]));

			asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
				legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
			asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
				legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);

			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
			return 2;
		} else {
			ecov = get_edge_cov_range(g, e, uni_cov_local);
			if (ecov.hi < __max(rcov[0].lo + rcov[2].lo, rcov[1].lo + rcov[3].lo)) {
				__VERBOSE("[Large Bridge] Bridge coverage is too small to split %ld(%ld)\n", e, e_rc);
				return 0;
			}
			__VERBOSE("[Large Bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
				g->edges[legs[0]].rc_id, legs[0], e, e_rc,
				legs[2], g->edges[legs[2]].rc_id);
			__VERBOSE("[Large Bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
				g->edges[legs[1]].rc_id, legs[1], e, e_rc,
				legs[3], g->edges[legs[3]].rc_id);
			//printf("%.3f %.3f\n",
			//	get_barcode_ratio(g, legs[0], legs[2]),
			//	get_barcode_ratio_unique(g, legs[0], legs[2]));
			//printf("%.3f %.3f\n",
			//	get_barcode_ratio(g, legs[1], legs[3]),
			//	get_barcode_ratio_unique(g, legs[1], legs[3]));

			asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
				legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
			asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
				legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);

			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
			return 2;
		}
	} else if ((cnt = (check_large_pair_superior(g, legs[0], legs[3], legs[2]) +
			check_large_pair_superior(g, legs[1], legs[2], legs[3]))) >= 1) {
		if (cnt == 2) {
			if (!__check_coverage(fcov[0], fcov[3], rcov[0], rcov[3]) ||
				!__check_coverage(fcov[1], fcov[2], rcov[1], rcov[2])) {
				__VERBOSE("Incompatible coverage range %ld(%ld)\n", e, e_rc);
				return 0;
			}
			__VERBOSE("[Large bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
				g->edges[legs[0]].rc_id, legs[0], e, e_rc,
				legs[3], g->edges[legs[3]].rc_id);
			__VERBOSE("[Large bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
				g->edges[legs[1]].rc_id, legs[1], e, e_rc,
				legs[2], g->edges[legs[2]].rc_id);
			// printf("%.3f %.3f\n",
			// 	get_barcode_ratio(g, legs[0], legs[3]),
			// 	get_barcode_ratio_unique(g, legs[0], legs[3]));
			// printf("%.3f %.3f\n",
			// 	get_barcode_ratio(g, legs[1], legs[2]),
			// 	get_barcode_ratio_unique(g, legs[1], legs[2]));

			asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
				legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
			asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
				legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);

			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
			return 2;
		} else {
			ecov = get_edge_cov_range(g, e, uni_cov_local);
			if (ecov.hi < __max(rcov[0].lo + rcov[2].lo, rcov[1].lo + rcov[3].lo)) {
				__VERBOSE("[Large Bridge] Bridge coverage is too small to split %ld(%ld)\n", e, e_rc);
				return 0;
			}
			__VERBOSE("[Large bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
				g->edges[legs[0]].rc_id, legs[0], e, e_rc,
				legs[3], g->edges[legs[3]].rc_id);
			__VERBOSE("[Large bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
				g->edges[legs[1]].rc_id, legs[1], e, e_rc,
				legs[2], g->edges[legs[2]].rc_id);
			// printf("%.3f %.3f\n",
			// 	get_barcode_ratio(g, legs[0], legs[3]),
			// 	get_barcode_ratio_unique(g, legs[0], legs[3]));
			// printf("%.3f %.3f\n",
			// 	get_barcode_ratio(g, legs[1], legs[2]),
			// 	get_barcode_ratio_unique(g, legs[1], legs[2]));

			asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
				legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
			asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
				legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);

			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
			return 2;
		}
	}
	return 0;
}

gint_t check_2_2_medium_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t e_rc, v, v_rc, u, u_rc;
	int i, k, flag, cnt;
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
		if (g->edges[legs[i]].seq_len < MIN_CONTIG_READPAIR)
			return 0;
	}
	double uni_cov_local = callibrate_uni_cov(g, legs, 4, uni_cov);
	/* calculate the coverage range for all legs */
	struct cov_range_t rcov[4], ecov;
	double fcov[4];
	for (i = 0; i < 4; ++i) {
		fcov[i] = __get_edge_cov(g->edges + legs[i], g->ksize) / uni_cov_local;
		rcov[i] = convert_cov_range(fcov[i]);
	}
	if ((cnt = (check_medium_pair_superior(g, legs[0], legs[2], legs[3]) +
			check_medium_pair_superior(g, legs[1], legs[3], legs[2]))) >= 1) {
		if (cnt == 2) {
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
		} else {
			ecov = get_edge_cov_range(g, e, uni_cov_local);
			if (ecov.hi < __max(rcov[0].lo + rcov[2].lo, rcov[1].lo + rcov[3].lo)) {
				__VERBOSE("[Med Bridge] Bridge coverage is too small to split %ld(%ld)\n", e, e_rc);
				return 0;
			}
			asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
				legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
			asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
				legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);

			__VERBOSE("[Med Bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
				g->edges[legs[0]].rc_id, legs[0], e, e_rc,
				legs[2], g->edges[legs[2]].rc_id);
			__VERBOSE("[Med Bridge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
				g->edges[legs[1]].rc_id, legs[1], e, e_rc,
				legs[3], g->edges[legs[3]].rc_id);

			asm_remove_edge(g, e);
			asm_remove_edge(g, e_rc);
			return 2;
		}
	} else if ((cnt = (check_large_pair_superior(g, legs[0], legs[3], legs[2]) +
			check_large_pair_superior(g, legs[1], legs[2], legs[3]))) >= 1) {
		if (cnt == 2) {
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
		} else {
			ecov = get_edge_cov_range(g, e, uni_cov_local);
			if (ecov.hi < __max(rcov[0].lo + rcov[2].lo, rcov[1].lo + rcov[3].lo)) {
				__VERBOSE("[Med Bridge] Bridge coverage is too small to split %ld(%ld)\n", e, e_rc);
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
	}
	return 0;
}

gint_t check_n_m_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t e_rc, v, v_rc, u, u_rc, ei;
	e_rc = g->edges[e].rc_id;
	v = g->edges[e].target;
	v_rc = g->nodes[v].rc_id;
	u = g->edges[e].source;
	u_rc = g->nodes[u].rc_id;
	/* condition for n-m edge */
	if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
		(g->nodes[u_rc].deg < 2 && g->nodes[v].deg < 2))
		return 0;

	/* Callibrate uni_coverage */
	gint_t *legs, *legs1, *legs2;
	gint_t n_leg1, n_leg2, n_leg, i, resolve, ret, e1, e2, et1;
	legs = alloca((g->nodes[v].deg + g->nodes[u_rc].deg) * sizeof(gint_t));
	legs1 = alloca(g->nodes[v].deg * sizeof(gint_t));
	legs2 = alloca(g->nodes[u_rc].deg * sizeof(gint_t));
	n_leg = n_leg1 = n_leg2 = 0;
	for (i = 0; i < g->nodes[u_rc].deg; ++i) {
		ei = g->nodes[u_rc].adj[i];
		if (g->edges[ei].seq_len < MIN_CONTIG_READPAIR)
			continue;
		legs[n_leg++] = ei;
		legs1[n_leg1++] = ei;
	}
	for (i = 0; i < g->nodes[v].deg; ++i) {
		ei = g->nodes[v].adj[i];
		if (g->edges[ei].seq_len < MIN_CONTIG_READPAIR)
			continue;
		legs[n_leg++] = ei;
		legs2[n_leg2++] = ei;
	}

	double uni_cov_local, fcov1, fcov2, e_cov;
	struct cov_range_t rcov1, rcov2, e_rcov;
	uni_cov_local = callibrate_uni_cov(g, legs, n_leg, uni_cov);
	e_cov = __get_edge_cov(g->edges + e, g->ksize) / uni_cov_local;
	uint64_t add_count, sub_count;
	sub_count = 0;
	ret = 0;
	do {
		resolve = 0;
		for (i = 0; i < n_leg1; ++i) {
			e1 = legs1[i];
			fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
			rcov1 = convert_cov_range(fcov1);
			e2 = bc_find_pair(g, e1, legs2, n_leg2);
			if (e2 < 0)
				continue;
			fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
			rcov2 = convert_cov_range(fcov2);
			if (!__check_coverage(fcov1, fcov2, rcov1, rcov2))
				continue;
			et1 = bc_find_pair(g, e2, legs1, n_leg1);
			if (et1 != -1 && et1 != e1) {
				__VERBOSE("[n-m Edge] Not best pair: (%ld, %ld) <-> %ld\n", e1, et1, e2);
				continue;
			}
			add_count = g->edges[e].count * (fcov1 + fcov2) / 2.0 / e_cov;
			__VERBOSE("[n-m Edge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
				g->edges[e1].rc_id, e1, e, e_rc, e2, g->edges[e2].rc_id);
			gint_t etmp = g->edges[e1].rc_id;
			asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
				e2, g->edges[e2].rc_id, add_count);
			n_leg1 -= remove_array_element(legs1, n_leg1, e1);
			n_leg2 -= remove_array_element(legs2, n_leg2, e2);
			sub_count += add_count;
			resolve = 1;
			break;
		}
		ret += resolve;
	} while (resolve);
	g->edges[e].count -= sub_count;
	g->edges[e_rc].count -= sub_count;
	if (g->nodes[u_rc].deg == 1 && g->nodes[v].deg == 1) {
		e1 = g->nodes[u_rc].adj[0];
		e2 = g->nodes[v].adj[0];
		fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
		fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
		rcov1 = convert_cov_range(fcov1);
		rcov2 = convert_cov_range(fcov2);
		e_cov = __get_edge_cov(g->edges + e, g->ksize) / uni_cov_local;
		e_rcov = convert_cov_range(e_cov);
		if (g->edges[e1].seq_len >= MIN_CONTIG_READPAIR &&
			g->edges[e2].seq_len >= MIN_CONTIG_READPAIR) {
			if (check_medium_pair_positive(g, e1, e2) &&
				__check_coverage(fcov1, fcov2, rcov1, rcov2) &&
				__check_coverage(fcov1, e_cov, rcov1, e_rcov) &&
				__check_coverage(fcov2, e_cov, rcov2, e_rcov)) {
				__VERBOSE("n-m Edge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
					g->edges[e1].rc_id, e1, e, e_rc, e2, g->edges[e2].rc_id);
				asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
					e2, g->edges[e2].rc_id, g->edges[e].count);
				++ret;
			} else {
				asm_remove_edge(g, e);
				asm_remove_edge(g, e_rc);
			}
		} else {
			if (check_medium_pair_positive(g, e1, e2) &&
				__check_coverage(fcov1, fcov2, rcov1, rcov2) &&
				__check_coverage(fcov1, e_cov, rcov1, e_rcov) &&
				__check_coverage(fcov2, e_cov, rcov2, e_rcov)) {
				__VERBOSE("n-m Edge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
					g->edges[e1].rc_id, e1, e, e_rc, e2, g->edges[e2].rc_id);
				asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
					e2, g->edges[e2].rc_id, g->edges[e].count);
				++ret;
			} else {
				asm_remove_edge(g, e);
				asm_remove_edge(g, e_rc);
			}
		}
	} else {
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
	}
	return ret;
}

gint_t check_n_m_node(struct asm_graph_t *g, gint_t u, double uni_cov)
{
	gint_t u_rc;
	u_rc = g->nodes[u].rc_id;
	if (u == u_rc)
		return 0;
	if (g->nodes[u].deg == 0 || g->nodes[u_rc].deg == 0)
		return 0;
	gint_t *legs, *legs1, *legs2;
	gint_t n_leg1, n_leg2, n_leg, i, e, resolve, ret, e1, e2, et1;
	legs = alloca((g->nodes[u].deg + g->nodes[u_rc].deg) * sizeof(gint_t));
	legs1 = alloca(g->nodes[u].deg * sizeof(gint_t));
	legs2 = alloca(g->nodes[u_rc].deg * sizeof(gint_t));
	n_leg = n_leg1 = n_leg2 = 0;
	for (i = 0; i < g->nodes[u].deg; ++i) {
		e = g->nodes[u].adj[i];
		if (g->edges[e].seq_len < MIN_CONTIG_READPAIR)
			continue;
		legs[n_leg++] = e;
		legs1[n_leg1++] = e;
	}
	for (i = 0; i < g->nodes[u_rc].deg; ++i) {
		e = g->nodes[u_rc].adj[i];
		if (g->edges[e].seq_len < MIN_CONTIG_READPAIR)
			continue;
		legs[n_leg++] = e;
		legs2[n_leg2++] = e;
	}
	double uni_cov_local, fcov1, fcov2;
	struct cov_range_t rcov1, rcov2;
	uni_cov_local = callibrate_uni_cov(g, legs, n_leg, uni_cov);
	ret = 0;
	do {
		resolve = 0;
		for (i = 0; i < n_leg1; ++i) {
			e1 = legs1[i];
			fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
			rcov1 = convert_cov_range(fcov1);
			e2 = bc_find_pair(g, e1, legs2, n_leg2);
			if (e2 < 0)
				continue;
			fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
			rcov2 = convert_cov_range(fcov2);
			if (!__check_coverage(fcov1, fcov2, rcov1, rcov2))
				continue;
			et1 = bc_find_pair(g, e2, legs1, n_leg1);
			if (et1 != -1 && et1 != e1) {
				__VERBOSE("[n-m Node] Not best pair: (%ld, %ld) <-> %ld\n", e1, et1, e2);
				continue;
			}
			__VERBOSE("[n-m Node] Join %ld(%ld) <-> %ld(%ld)\n",
				g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
			asm_join_edge(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
			n_leg1 -= remove_array_element(legs1, n_leg1, e1);
			n_leg2 -= remove_array_element(legs2, n_leg2, e2);
			resolve = 1;
			break;
		}
		ret += resolve;
	} while (resolve);
	if (g->nodes[u_rc].deg == 1 && g->nodes[u].deg == 1) {
		e1 = g->nodes[u].adj[0];
		e2 = g->nodes[u_rc].adj[0];
		fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
		fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
		rcov1 = convert_cov_range(fcov1);
		rcov2 = convert_cov_range(fcov2);
		if (g->edges[e1].seq_len >= MIN_CONTIG_READPAIR &&
			g->edges[e2].seq_len >= MIN_CONTIG_READPAIR) {
			if (check_medium_pair_positive(g, e1, e2) && __check_coverage(fcov1, fcov2, rcov1, rcov2)) {
				__VERBOSE("[n-m Node] Join %ld(%ld) <-> %ld(%ld)\n",
					g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				asm_join_edge(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				++ret;
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
			if (__check_coverage(fcov1, fcov2, rcov1, rcov2)) {
				__VERBOSE("[n-m Node] Join %ld(%ld) <-> %ld(%ld)\n",
					g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				asm_join_edge(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				++ret;
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
	}
	return ret;
}

gint_t join_n_m_small_jungle(struct asm_graph_t *g, khash_t(gint) *set_e,
					khash_t(gint) *set_leg, double uni_cov)
{
	gint_t *legs = alloca(kh_size(set_leg) * sizeof(gint_t));
	gint_t n_leg, ret, resolve, i;
	n_leg = get_array_legs(g, legs, set_e, set_leg);
	double uni_cov_local = callibrate_uni_cov(g, legs, n_leg, uni_cov);
	ret = 0;
	gint_t mpath_seq = 0x100, lpath_seq;
	gint_t *path_seq = malloc(mpath_seq * sizeof(gint_t));
	khiter_t k;
	do {
		resolve = 0;
		gint_t e1, e2, et1;
		double fcov1, fcov2, gap_size;
		struct cov_range_t rcov1, rcov2;
		for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
			if (!kh_exist(set_leg, k))
				continue;
			e1 = kh_key(set_leg, k);
			fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
			rcov1 = convert_cov_range(fcov1);
			e2 = bc_find_pair_check_path(g, set_e, e1, set_leg);
			if (e2 < 0)
				continue;
			fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
			rcov2 = convert_cov_range(fcov2);
			if (!__check_coverage(fcov1, fcov2, rcov1, rcov2))
				continue;
			et1 = bc_find_pair_check_path(g, set_e, e2, set_leg);
			if (et1 != -1 && e1 != et1) {
				__VERBOSE("[Small Jungle] Not best pair (%ld, %ld) <-> %ld\n",
					e1, et1, e2);
				continue;
			}
			gap_size = get_dist(g, set_e, &path_seq, &mpath_seq,
				&lpath_seq, g->nodes[g->edges[e1].source].rc_id,
				g->edges[e2].source);
			__VERBOSE("[Small Jungle] Join %ld(%ld) <-> %ld(%ld)\n",
				g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
			if (gap_size < 1000)
				join_edge_path(g, g->edges[e1].rc_id, e2,
					path_seq, lpath_seq, uni_cov_local);
			else
				asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
					e2, g->edges[e2].rc_id, gap_size);
			/* remove legs */
			kh_del(gint, set_leg, kh_get(gint, set_leg, e1));
			kh_del(gint, set_leg, kh_get(gint, set_leg, e2));
			// n_leg -= remove_array_element(legs, n_leg, e1);
			// n_leg -= remove_array_element(legs, n_leg, e2);
			++resolve;
		}
		ret += resolve;
	} while (resolve);
	free(path_seq);
	return ret;
}

gint_t get_contig_array(gint_t *legs, khash_t(gint) *set_leg,
			khash_t(gint) *set_self, khash_t(gint) *set_e)
{
	khiter_t k;
	gint_t e, ret = 0;
	for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
		if (kh_exist(set_leg, k)) {
			e = kh_key(set_leg, k);
			legs[ret++] = e;
			kh_del(gint, set_e, kh_get(gint, set_e, e));
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

gint_t join_n_m_complex_jungle(struct asm_graph_t *g, khash_t(gint) *set_e,
		khash_t(gint) *set_leg, khash_t(gint) *set_self, double uni_cov)
{
	khiter_t k;
	// for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
	// 	if (!kh_exist(set_leg, k))
	// 		continue;
	// 	__VERBOSE("[Complex Jungle] legs = %ld\n", kh_key(set_leg, k));
	// }
	// for (k = kh_begin(set_self); k != kh_end(set_self); ++k) {
	// 	if (!kh_exist(set_self, k))
	// 		continue;
	// 	__VERBOSE("[Complex Jungle] self = %ld\n", kh_key(set_self, k));
	// }
	gint_t *contigs = alloca(kh_size(set_leg) * sizeof(gint_t));
	gint_t n_contig = get_contig_array(contigs, set_leg, set_self, set_e);
	gint_t ret, resolve;
	ret = 0;
	double uni_cov_local = callibrate_uni_cov(g, contigs, n_contig, uni_cov);
	gint_t mpath_seq = 0x100, lpath_seq;
	gint_t *path_seq = malloc(mpath_seq * sizeof(gint_t));
	do {
		resolve = 0;
		gint_t e1, e2, et1, ec;
		double fcov1, fcov2, gap_size;
		struct cov_range_t rcov1, rcov2;
		int flag, hash_ret;
		for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
			if (!kh_exist(set_leg, k))
				continue;
			e1 = kh_key(set_leg, k);
			flag = 1;
			fcov1 = __get_edge_cov(g->edges + e1, g->ksize) / uni_cov_local;
			rcov1 = convert_cov_range(fcov1);
			e2 = bc_find_pair_check_path(g, set_e, e1, set_leg);
			if (e2 >= 0) {
				et1 = bc_find_pair_check_path(g, set_e, e2, set_leg);
				if (et1 != -1 && e1 != et1) {
					__VERBOSE("[Complex Jungle] Not best pair (%ld, %ld) <-> %ld\n",
						e1, et1, e2);
					flag = 0;
				}
			}
			ec = bc_find_alter_check_path(g, set_e, e1, e2, set_self);
			if (ec >= 0) {
				fcov2 = __get_edge_cov(g->edges + ec, g->ksize) / uni_cov_local;
				rcov2 = convert_cov_range(fcov2);
				if (!__check_coverage(fcov1, fcov2, rcov1, rcov2))
					continue;
				et1 = bc_find_pair_check_path(g, set_e, ec, set_leg);
				if (et1 != -1 && e1 != et1) {
					__VERBOSE("[Complex Jungle] Not best pair (%ld, %ld) <-> %ld\n",
						e1, et1, e2);
				}
				gap_size = get_dist(g, set_e, &path_seq, &mpath_seq,
					&lpath_seq, g->nodes[g->edges[e1].source].rc_id,
					g->edges[ec].source);
				__VERBOSE("Complex Jungle] Join %ld(%ld) <-> %ld(%ld)\n",
						g->edges[e1].rc_id, e1, ec, g->edges[ec].rc_id);
				if (gap_size < 1000)
					join_edge_path(g, g->edges[e1].rc_id, e2,
						path_seq, lpath_seq, uni_cov_local);
				else
					asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
						ec, g->edges[ec].rc_id, gap_size);
				kh_del(gint, set_leg, kh_get(gint, set_leg, e1));
				kh_del(gint, set_self, kh_get(gint, set_self, ec));
				kh_del(gint, set_self, kh_get(gint, set_self, g->edges[ec].rc_id));
				kh_put(gint, set_leg, g->edges[ec].rc_id, &hash_ret);
				++resolve;
			} else if (e2 >= 0 && flag) {
				fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
				rcov2 = convert_cov_range(fcov2);
				if (!__check_coverage(fcov1, fcov2, rcov1, rcov2))
					continue;
				gap_size = get_dist(g, set_e, &path_seq, &mpath_seq,
					&lpath_seq, g->nodes[g->edges[e1].source].rc_id,
					g->edges[e2].source);
				__VERBOSE("[Small Jungle] Join %ld(%ld) <-> %ld(%ld)\n",
					g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				if (gap_size < 1000)
					join_edge_path(g, g->edges[e1].rc_id, e2,
						path_seq, lpath_seq, uni_cov_local);
				else
					asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
						e2, g->edges[e2].rc_id, gap_size);
				/* remove legs */
				kh_del(gint, set_leg, kh_get(gint, set_leg, e1));
				kh_del(gint, set_leg, kh_get(gint, set_leg, e2));
				// n_leg -= remove_array_element(legs, n_leg, e1);
				// n_leg -= remove_array_element(legs, n_leg, e2);
				++resolve;
			}
		}
	} while (resolve);
	free(path_seq);
	return 0;
}

/*************************** Iterate regions **********************************/

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

gint_t collapse_n_m_bridge(struct asm_graph_t *g)
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
				ret = check_n_m_bridge(g, e, uni_cov);
				cnt_local += ret;
			} else {
				++cnt_local;
			}
		}
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of joined n-m bridge: %ld\n", cnt);
	return cnt;
}

gint_t collapse_n_m_node(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t u, cnt_local, cnt;
	cnt = 0;
	do {
		cnt_local = 0;
		for (u = 0; u < g->n_v; ++u)
			cnt_local += check_n_m_node(g, u, uni_cov);
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of joined n-m node: %ld\n", cnt);
	return cnt;
}

/*************************** Process entry point ******************************/

void resolve_n_m_simple(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	gint_t cnt = 0, cnt_local;
	do {
		cnt_local = 0;
		// cnt_local += collapse_1_1_small_jungle(g0);
		cnt_local += collapse_2_2_large_bridge(g0);
		cnt_local += collapse_2_2_medium_bridge(g0);
		cnt_local += collapse_n_m_node(g0);
		cnt_local += collapse_n_m_bridge(g0);
		cnt += cnt_local;
	} while (cnt_local);
	asm_condense(g0, g);
}

void resolve_complex(struct asm_graph_t *g, struct asm_graph_t *gd)
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
		if (kh_get(gint, visited, e) != kh_end(visited) || len < MIN_CONTIG_BARCODE)
			continue;
		find_region(g, e, MIN_CONTIG_BARCODE, MAX_EDGE_COUNT, uni_cov, set_v, set_e);
		if (kh_size(set_e) < MAX_EDGE_COUNT) {
			kh_merge_set(visited, set_e);
			detect_leg(g, MIN_LONG_CONTIG, MAX_MOLECULE_LEN,
					set_v, set_e, set_leg, set_self);
			n_leg = kh_size(set_leg);
			n_self = kh_size(set_self);
			if (n_self == 0 && n_leg >= 2) {
				ret += join_n_m_small_jungle(g, set_e, set_leg, uni_cov);
			} else if (n_self + n_leg >= 2) {
				ret += join_n_m_complex_jungle(g, set_e, set_leg, set_self, uni_cov);
			}
		}
		kh_clear(gint, set_leg);
		kh_clear(gint, set_e);
		kh_clear(gint, set_v);
		kh_clear(gint, set_self);
	}
	__VERBOSE("Number of joined pair(s) through jungle: %ld\n", ret);
	asm_condense(g, gd);
}