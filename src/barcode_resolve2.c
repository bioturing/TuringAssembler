#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "assembly_graph.h"
#include "buffer_file_wrapper.h"
#include "io_utils.h"
#include "khash.h"
#include "kmer.h"
#include "sort_read.h"
#include "radix_sort.h"
#include "resolve.h"
#include "time_utils.h"
#include "utils.h"
#include "utils.h"
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

/*************************** Helper check functions ***************************/

static uint32_t count_shared_bc(struct barcode_hash_t *t1, struct barcode_hash_t *t2)
{
	uint32_t i, k, ret = 0;
	for (i = 0; i < t1->size; ++i) {
		if (t1->keys[i] == (uint64_t)-1)
			continue;
		k = barcode_hash_get(t2, t1->keys[i]);
		ret += (int)(k != BARCODE_HASH_END(t2));
	}
	return ret;
}

static inline struct barcode_hash_t *get_max_barcode_set(struct asm_graph_t *g, gint_t e, uint32_t len)
{
	if (len < CONTIG_USE_BARCODE)
		return NULL;
	if (len < CONTIG_LEVEL_0)
		return g->edges[e].barcodes;
	if (len < CONTIG_LEVEL_1)
		return g->edges[e].barcodes + 1;
	return g->edges[e].barcodes + 2;
}

static int check_barcode_positive(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	uint32_t len, shared;
	struct barcode_hash_t *h1, *h2;
	len = __min(g->edges[e1].seq_len, g->edges[e2].seq_len);
	h1 = get_max_barcode_set(g, e1, len);
	h2 = get_max_barcode_set(g, e2, len);
	if (h1 == NULL || h2 == NULL)
		return -1;
	shared = count_shared_bc(h1, h2);
	double ratio = shared * 1.0 / __min(h1->n_item, h2->n_item);
	return (ratio + EPS > MIN_BARCODE_RATIO);
}

static int check_barcode_superior(struct asm_graph_t *g, gint_t e1, gint_t e2, gint_t e2a)
{
	struct barcode_hash_t *h1, *h2, *h2a;
	h1 = get_max_barcode_set(g, e1, g->edges[e1].seq_len);
	uint32_t len2 = __min(g->edges[e2].seq_len, g->edges[e2a].seq_len);
	h2 = get_max_barcode_set(g, e2, len2);
	h2a = get_max_barcode_set(g, e2a, len2);
	if (h1 == NULL || h2 == NULL || h2a == NULL)
		return -1;
	uint32_t share_1_2, share_1_2a, share_1_2_2a, i, k2, k2a;
	share_1_2 = share_1_2a = share_1_2_2a = 0;
	for (i = 0; i < h1->size; ++i) {
		if (h1->keys[i] == (uint64_t)-1)
			continue;
		k2 = barcode_hash_get(h2, h1->keys[i]);
		k2a = barcode_hash_get(h2a, h1->keys[i]);
		share_1_2 += (k2 != BARCODE_HASH_END(h2));
		share_1_2a += (k2a != BARCODE_HASH_END(h2a));
		share_1_2_2a += (k2 != BARCODE_HASH_END(h2) &&
					k2a != BARCODE_HASH_END(h2a));
	}
	if (share_1_2 > share_1_2a * 2)
		return 1;
	// sub_share_1_2 = share_1_2 - share_1_2_2a;
	// sub_share_1_2a = share_1_2a - share_1_2_2a;
	// if (sub_share_1_2 > sub_share_1_2a * 2)
	// 	return 1;
	return 0;
}

static int check_barcode_greater(struct asm_graph_t *g, gint_t e1, gint_t e2, gint_t e2a)
{
	struct barcode_hash_t *h1, *h2, *h2a;
	h1 = get_max_barcode_set(g, e1, g->edges[e1].seq_len);
	uint32_t len2 = __min(g->edges[e2].seq_len, g->edges[e2a].seq_len);
	h2 = get_max_barcode_set(g, e2, len2);
	h2a = get_max_barcode_set(g, e2a, len2);
	if (h1 == NULL || h2 == NULL || h2a == NULL)
		return -1;
	uint32_t share_1_2, share_1_2a, share_1_2_2a, i, k2, k2a;
	share_1_2 = share_1_2a = share_1_2_2a = 0;
	for (i = 0; i < h1->size; ++i) {
		if (h1->keys[i] == (uint64_t)-1)
			continue;
		k2 = barcode_hash_get(h2, h1->keys[i]);
		k2a = barcode_hash_get(h2a, h1->keys[i]);
		share_1_2 += (k2 != BARCODE_HASH_END(h2));
		share_1_2a += (k2a != BARCODE_HASH_END(h2a));
		share_1_2_2a += (k2 != BARCODE_HASH_END(h2) &&
					k2a != BARCODE_HASH_END(h2a));
	}
	if (share_1_2 > share_1_2a)
		return 1;
	// sub_share_1_2 = share_1_2 - share_1_2_2a;
	// sub_share_1_2a = share_1_2a - share_1_2_2a;
	// if (sub_share_1_2 > sub_share_1_2a * 2)
	// 	return 1;
	return 0;
}

/* end check function */

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
			return g->edges[e].seq_len - g->ksize;
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
			return g->edges[e].seq_len - g->ksize;
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
		if (e == se || e == g->edges[se].rc_id)
			continue;
		if (check_barcode_positive(g, se, e)) {
			if (ret_e == -1 || check_barcode_greater(g, se, e, ret_e)) {
				sec_e = ret_e;
				ret_e = e;
			} else if (sec_e == -1 || check_barcode_greater(g, se, e, sec_e)) {
				sec_e = e;
			}
		}
	}
	if (ret_e == -1)
		return -1;
	if (sec_e != -1 && !check_barcode_superior(g, se, ret_e, sec_e)) {
		// printf("ret_e = %ld; sec_e = %ld\n", ret_e, sec_e);
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
		if (e == se || e == g->edges[se].rc_id)
			continue;
		if (check_barcode_positive(g, se, e)) {
			if (get_dist_simple(g, set_e,
					g->nodes[g->edges[se].source].rc_id,
					g->edges[e].source) != -1) {
				if (ret_e == -1 || check_barcode_greater(g, se, e, ret_e)) {
					sec_e = ret_e;
					ret_e = e;
				} else if (sec_e == -1 || check_barcode_greater(g, se, e, sec_e)) {
					sec_e = e;
				}
			}
		}
	}
	if (ret_e == -1) {
		return -1;
	}
	if (sec_e != -1 && !check_barcode_superior(g, se, ret_e, sec_e)) {
		return -2;
	}
	return ret_e;
}

static gint_t bc_find_alter_check_path(struct asm_graph_t *g, khash_t(gint) *set_e,
		gint_t se, gint_t de, khash_t(gint) *set_candidate)
{
	gint_t e, sec_e, ret_e;
	khiter_t k;
	ret_e = de;
	sec_e = -1;
	for (k = kh_begin(set_candidate); k != kh_end(set_candidate); ++k) {
		if (!kh_exist(set_candidate, k))
			continue;
		e = kh_key(set_candidate, k);
		if (e == se || e == g->edges[se].rc_id)
			continue;
		if (check_barcode_positive(g, se, e)) {
			if (get_dist_simple(g, set_e,
						g->nodes[g->edges[se].source].rc_id,
						g->edges[e].source) != -1) {
				if (ret_e < 0 || check_barcode_greater(g, se, e, ret_e)) {
					sec_e = ret_e;
					ret_e = e;
				} else if (sec_e < 0 || check_barcode_greater(g, se, e, sec_e)) {
					sec_e = e;
				}
			}
		}
	}
	if (sec_e >= 0 && !check_barcode_superior(g, se, ret_e, sec_e)) {
		return -2;
	}
	if (de >= 0 && ret_e == de) {
		return -3;
	}
	if (ret_e < 0) {
		return -1;
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

int check_2_2_high_strict_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
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

	/* calculate the coverage range for all legs */
	double uni_cov_local = callibrate_uni_cov(g, legs, 4, uni_cov);
	struct cov_range_t rcov[4], ecov;
	double fcov[4], ratio[4];
	for (i = 0; i < 4; ++i) {
		fcov[i] = __get_edge_cov(g->edges + legs[i], g->ksize) / uni_cov_local;
		rcov[i] = convert_cov_range(fcov[i]);
	}

	if (check_barcode_superior(g, legs[0], legs[2], legs[3]) == 1) {
		if (check_barcode_superior(g, legs[1], legs[3], legs[2]) == 0 ||
			check_barcode_superior(g, legs[2], legs[0], legs[1]) == 0 ||
			check_barcode_superior(g, legs[3], legs[1], legs[0]) == 0) {
			__VERBOSE("[High strict 1] Not enough condition for superior pairs %ld %ld\n", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[2]) == 0 ||
			check_barcode_positive(g, legs[1], legs[3]) == 0) {
			__VERBOSE("[High strict 1] Contradict condition for positive pairs %ld %ld\n", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[2], rcov[0], rcov[2]) ||
			!__check_coverage(fcov[1], fcov[3], rcov[1], rcov[3])) {
			__VERBOSE("[High strict 1] Incompatible coverage range %ld(%ld)\n", e, e_rc);
			return 0;
		}
		__VERBOSE("[High strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);
		__VERBOSE("[High strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;
	} else if (check_barcode_superior(g, legs[0], legs[3], legs[2]) == 1) {
		if (check_barcode_superior(g, legs[1], legs[2], legs[3]) == 0 ||
			check_barcode_superior(g, legs[2], legs[1], legs[0]) == 0 ||
			check_barcode_superior(g, legs[3], legs[0], legs[1]) == 0) {
			__VERBOSE("[High strict 2] Not enough condition for superior pairs %ld %ld\n", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[3]) == 0 ||
			check_barcode_positive(g, legs[1], legs[2]) == 0) {
			__VERBOSE("[High strict 2] Contradict for positive pairs %ld %ld\n", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[3], rcov[0], rcov[3]) ||
			!__check_coverage(fcov[1], fcov[2], rcov[1], rcov[2])) {
			__VERBOSE("[High strict 2] Incompatible coverage range %ld(%ld)\n", e, e_rc);
			return 0;
		}
		__VERBOSE("[High strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);
		__VERBOSE("[High strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;
	}
	return 0;
}

int check_2_2_med_strict_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
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

	/* calculate the coverage range for all legs */
	double uni_cov_local = callibrate_uni_cov(g, legs, 4, uni_cov);
	struct cov_range_t rcov[4], ecov;
	double fcov[4], ratio[4];
	for (i = 0; i < 4; ++i) {
		fcov[i] = __get_edge_cov(g->edges + legs[i], g->ksize) / uni_cov_local;
		rcov[i] = convert_cov_range(fcov[i]);
	}
	if (check_barcode_superior(g, legs[0], legs[2], legs[3]) == 1 ||
		check_barcode_superior(g, legs[1], legs[3], legs[2]) == 1 ||
		check_barcode_superior(g, legs[2], legs[0], legs[1]) == 1 ||
		check_barcode_superior(g, legs[3], legs[1], legs[0]) == 1) {
		if (check_barcode_greater(g, legs[0], legs[3], legs[2]) == 1 ||
			check_barcode_greater(g, legs[1], legs[2], legs[3]) == 1 ||
			check_barcode_greater(g, legs[2], legs[1], legs[0]) == 1 ||
			check_barcode_greater(g, legs[3], legs[0], legs[1]) == 1) {
			__VERBOSE("[Med strict 1] Contradict pair chosen %ld %ld\n", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[2]) == 0 ||
			check_barcode_positive(g, legs[1], legs[3]) == 0) {
			__VERBOSE("[Med strict 1] Contradict condition for positive pairs %ld %ld\n", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[2], rcov[0], rcov[2]) ||
			!__check_coverage(fcov[1], fcov[3], rcov[1], rcov[3])) {
			__VERBOSE("[Med strict 1] Incompatible coverage range %ld(%ld)\n", e, e_rc);
			return 0;
		}
		__VERBOSE("[Med strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);
		__VERBOSE("[Med strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;
	} else if (check_barcode_superior(g, legs[0], legs[3], legs[2]) == 1 ||
		check_barcode_superior(g, legs[1], legs[2], legs[3]) == 1 ||
		check_barcode_superior(g, legs[2], legs[1], legs[0]) == 1 ||
		check_barcode_superior(g, legs[3], legs[0], legs[1]) == 1) {
		if (check_barcode_greater(g, legs[0], legs[2], legs[3]) == 1 ||
			check_barcode_greater(g, legs[1], legs[3], legs[2]) == 1 ||
			check_barcode_greater(g, legs[2], legs[0], legs[1]) == 1 ||
			check_barcode_greater(g, legs[3], legs[1], legs[0]) == 1) {
			__VERBOSE("[Med strict 2] Constradict pair chosen %ld %ld\n", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[3]) == 0 ||
			check_barcode_positive(g, legs[1], legs[2]) == 0) {
			__VERBOSE("[Med strict 2] Contradict condition for positive pairs %ld %ld\n", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[3], rcov[0], rcov[3]) ||
			!__check_coverage(fcov[1], fcov[2], rcov[1], rcov[2])) {
			__VERBOSE("[Med strict 2] Incompatible coverage range %ld(%ld)\n", e, e_rc);
			return 0;
		}
		__VERBOSE("[Med strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);
		__VERBOSE("[Med strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;
	}
	return 0;
}

int check_2_2_low_strict_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
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

	/* calculate the coverage range for all legs */
	double uni_cov_local = callibrate_uni_cov(g, legs, 4, uni_cov);
	struct cov_range_t rcov[4], ecov;
	double fcov[4], ratio[4];
	for (i = 0; i < 4; ++i) {
		fcov[i] = __get_edge_cov(g->edges + legs[i], g->ksize) / uni_cov_local;
		rcov[i] = convert_cov_range(fcov[i]);
	}
	if (check_barcode_positive(g, legs[0], legs[2]) == 1 ||
		check_barcode_positive(g, legs[1], legs[3]) == 1) {
		if (check_barcode_positive(g, legs[0], legs[3]) == 1 ||
			check_barcode_positive(g, legs[1], legs[2]) == 1) {
			__VERBOSE("[Low strict 1] Contradict pair chosen %ld %ld\n", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[2]) == 0 ||
			check_barcode_positive(g, legs[1], legs[3]) == 0) {
			__VERBOSE("[Med strict 1] Contradict condition for positive pairs %ld %ld\n", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[2], rcov[0], rcov[2]) ||
			!__check_coverage(fcov[1], fcov[3], rcov[1], rcov[3])) {
			__VERBOSE("[Low strict 1] Incompatible coverage range %ld(%ld)\n", e, e_rc);
			return 0;
		}
		__VERBOSE("[Low strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);
		__VERBOSE("[Low strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;
	} else if (check_barcode_positive(g, legs[0], legs[3]) == 1 ||
		check_barcode_positive(g, legs[1], legs[2]) == 1) {
		if (check_barcode_positive(g, legs[0], legs[2]) == 1 ||
			check_barcode_positive(g, legs[1], legs[3]) == 1) {
			__VERBOSE("[Low strict 2] Contradict pair chosen %ld %ld\n", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[3]) == 0 ||
			check_barcode_positive(g, legs[1], legs[2]) == 0) {
			__VERBOSE("[Low strict 2] Contradict condition for positive pairs %ld %ld\n", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[3], rcov[0], rcov[3]) ||
			!__check_coverage(fcov[1], fcov[2], rcov[1], rcov[2])) {
			__VERBOSE("[Low strict 2] Incompatible coverage range %ld(%ld)\n", e, e_rc);
			return 0;
		}
		__VERBOSE("[Low strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);
		__VERBOSE("[Low strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
			g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);
		asm_join_edge3(g, g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id, g->edges[e].count / 2);
		asm_join_edge3(g, g->edges[legs[1]].rc_id, legs[1], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id, g->edges[e].count / 2);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 2;
	}
	return 0;
}

int check_n_m_bridge(struct asm_graph_t *g, gint_t e, double uni_cov)
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
	legs1 = alloca(g->nodes[u_rc].deg * sizeof(gint_t));
	legs2 = alloca(g->nodes[v].deg * sizeof(gint_t));
	n_leg = n_leg1 = n_leg2 = 0;
	for (i = 0; i < g->nodes[u_rc].deg; ++i) {
		ei = g->nodes[u_rc].adj[i];
		if (g->edges[ei].seq_len < CONTIG_USE_BARCODE)
			continue;
		legs[n_leg++] = ei;
		legs1[n_leg1++] = ei;
	}
	for (i = 0; i < g->nodes[v].deg; ++i) {
		ei = g->nodes[v].adj[i];
		if (g->edges[ei].seq_len < CONTIG_USE_BARCODE)
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
	if (sub_count <= g->edges[e].count) {
		g->edges[e].count -= sub_count;
		g->edges[e_rc].count -= sub_count;
	} else {
		g->edges[e].count = g->edges[e_rc].count = 0;
	}
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
			if (check_barcode_positive(g, e1, e2) &&
				__check_coverage(fcov1, fcov2, rcov1, rcov2) &&
				__check_coverage(fcov1, e_cov, rcov1, e_rcov) &&
				__check_coverage(fcov2, e_cov, rcov2, e_rcov)) {
			// if (check_barcode_positive(g, e1, e2) &&
			// 	__check_coverage(fcov1, fcov2, rcov1, rcov2)) {
				__VERBOSE("n-m Edge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)\n",
					g->edges[e1].rc_id, e1, e, e_rc, e2, g->edges[e2].rc_id);
				asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
					e2, g->edges[e2].rc_id, g->edges[e].count);
				++ret;
			}
		}
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
	} else if (g->nodes[u_rc].deg + g->nodes[v].deg == 1) {
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
	}
	return ret;
}

static inline void isolate_edge(struct asm_graph_t *g, gint_t e)
{
	asm_remove_node_adj(g, g->edges[e].source, e);
	gint_t n = asm_create_node(g);
	g->edges[e].source = n;
	g->nodes[n].adj = malloc(sizeof(gint_t));
	g->nodes[n].adj[0] = e;
	g->nodes[n].deg = 1;
	g->edges[g->edges[e].rc_id].target = g->nodes[n].rc_id;
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
	legs1 = alloca(g->nodes[u_rc].deg * sizeof(gint_t));
	legs2 = alloca(g->nodes[u].deg * sizeof(gint_t));
	n_leg = n_leg1 = n_leg2 = 0;
	for (i = 0; i < g->nodes[u_rc].deg; ++i) {
		e = g->nodes[u_rc].adj[i];
		if (g->edges[e].seq_len < CONTIG_USE_BARCODE)
			continue;
		legs[n_leg++] = e;
		legs1[n_leg1++] = e;
	}
	for (i = 0; i < g->nodes[u].deg; ++i) {
		e = g->nodes[u].adj[i];
		if (g->edges[e].seq_len < CONTIG_USE_BARCODE)
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
			if (check_barcode_positive(g, e1, e2) && __check_coverage(fcov1, fcov2, rcov1, rcov2)) {
				__VERBOSE("[n-m Node] Join %ld(%ld) <-> %ld(%ld)\n",
					g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				asm_join_edge(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				++ret;
			} else {
				isolate_edge(g, e1);
			}
		} else {
			if (__check_coverage(fcov1, fcov2, rcov1, rcov2)) {
				__VERBOSE("[n-m Node] Join %ld(%ld) <-> %ld(%ld)\n",
					g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				asm_join_edge(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				++ret;
			} else {
				isolate_edge(g, e1);
			}
		}
	}
	return ret;
}

gint_t join_1_1_small_jungle(struct asm_graph_t *g, khash_t(gint) *set_e,
					khash_t(gint) *set_leg, double uni_cov)
{
	gint_t *legs = alloca(kh_size(set_leg) * sizeof(gint_t));
	gint_t n_leg, e1, e2, e;
	khiter_t k;
	n_leg = get_array_legs(g, legs, set_e, set_leg);
	e1 = legs[0];
	e2 = legs[1];
	uint32_t gap_len = 0;
	double fcov;
	struct cov_range_t rcov;
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (!kh_exist(set_e, k))
			continue;
		e = kh_key(set_e, k);
		fcov = __get_edge_cov(g->edges + e, g->ksize) / uni_cov;
		rcov = convert_cov_range(fcov);
		gap_len += rcov.lo * (g->edges[e].seq_len - g->ksize);
	}
	asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id, gap_len / 2);
	return 1;
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
		gint_t e1, e2, et1, gap_size;
		double fcov1, fcov2;
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
			__VERBOSE("[Small Jungle] Join %ld(%ld) <-> %ld(%ld)\n",
				g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
			asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
				e2, g->edges[e2].rc_id, 50);
			// gap_size = get_dist(g, set_e, &path_seq, &mpath_seq,
			// 	&lpath_seq, g->nodes[g->edges[e1].source].rc_id,
			// 	g->edges[e2].source);
			// if (gap_size < 1000)
			// 	join_edge_path(g, g->edges[e1].rc_id, e2,
			// 		path_seq, lpath_seq, uni_cov_local);
			// else
			// 	asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
			// 		e2, g->edges[e2].rc_id, gap_size);
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
			kh_del(gint, set_e, kh_get(gint, set_e, e));
		}
	}
	return ret;
}

gint_t join_n_m_complex_jungle(struct asm_graph_t *g, khash_t(gint) *set_e,
		khash_t(gint) *set_leg, khash_t(gint) *set_self, double uni_cov)
{
	khiter_t k;
	for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
		if (!kh_exist(set_leg, k))
			continue;
		__VERBOSE("[Complex Jungle] legs = %ld\n", kh_key(set_leg, k));
	}
	for (k = kh_begin(set_self); k != kh_end(set_self); ++k) {
		if (!kh_exist(set_self, k))
			continue;
		__VERBOSE("[Complex Jungle] self = %ld\n", kh_key(set_self, k));
	}
	gint_t *contigs = alloca((kh_size(set_leg) + kh_size(set_self)) * sizeof(gint_t));
	gint_t n_contig = get_contig_array(contigs, set_leg, set_self, set_e);
	gint_t ret, resolve;
	ret = 0;
	double uni_cov_local = callibrate_uni_cov(g, contigs, n_contig, uni_cov);
	gint_t mpath_seq = 0x100, lpath_seq;
	gint_t *path_seq = malloc(mpath_seq * sizeof(gint_t));
	do {
		resolve = 0;
		gint_t e1, e2, et1, ec, gap_size;
		double fcov1, fcov2;
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
			__VERBOSE("e2 = %ld; ec = %ld\n", e2, ec);
			if (ec >= 0) {
				fcov2 = __get_edge_cov(g->edges + ec, g->ksize) / uni_cov_local;
				rcov2 = convert_cov_range(fcov2);
				if (!__check_coverage(fcov1, fcov2, rcov1, rcov2))
					continue;
				et1 = bc_find_pair_check_path(g, set_e, ec, set_leg);
				if (et1 != -1 && e1 != et1) {
					__VERBOSE("[Complex Jungle] Not best pair (%ld, %ld) <-> %ld\n",
						e1, et1, ec);
					continue;
				}
				// gap_size = get_dist(g, set_e, &path_seq, &mpath_seq,
				// 	&lpath_seq, g->nodes[g->edges[e1].source].rc_id,
				// 	g->edges[ec].source);
				// assert(gap_size != -1);
				__VERBOSE("[Complex Jungle] Join %ld(%ld) <-> %ld(%ld)\n",
						g->edges[e1].rc_id, e1, ec, g->edges[ec].rc_id);
				asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
					ec, g->edges[ec].rc_id, 50);
				// __VERBOSE("Join to self loop, distance = %ld\n", gap_size);
				// if (gap_size < 1000)
				// 	join_edge_path(g, g->edges[e1].rc_id, ec,
				// 		path_seq, lpath_seq, uni_cov_local);
				// else
				// 	asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
				// 		ec, g->edges[ec].rc_id, gap_size);
				kh_del(gint, set_leg, kh_get(gint, set_leg, e1));
				kh_del(gint, set_self, kh_get(gint, set_self, ec));
				kh_del(gint, set_self, kh_get(gint, set_self, g->edges[ec].rc_id));
				kh_put(gint, set_leg, g->edges[ec].rc_id, &hash_ret);
				++resolve;
			} else if (e2 >= 0 && flag && ec != -2) {
				fcov2 = __get_edge_cov(g->edges + e2, g->ksize) / uni_cov_local;
				rcov2 = convert_cov_range(fcov2);
				if (!__check_coverage(fcov1, fcov2, rcov1, rcov2))
					continue;
				// gap_size = get_dist(g, set_e, &path_seq, &mpath_seq,
				// 	&lpath_seq, g->nodes[g->edges[e1].source].rc_id,
				// 	g->edges[e2].source);
				// assert(gap_size != -1);
				__VERBOSE("[Complex Jungle] Join %ld(%ld) <-> %ld(%ld)\n",
					g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
					e2, g->edges[e2].rc_id, 500);
				// __VERBOSE("Join two legs, distance = %ld\n", gap_size);
				// if (gap_size < 1000)
				// 	join_edge_path(g, g->edges[e1].rc_id, e2,
				// 		path_seq, lpath_seq, uni_cov_local);
				// else
				// 	asm_join_edge_with_gap(g, g->edges[e1].rc_id, e1,
				// 		e2, g->edges[e2].rc_id, gap_size);
				/* remove legs */
				kh_del(gint, set_leg, kh_get(gint, set_leg, e1));
				kh_del(gint, set_leg, kh_get(gint, set_leg, e2));
				// n_leg -= remove_array_element(legs, n_leg, e1);
				// n_leg -= remove_array_element(legs, n_leg, e2);
				++resolve;
			}
		}
		ret += resolve;
	} while (resolve);
	free(path_seq);
	return ret;
}

/*************************** long loop ****************************************/

static inline int check_long_loop(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	gint_t u, v, u_rc, v_rc, j, e_return, e_return_rc, e_rc, e1, e2;
	int flag1, flag2, flag3;
	u = g->edges[e].source;
	v = g->edges[e].target;
	e_rc = g->edges[e].rc_id;
	u_rc = g->nodes[u].rc_id;
	v_rc = g->nodes[v].rc_id;
	if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
		g->nodes[u_rc].deg > 2 || g->nodes[v].deg > 2)
		return 0;
	e2 = e_return = -1;
	for (j = 0; j < g->nodes[v].deg; ++j) {
		if (g->edges[g->nodes[v].adj[j]].target == u)
			e_return = g->nodes[v].adj[j];
		else
			e2 = g->nodes[v].adj[j];
	}
	if (e_return == -1 || e2 == -1)
		return 0;
	e1 = e_return_rc = -1;
	for (j = 0; j < g->nodes[u_rc].deg; ++j) {
		if (g->edges[g->nodes[u_rc].adj[j]].target == v_rc)
			e_return_rc = g->nodes[u_rc].adj[j];
		else
			e1 = g->nodes[u_rc].adj[j];
	}
	if (e_return_rc != g->edges[e_return].rc_id) {
		__VERBOSE("Something happens\n");
		return 0;
	}
	double fcov_e, fcov_e_return;
	fcov_e = __get_edge_cov(g->edges + e, g->ksize) / uni_cov;
	fcov_e_return = __get_edge_cov(g->edges + e_return, g->ksize) / uni_cov;
	struct cov_range_t rcov_e, rcov_e_return;
	rcov_e = convert_cov_range(fcov_e);
	rcov_e_return = convert_cov_range(fcov_e_return);
	int rep = __min(rcov_e.hi - 1, rcov_e_return.hi);
	// rep = __min(rep, 2);
	if (rep <= 0) {
		if (g->edges[e_return].seq_len < MIN_NOTICE_LEN || fcov_e < fcov_e_return) {
			__VERBOSE("[Loop] Remove edge %ld(%ld)\n", e_return, e_return_rc);
			asm_remove_edge(g, e_return);
			asm_remove_edge(g, e_return_rc);
			return 1;
		} else {
			rep = 1;
		}
	}
	__VERBOSE("[Loop] Unroll %ld(%ld) <-> %ld(%ld) <-> %ld(%ld) rep = %d\n",
		e, e_rc, e_return, e_return_rc, e, e_rc, rep);
	asm_unroll_loop_forward(g, e, e_return, rep);
	asm_unroll_loop_forward(g, e_rc, e_return_rc, rep);
	asm_remove_edge(g, e_return);
	asm_remove_edge(g, e_return_rc);

	flag1 = flag2 = flag3 = 0;
	if (e1 != -1) {
		if (g->edges[e1].seq_len >= CONTIG_USE_BARCODE &&
			g->edges[e].seq_len >= CONTIG_USE_BARCODE)
			flag1 = check_barcode_positive(g, e1, e);
		else
			flag1 = 1;
	}

	if (e2 != -1) {
		if (g->edges[e2].seq_len >= CONTIG_USE_BARCODE &&
			g->edges[e].seq_len >= CONTIG_USE_BARCODE)
			flag2 = check_barcode_positive(g, e2, e_rc);
		else
			flag2 = 1;
	}

	if (e1 != -1 && e2 != -1 &&
		g->edges[e1].seq_len >= CONTIG_USE_BARCODE &&
		g->edges[e2].seq_len >= CONTIG_USE_BARCODE)
		flag3 = check_barcode_positive(g, e1, e2);

	if ((flag1 && flag2) || (flag3 && (flag1 || flag2 || g->edges[e].seq_len < MIN_NOTICE_LEN))) {
		__VERBOSE("[Loop] Join\n");
		asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
				e2, g->edges[e2].rc_id, g->edges[e].count);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 1;
	} else {
		__VERBOSE("[Loop] Break\n");
		if (!flag1)
			isolate_edge(g, e);
		if (!flag2)
			isolate_edge(g, e2);
		return 0;
	}
}

/*************************** Iterate regions **********************************/

gint_t resolve_2_2_bridge_high_strict(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t e, cnt, cnt_local, ret;
	cnt = 0;
	do {
		cnt_local = 0;
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			ret = check_long_loop(g, e, uni_cov);
			if (ret == 0) {
				ret = check_2_2_high_strict_bridge(g, e, uni_cov);
				cnt_local += ret;
			} else {
				++cnt_local;
			}
		}
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of joined 2-2 bridge (using barcode highly strict): %ld\n", cnt);
	return cnt;
}

gint_t resolve_2_2_bridge_med_strict(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t e, cnt, cnt_local, ret;
	cnt = 0;
	do {
		cnt_local = 0;
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			ret = check_long_loop(g, e, uni_cov);
			if (ret == 0) {
				ret = check_2_2_med_strict_bridge(g, e, uni_cov);
				cnt_local += ret;
			} else {
				++cnt_local;
			}
		}
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of joined 2-2 bridge (using barcode med strict): %ld\n", cnt);
	return cnt;
}

gint_t resolve_2_2_bridge_low_strict(struct asm_graph_t *g)
{
	double uni_cov = get_genome_coverage(g);
	gint_t e, cnt, cnt_local, ret;
	cnt = 0;
	do {
		cnt_local = 0;
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			ret = check_long_loop(g, e, uni_cov);
			if (ret == 0) {
				ret = check_2_2_low_strict_bridge(g, e, uni_cov);
				cnt_local += ret;
			} else {
				++cnt_local;
			}
		}
		cnt += cnt_local;
	} while (cnt_local);
	__VERBOSE("Number of joined 2-2 bridge (using barcode med strict): %ld\n", cnt);
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
			ret = check_long_loop(g, e, uni_cov);
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

void resolve_simple_complex(struct asm_graph_t *g)
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
			if (n_self == 0 && n_leg == 2)
				ret += join_1_1_small_jungle(g, set_e, set_leg, uni_cov);
		}
		kh_clear(gint, set_leg);
		kh_clear(gint, set_e);
		kh_clear(gint, set_v);
		kh_clear(gint, set_self);
	}
	__VERBOSE("Number of joined 1-1 pair(s) through jungle: %ld\n", ret);
}

void resolve_n_m_simple(struct asm_graph_t *g0, struct asm_graph_t *g)
{
	gint_t cnt = 0, cnt_local;
	do {
		cnt_local = 0;
		cnt_local += resolve_2_2_bridge_high_strict(g0);
		cnt_local += resolve_2_2_bridge_med_strict(g0);
		cnt_local += resolve_2_2_bridge_low_strict(g0);
		// cnt_local += collapse_1_1_small_jungle(g0);
		// cnt_local += collapse_2_2_large_bridge(g0);
		// cnt_local += collapse_2_2_medium_bridge(g0);
		// cnt_local += collapse_n_m_node(g0);
		// cnt_local += collapse_n_m_bridge(g0);
		cnt += cnt_local;
	} while (cnt_local);
	// resolve_simple_complex(g0);
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

static inline void dump_edge_fasta(char **seq, int *m_seq, struct asm_edge_t *e)
{
	int len, i, k;
	len = e->seq_len + (e->seq_len / 80) + (e->seq_len % 80 != 0);
	if (*m_seq < len + 1) {
		*m_seq = len + 1;
		*seq = realloc(*seq, *m_seq);
	}
	for (i = k = 0; i < (int)e->seq_len; ++i) {
		(*seq)[k++] = nt4_char[__binseq_get(e->seq, i)];
		if ((i + 1) % 80 == 0)
			(*seq)[k++] = '\n';
	}
	if (e->seq_len % 80 != 0)
		(*seq)[k++] = '\n';
	(*seq)[k] = '\0';
}

void write_fasta_seq(struct asm_graph_t *g, const char *fasta_path)
{
	FILE *fp = xfopen(fasta_path, "wb");
	gint_t e;
	char *seq = NULL;
	int m_seq = 0;
	for (e = 0; e < g->n_e; ++e) {
		dump_edge_fasta(&seq, &m_seq, g->edges + e);
		fprintf(fp, ">SEQ_%ld\n%s", e, seq);
	}
	free(seq);
	fclose(fp);
}

struct read_index_t {
	int64_t r1_offset;
	int64_t r2_offset;
	int64_t r1_len;
	int64_t r2_len;
};

#define read_index_get_key(p) ((p).r1_offset)
RS_IMPL(read_index, struct read_index_t, 64, 8, read_index_get_key);

KHASH_MAP_INIT_INT64(bcpos, struct read_index_t);

void construct_read_index(struct read_path_t *rpath, khash_t(bcpos) *h)
{
	FILE *fp = xfopen(rpath->idx_path, "rb");
	size_t byte_read;
	khint_t k;
	int ret;
	uint64_t barcode;
	char *buf = alloca(40);
	while ((byte_read = fread(buf, 1, 40, fp))) {
		if (byte_read != 40)
			__ERROR("Corrupted barcode in read index file");
		barcode = unpack_int64((uint8_t *)buf);
		k = kh_put(bcpos, h, barcode, &ret);
		if (ret != 1)
			__ERROR("Insert barcode failed");
		kh_value(h, k).r1_offset = unpack_int64((uint8_t *)buf + 8);
		kh_value(h, k).r2_offset = unpack_int64((uint8_t *)buf + 16);
		kh_value(h, k).r1_len = unpack_int64((uint8_t *)buf + 24);
		kh_value(h, k).r2_len = unpack_int64((uint8_t *)buf + 32);
	}
	fclose(fp);
}

void filter_read(struct read_path_t *ref, khash_t(bcpos) *dict,
		struct read_path_t *ans, uint64_t *shared, int n_shared)
{
	struct read_index_t *pos;
	pos = malloc(n_shared * sizeof(struct read_index_t));
	int i;
	khiter_t k;
	for (i = 0; i < n_shared; ++i) {
		k = kh_get(bcpos, dict, shared[i]);
		pos[i] = kh_value(dict, k);
	}
	rs_sort(read_index, pos, pos + n_shared);
	struct buffered_file_t fo1, fo2;
	bf_open(&fo1, ans->R1_path, "wb", SIZE_16MB);
	bf_open(&fo2, ans->R2_path, "wb", SIZE_16MB);
	FILE *fi1 = xfopen(ref->R1_path, "rb");
	FILE *fi2 = xfopen(ref->R2_path, "rb");
	int64_t m_buf, len;
	m_buf = 0x100;
	char *buf = malloc(m_buf);
	for (i = 0; i < n_shared; ++i) {
		len = __max(pos[i].r1_len, pos[i].r2_len);
		if (len > m_buf) {
			m_buf = len;
			buf = realloc(buf, m_buf);
		}
		fseek(fi1, pos[i].r1_offset, SEEK_SET);
		xfread(buf, 1, pos[i].r1_len, fi1);
		bf_write(&fo1, buf, pos[i].r1_len);
		fseek(fi2, pos[i].r2_offset, SEEK_SET);
		xfread(buf, 1, pos[i].r2_len, fi2);
		bf_write(&fo2, buf, pos[i].r2_len);
	}
	fclose(fi1);
	fclose(fi2);
	bf_close(&fo1);
	bf_close(&fo2);
	free(buf);
	free(pos);
}

void get_local_reads(struct read_path_t *reads, struct read_path_t *rpath,
			khash_t(bcpos) *dict, struct asm_graph_t *g,
			gint_t e1, gint_t e2, const char *prefix)
{
	char path[MAX_PATH];
	sprintf(path, "%s/R1.sub.fq", prefix);
	rpath->R1_path = strdup(path);
	sprintf(path, "%s/R2.sub.fq", prefix);
	rpath->R2_path = strdup(path);
	struct barcode_hash_t *h1, *h2;
	h1 = g->edges[e1].barcodes + 2;
	h2 = g->edges[e2].barcodes + 2;
	int n_shared, m_shared;
	n_shared = 0;
	m_shared = 0x100;
	uint64_t *shared = malloc(m_shared * sizeof(uint64_t));
	kmint_t i, k;
	for (i = 0; i < h1->size; ++i) {
		if (h1->keys[i] != (uint64_t)-1) {
			if (n_shared == m_shared) {
				m_shared <<= 1;
				shared = realloc(shared, m_shared * sizeof(uint64_t));
			}
			shared[n_shared++] = h1->keys[i];
		}
	}
	for (i = 0; i < h2->size; ++i) {
		if (h2->keys[i] == (uint64_t)-1)
			continue;
		k = barcode_hash_get(h1, h2->keys[i]);
		if (k == BARCODE_HASH_END(h1)) {
			if (n_shared == m_shared) {
				m_shared <<= 1;
				shared = realloc(shared, m_shared * sizeof(uint64_t));
			}
			shared[n_shared++] = h2->keys[i];
		}
	}
	filter_read(reads, dict, rpath, shared, n_shared);
	free(shared);
}

static void add_2_2_path(struct asm_graph_t *g, gint_t e1, gint_t e2, gint_t e, FILE *fp)
{
	uint32_t l1, l2, p1, p2, i, l, e_len, e1_len;
	gint_t e1_rc = g->edges[e1].rc_id;
	e_len = g->edges[e].seq_len;
	e1_len = g->edges[e1].seq_len;
	l1 = __min(g->edges[e1].seq_len, 200);
	l2 = __min(g->edges[e2].seq_len, 200);
	char *seq = malloc(l1 + l2 + e_len - g->ksize * 2 + 1);
	p1 = l1;
	p2 = p1 + g->edges[e].seq_len - 2 * g->ksize + 1;
	l = 0;
	for (i = 0; i < l1; ++i)
		seq[l++] = nt4_char[__binseq_get(g->edges[e1_rc].seq, e1_len - l1 + i)];
	for (i = g->ksize; i < e_len; ++i) {
		seq[l++] = nt4_char[__binseq_get(g->edges[e].seq, i)];
	}
	for (i = g->ksize; i < l2; ++i) {
		seq[l++] = nt4_char[__binseq_get(g->edges[e2].seq, i)];
	}
	struct pair_contig_t key = (struct pair_contig_t){e1, e2};
	int ret;
	khint_t k = kh_put(pair_contig_count, g->candidates, key, &ret);
	kh_value(g->candidates, k) = (struct contig_count_t){(int)0, (int)0};
	seq[l] = '\0';
	fprintf(fp, ">QRY_%ld_%ld_%u_%u_\n", e1, e2, p1, p2);
	fprintf(fp, "%s\n", seq);
	free(seq);
}

static void list_2_2_bridge(struct asm_graph_t *g, gint_t e, FILE *fp)
{
	gint_t e_rc, v, v_rc, u, u_rc;
	int i, k, flag, cnt;
	e_rc = g->edges[e].rc_id;
	v = g->edges[e].target;
	v_rc = g->nodes[v].rc_id;
	u = g->edges[e].source;
	u_rc = g->nodes[u].rc_id;
	/* condition for 2-2 edge */
	// if (g->nodes[u].deg != 1 || g->nodes[v_rc].deg != 1 ||
	// 	g->nodes[u_rc].deg != 2 || g->nodes[v].deg != 2)
	// 	return 0;
	for (i = 0; i < 2; ++i) {
		for (k = 0; k < 2; ++k) {
			add_2_2_path(g, g->nodes[u_rc].adj[i], g->nodes[v].adj[k], e, fp);
		}
	}
}

static int check_2_2_work(struct asm_graph_t *g, gint_t e)
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
	for (i = 0; i < 2; ++i) {
		for (k = 0; k < 2; ++k) {
			struct pair_contig_t key = (struct pair_contig_t){g->nodes[u_rc].adj[i], g->nodes[v].adj[k]};
			khint_t it = kh_get(pair_contig_count, g->candidates, key);
			if (it == kh_end(g->candidates))
				fprintf(stdout, "wowowowo\n");
			else
				fprintf(stdout, "%ld\t%ld\t%d\t%d\n", key.e1, key.e2,
					kh_value(g->candidates, it).n_read, kh_value(g->candidates, it).n_pair);
		}
	}
	return 0;
}

void resolve_local(struct opt_proc_t *opt, struct read_path_t *read_path,
				struct asm_graph_t *g, const char *work_dir)
{
	char path[MAX_PATH];
	sprintf(path, "%s/candidate.fasta", work_dir);
	FILE *fp = xfopen(path, "wb");
	gint_t e;
	g->candidates = kh_init(pair_contig_count);
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		list_2_2_bridge(g, e, fp);
	}
	fclose(fp);
	construct_aux_info(opt, g, read_path, path, ASM_BUILD_CANDIDATE);
	int resolved = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		resolved += check_2_2_work(g, e);
	}
}

static inline void get_first_edge_kmer(uint8_t *buf, int ksize, int word_size,
								uint32_t *seq)
{
	uint32_t c;
	int i;
	memset(buf, 0, word_size);
	for (i = 0; i < ksize; ++i) {
		c = __binseq_get(seq, i);
		km_shift_append(buf, ksize, word_size, c);
	}
}

static int find_path_label_helper(struct asm_graph_t *g,
		uint32_t *seq, uint32_t seq_len, gint_t **path, int *path_len,
						gint_t e, uint32_t e_begin)
{
	uint32_t i, k, cseq, cedge, ksize, word_size;
	uint8_t *kseq, *kedge;
	ksize = g->ksize + 1;
	word_size = (ksize + 3) >> 2;
	kseq = alloca(word_size);
	kedge = alloca(word_size);
	*path = malloc(sizeof(gint_t));
	*path_len = 1;
	(*path)[0] = e;
	for (i = 0, k = e_begin; i < seq_len; ++i) {
		cseq = __binseq_get(seq, i);
		km_shift_append(kseq, ksize, word_size, cseq);
		if (k == g->edges[e].seq_len) {
			if (i + 1 < ksize) {
				free(*path);
				*path = NULL;
				*path_len = 0;
				return 0;
			}
			gint_t v = g->edges[e].target, j, next_e;
			e = -1;
			for (j = 0; j < g->nodes[v].deg; ++j) {
				next_e = g->nodes[v].adj[j];
				get_first_edge_kmer(kedge, ksize, word_size,
					g->edges[next_e].seq);
				if (km_cmp(kseq, kedge, word_size) == 0) {
					e = next_e;
					break;
				}
			}
			if (e == -1) {
				free(*path);
				*path = NULL;
				*path_len = 0;
				return 0;
			}
			*path = realloc(*path, (*path_len + 1) * sizeof(gint_t));
			(*path)[*path_len++] = e;
			k = ksize;
		} else {
			cedge = __binseq_get(g->edges[e].seq, k);
			km_shift_append(kedge, ksize, word_size, cedge);
			++k;
		}
		if (i + 1 >= ksize) {
			if (km_cmp(kseq, kedge, word_size) != 0) {
				free(*path);
				*path = NULL;
				*path_len = 0;
				return 0;
			}
		}
	}
	return 1;
}

int find_path_label(struct asm_graph_t *g, uint32_t *seq, uint32_t seq_len, gint_t **path, int *path_len)
{
	/* get first (k + 1)-mer of seq */
	uint32_t i, c, ksize, word_size;
	int ret;
	uint8_t *kseq, *kedge;
	ksize = g->ksize + 1;
	word_size = (ksize + 3) >> 2;
	kseq = alloca(word_size);
	kedge = alloca(word_size);
	gint_t e;
	get_first_edge_kmer(kseq, ksize, word_size, seq);
	for (e = 0; e < g->n_e; ++e) {
		for (i = 0; i < g->edges[e].seq_len; ++i) {
			c = __binseq_get(g->edges[e].seq, i);
			km_shift_append(kedge, ksize, word_size, c);
			if (i + 1 >= ksize && km_cmp(kseq, kedge, word_size) == 0) {
				ret = find_path_label_helper(g, seq, seq_len, path, path_len, e, i + 1 - ksize);
				if (ret)
					return 1;
			}
		}
	}
	return 0;
}

void find_path_local(struct asm_graph_t *g0, struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	gint_t *epath = NULL;
	int epath_len = 0, i, ret;
	ret = find_path_label(g, g0->edges[e1].seq, g->edges[e1].seq_len, &epath, &epath_len);
	if (ret) {
		fprintf(stdout, "e1: ");
		for (i = 0; i < epath_len; ++i) {
			fprintf(stdout, "%ld, ", epath[i]);
			fprintf(stdout, "\n");
		}
	} else {
		fprintf(stdout, "e1: NULL\n");
	}
}

void test_local_assembly(struct opt_proc_t *opt, struct asm_graph_t *g,
							gint_t e1, gint_t e2)
{
	struct read_path_t read_sorted_path, local_read_path;
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
		sort_read(opt, &read_sorted_path);
	}
	khash_t(bcpos) *dict = kh_init(bcpos);
	construct_read_index(&read_sorted_path, dict);

	char work_dir[MAX_PATH];
	sprintf(work_dir, "%s/local_assembly_%ld_%ld", opt->out_dir, e1, e2);
	mkdir(work_dir, 0755);
	get_local_reads(&read_sorted_path, &local_read_path, dict, g, e1, e2, work_dir);
	struct asm_graph_t lg, lg1;
	build_local_assembly_graph(g->ksize, opt->n_threads, opt->mmem, 1,
		&(local_read_path.R1_path), &(local_read_path.R2_path), work_dir,
		&lg, g, e1, e2);
	build_0_1(&lg, &lg1);
	save_graph_info(work_dir, &lg1, "local");
	find_path_local(g, &lg1, e1, e2);
	// resolve_local(opt, &local_read_path, &lg1, work_dir);
}


void resolve_n_m_local(struct opt_proc_t *opt, struct read_path_t *rpath,
				struct asm_graph_t *g0, struct asm_graph_t *g1)
{
	// char path[MAX_PATH];
	// strcpy(path, opt->out_dir);
	// strcat(path, "/level4_working/");
	// mkdir(path, 0755);
	// strcat(path, "ref.fasta");
	// construct_fasta(g0, path);
	// construct_aux_info(opt, g0, rpath, path, 0);
}
