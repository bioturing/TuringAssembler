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
#include "radix_sort.h"
#include "resolve.h"
#include "time_utils.h"
#include "utils.h"
#include "utils.h"
#include "verbose.h"
#include "barcode_resolve2.h"
#include "process.h"
#include "barcode_builder.h"
#include "unit_test.h"
#include "resolve_big.h"
#include "build_hash_table.h"

#include "kmer_build.h"
#define MIN_EXCLUDE_BARCODE_CONTIG_LEN 6000
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

void destroy_read_path(struct read_path_t *reads)
{
	free(reads->R1_path);
	free(reads->R2_path);
	if (reads->idx_path != NULL)
		free(reads->idx_path);
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
		if (t1->keys[i] == (uint64_t) -1)
			continue;
		k = barcode_hash_get(t2, t1->keys[i]);
		ret += (int) (k != BARCODE_HASH_END(t2));
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
		if (h1->keys[i] == (uint64_t) -1)
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
		if (h1->keys[i] == (uint64_t) -1)
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

/*      legs0                legs2
 * ==========O            O=========
 *            \          /
 *             \        /
 *              O======O
 *             /        \
 *            /          \
 * ==========O            O=========
 *      lengs1               legs3
 */
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
			log_debug("[High strict 1] Not enough condition for superior pairs %ld %ld", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[2]) == 0 ||
			check_barcode_positive(g, legs[1], legs[3]) == 0) {
			log_debug("[High strict 1] Contradict condition for positive pairs %ld %ld", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[2], rcov[0], rcov[2]) ||
			!__check_coverage(fcov[1], fcov[3], rcov[1], rcov[3])) {
			log_debug("[High strict 1] Incompatible coverage range %ld(%ld)", e, e_rc);
			return 0;
		}
		log_debug("[High strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);
		log_debug("[High strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
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
			log_debug("[High strict 2] Not enough condition for superior pairs %ld %ld", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[3]) == 0 ||
			check_barcode_positive(g, legs[1], legs[2]) == 0) {
			log_debug("[High strict 2] Contradict for positive pairs %ld %ld", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[3], rcov[0], rcov[3]) ||
			!__check_coverage(fcov[1], fcov[2], rcov[1], rcov[2])) {
			log_debug("[High strict 2] Incompatible coverage range %ld(%ld)", e, e_rc);
			return 0;
		}
		log_debug("[High strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);
		log_debug("[High strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
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
			log_debug("[Med strict 1] Contradict pair chosen %ld %ld", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[2]) == 0 ||
			check_barcode_positive(g, legs[1], legs[3]) == 0) {
			log_debug("[Med strict 1] Contradict condition for positive pairs %ld %ld", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[2], rcov[0], rcov[2]) ||
			!__check_coverage(fcov[1], fcov[3], rcov[1], rcov[3])) {
			log_debug("[Med strict 1] Incompatible coverage range %ld(%ld)", e, e_rc);
			return 0;
		}
		log_debug("[Med strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);
		log_debug("[Med strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
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
			log_debug("[Med strict 2] Constradict pair chosen %ld %ld", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[3]) == 0 ||
			check_barcode_positive(g, legs[1], legs[2]) == 0) {
			log_debug("[Med strict 2] Contradict condition for positive pairs %ld %ld", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[3], rcov[0], rcov[3]) ||
			!__check_coverage(fcov[1], fcov[2], rcov[1], rcov[2])) {
			log_debug("[Med strict 2] Incompatible coverage range %ld(%ld)", e, e_rc);
			return 0;
		}
		log_debug("[Med strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);
		log_debug("[Med strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
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
			log_debug("[Low strict 1] Contradict pair chosen %ld %ld", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[2]) == 0 ||
			check_barcode_positive(g, legs[1], legs[3]) == 0) {
			log_debug("[Med strict 1] Contradict condition for positive pairs %ld %ld", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[2], rcov[0], rcov[2]) ||
			!__check_coverage(fcov[1], fcov[3], rcov[1], rcov[3])) {
			log_debug("[Low strict 1] Incompatible coverage range %ld(%ld)", e, e_rc);
			return 0;
		}
		log_debug("[Low strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[2], g->edges[legs[2]].rc_id);
		log_debug("[Low strict 1] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
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
			log_debug("[Low strict 2] Contradict pair chosen %ld %ld", e, e_rc);
			return 0;
		}
		if (check_barcode_positive(g, legs[0], legs[3]) == 0 ||
			check_barcode_positive(g, legs[1], legs[2]) == 0) {
			log_debug("[Low strict 2] Contradict condition for positive pairs %ld %ld", e, e_rc);
			return 0;
		}
		if (!__check_coverage(fcov[0], fcov[3], rcov[0], rcov[3]) ||
			!__check_coverage(fcov[1], fcov[2], rcov[1], rcov[2])) {
			log_debug("[Low strict 2] Incompatible coverage range %ld(%ld)", e, e_rc);
			return 0;
		}
		log_debug("[Low strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
			g->edges[legs[0]].rc_id, legs[0], e, e_rc,
			legs[3], g->edges[legs[3]].rc_id);
		log_debug("[Low strict 2] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
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
				log_debug("[n-m Edge] Not best pair: (%ld, %ld) <-> %ld", e1, et1, e2);
				continue;
			}
			add_count = g->edges[e].count * (fcov1 + fcov2) / 2.0 / e_cov;
			log_debug("[n-m Edge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
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
				log_debug("n-m Edge] Join %ld(%ld) <-> %ld(%ld) <-> %ld(%ld)",
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
				log_debug("[n-m Node] Not best pair: (%ld, %ld) <-> %ld", e1, et1, e2);
				continue;
			}
			log_debug("[n-m Node] Join %ld(%ld) <-> %ld(%ld)",
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
				log_debug("[n-m Node] Join %ld(%ld) <-> %ld(%ld)",
					g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				asm_join_edge(g, g->edges[e1].rc_id, e1, e2, g->edges[e2].rc_id);
				++ret;
			} else {
				isolate_edge(g, e1);
			}
		} else {
			if (__check_coverage(fcov1, fcov2, rcov1, rcov2)) {
				log_debug("[n-m Node] Join %ld(%ld) <-> %ld(%ld)",
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
				log_debug("[Small Jungle] Not best pair (%ld, %ld) <-> %ld",
					e1, et1, e2);
				continue;
			}
			log_debug("[Small Jungle] Join %ld(%ld) <-> %ld(%ld)",
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
		log_debug("Something happens");
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
			log_debug("[Loop] Remove edge %ld(%ld)", e_return, e_return_rc);
			asm_remove_edge(g, e_return);
			asm_remove_edge(g, e_return_rc);
			return 1;
		} else {
			rep = 1;
		}
	}
	log_debug("[Loop] Unroll %ld(%ld) <-> %ld(%ld) <-> %ld(%ld) rep = %d",
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
		log_debug("[Loop] Join");
		asm_join_edge3(g, g->edges[e1].rc_id, e1, e, e_rc,
				e2, g->edges[e2].rc_id, g->edges[e].count);
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		return 1;
	} else {
		log_debug("[Loop] Break");
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
	log_debug("Number of joined 2-2 bridge (using barcode highly strict): %ld", cnt);
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
	log_debug("Number of joined 2-2 bridge (using barcode med strict): %ld", cnt);
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
	log_debug("Number of joined 2-2 bridge (using barcode med strict): %ld", cnt);
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
	log_debug("Number of joined n-m bridge: %ld", cnt);
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
	log_debug("Number of joined n-m node: %ld", cnt);
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
	log_debug("Number of joined 1-1 pair(s) through jungle: %ld", ret);
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

void do_some_resolve_bridge(struct asm_graph_t *g)
{
	gint_t cnt = 0, cnt_local;
	do {
		cnt_local = 0;
		cnt_local += resolve_2_2_bridge_high_strict(g);
		cnt_local += resolve_2_2_bridge_med_strict(g);
		cnt += cnt_local;
	} while (cnt_local);
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
				// ret += join_n_m_complex_jungle(g, set_e, set_leg, set_self, uni_cov);
			}
		}
		kh_clear(gint, set_leg);
		kh_clear(gint, set_e);
		kh_clear(gint, set_v);
		kh_clear(gint, set_self);
	}
	log_debug("Number of joined pair(s) through jungle: %ld", ret);
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
	for (i = k = 0; i < (int) e->seq_len; ++i) {
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
		barcode = unpack_int64((uint8_t *) buf);
		k = kh_put(bcpos, h, barcode, &ret);
		if (ret != 1)
			__ERROR("Insert barcode failed");
		kh_value(h, k).r1_offset = unpack_int64((uint8_t *) buf + 8);
		kh_value(h, k).r2_offset = unpack_int64((uint8_t *) buf + 16);
		kh_value(h, k).r1_len = unpack_int64((uint8_t *) buf + 24);
		kh_value(h, k).r2_len = unpack_int64((uint8_t *) buf + 32);
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

void test_same_barcode(struct read_path_t *ref, khash_t(bcpos) *dict,
					   uint64_t *shared, int n_shared)
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
		check_data(buf, pos[i].r1_len);
		fseek(fi2, pos[i].r2_offset, SEEK_SET);
		xfread(buf, 1, pos[i].r2_len, fi2);
		check_data(buf, pos[i].r2_len);
	}
	fclose(fi1);
	fclose(fi2);
	free(buf);
	free(pos);
}

/**
 * @brief: gets reads to build the local assembly graph
 * @param reads: path to the original read files
 * @param rpath: path to the local reads after being constructed
 * @param dict: barcodes offset in file
 * @param g: the global graph
 * @param e1, e2: the two bridge contigs
 * @param prefix: prefix of the paths to the local reads
 * @return: 1 if everything runs smoothly, 0 if something supicious happens
 */
int get_reads_local_graph(struct read_path_t *reads, struct read_path_t *rpath,
			khash_t(bcpos) *dict, struct asm_graph_t *g,
			gint_t e1, gint_t e2, const char *prefix)
{
	char path[MAX_PATH];
	sprintf(path, "%s/R1.sub.fq", prefix);
	rpath->R1_path = strdup(path);
	sprintf(path, "%s/R2.sub.fq", prefix);
	rpath->R2_path = strdup(path);
	rpath->idx_path = NULL;
	struct barcode_hash_t *bc1, *bc2;
	bc1 = g->edges[e1].barcodes + 2;
	bc2 = g->edges[e2].barcodes + 2;
	khash_t(gint) *h1 = barcode_hash_2_khash(bc1);
	khash_t(gint) *h2 = barcode_hash_2_khash(bc2);

	int e1_rc = g->edges[e1].rc_id;
	int e2_rc = g->edges[e2].rc_id;
	khash_t(gint) *h_head = barcode_hash_2_khash(g->edges[e1_rc].barcodes + 2);
	khash_t(gint) *h_tail = barcode_hash_2_khash(g->edges[e2_rc].barcodes + 2);
	khash_t(gint) *h_exclude = get_union_bc(h_head, h_tail);
	kh_destroy(gint, h_head);
	kh_destroy(gint, h_tail);

	khash_t(gint) *h_shared = get_shared_bc(h1, h2);
	khash_t(gint) *h_union = get_union_bc(h1, h2);

	khash_t(gint) *include = kh_init(gint);
	for (int i = 0; i < g->n_e; ++i){
		int e = i;
		int e_rc = g->edges[e].rc_id;
		khash_t(gint) *h = barcode_hash_2_khash(g->edges[e].barcodes + 2);
		khash_t(gint) *h_rc = barcode_hash_2_khash(g->edges[e_rc].barcodes + 2);
		khash_t(gint) *h_u = get_union_bc(h, h_rc);

		if (kh_size(get_shared_bc(h_u, h_shared)) == 0)
			continue;

		for (khiter_t it = kh_begin(h_u); it != kh_end(h_u); ++it){
			if (!kh_exist(h_u, it))
				continue;
			uint64_t key = kh_key(h_u, it);
			if (!kh_gint_exist(h1, key) && !kh_gint_exist(h2, key))
				continue;
			if (kh_gint_exist(h_exclude, key))
				continue;
			kh_gint_add(include, key);
		}
		kh_destroy(gint, h);
		kh_destroy(gint, h_rc);
		kh_destroy(gint, h_u);
	}

	uint64_t *shared;
	int n;
	khash_2_arr(include, &shared, &n);
	log_debug("Keep %d barcodes in %d total", n, kh_size(h_union));
	filter_read(reads, dict, rpath, shared, n);
	free(shared);
	kh_destroy(gint, h1);
	kh_destroy(gint, h2);
	kh_destroy(gint, h_shared);
	kh_destroy(gint, h_exclude);
	kh_destroy(gint, h_union);

	int res = 1;
	if (check_file_empty(rpath->R1_path) || check_file_empty(rpath->R2_path))
		res = 0;
	return res;
}

/**
 * @brief: gets reads to build the estimated coverage for each edge in the local graph
 * @param reads: path to the original read files
 * @param rpath: path to the local reads after being constructed
 * @param dict: barcodes offset in file
 * @param g: the global graph
 * @param e1, e2: the two bridge contigs
 * @param prefix: prefix of the paths to the local reads
 * @return: 1 if everything runs smoothly, 0 if something supicious happens
 */
int get_reads_build_cov(struct read_path_t *reads, struct read_path_t *rpath,
			khash_t(bcpos) *dict, struct asm_graph_t *g,
			gint_t e1, gint_t e2, const char *prefix)
{
	char path[MAX_PATH];
	sprintf(path, "%s/R1.sub.fq", prefix);
	rpath->R1_path = strdup(path);
	sprintf(path, "%s/R2.sub.fq", prefix);
	rpath->R2_path = strdup(path);
	rpath->idx_path = NULL;
	struct barcode_hash_t *bc1, *bc2;
	bc1 = &(g->edges[e1].barcodes_cov);
	bc2 = &(g->edges[e2].barcodes_cov);

	struct barcode_hash_t *bc_head;
	struct barcode_hash_t *bc_tail;
	bc_head = g->edges[g->edges[e1].rc_id].barcodes + 2;
	bc_tail = g->edges[g->edges[e2].rc_id].barcodes + 2;

	khash_t(gint) *h1 = barcode_hash_2_khash(bc1);
	khash_t(gint) *h2 = barcode_hash_2_khash(bc2);

	khash_t(gint) *h_shared = get_shared_bc(h1, h2);
	uint64_t *shared;
	int n;
	khash_2_arr(h_shared, &shared, &n);
	log_debug("Keep %d barcodes in %d total", n, kh_size(h_shared));
	filter_read(reads, dict, rpath, shared, n);
	free(shared);
	kh_destroy(gint, h1);
	kh_destroy(gint, h2);
	kh_destroy(gint, h_shared);

	int res = 1;
	if (check_file_empty(rpath->R1_path) || check_file_empty(rpath->R2_path))
		res = 0;
	return res;
}

void get_local_reads(struct read_path_t *reads, struct read_path_t *rpath,
			khash_t(bcpos) *dict, struct asm_graph_t *g,
			gint_t e1, gint_t e2, const char *prefix)
{
	const int contig_level = 1;
	get_union_reads(reads, rpath, dict, g, e1, e2, prefix, contig_level);
}

/**
 * @brief: gets the shared reads of two contigs
 * @param reads: path to the original read files
 * @param rpath: path to the local reads after being constructed
 * @param dict: barcodes offset in file
 * @param g: the global graph
 * @param e1, e2: the two bridge contigs
 * @param prefix: prefix of the paths to the local reads
 * @param contig_level: the level of the barcodes that are being retrived,
 * 	based on CONTIG_LEVEL_0, CONTIG_LEVEL_1, CONTIG_LEVEL_2
 */
void get_shared_reads(struct read_path_t *reads, struct read_path_t *rpath,
		khash_t(bcpos) *dict, struct asm_graph_t *g, gint_t e1,
		gint_t e2, const char *prefix, int contig_level)
{
	char path[MAX_PATH];
	sprintf(path, "%s/R1.sub.fq", prefix);
	rpath->R1_path = strdup(path);
	sprintf(path, "%s/R2.sub.fq", prefix);
	rpath->R2_path = strdup(path);
	rpath->idx_path = NULL;
	struct barcode_hash_t *bc1, *bc2;
	bc1 = g->edges[e1].barcodes + contig_level;
	bc2 = g->edges[e2].barcodes + contig_level;

	khash_t(gint) *h1 = barcode_hash_2_khash(bc1);
	khash_t(gint) *h2 = barcode_hash_2_khash(bc2);
	khash_t(gint) *h_shared = get_shared_bc(h1, h2);
	uint64_t *shared;
	int n;
	khash_2_arr(h_shared, &shared, &n);
	filter_read(reads, dict, rpath, shared, n);
	free(shared);
	kh_destroy(gint, h1);
	kh_destroy(gint, h2);
	kh_destroy(gint, h_shared);
}

/**
 * @brief: gets the union reads of two contigs
 * @param reads: path to the original read files
 * @param rpath: path to the local reads after being constructed
 * @param dict: barcodes offset in file
 * @param g: the global graph
 * @param e1, e2: the two bridge contigs
 * @param prefix: prefix of the paths to the local reads
 * @param contig_level: the level of the barcodes that are being retrived,
 * 	based on CONTIG_LEVEL_0, CONTIG_LEVEL_1, CONTIG_LEVEL_2
 */
void get_union_reads(struct read_path_t *reads, struct read_path_t *rpath,
		khash_t(bcpos) *dict, struct asm_graph_t *g, gint_t e1,
		gint_t e2, const char *prefix, int contig_level)
{
	char path[MAX_PATH];
	sprintf(path, "%s/R1.sub.fq", prefix);
	rpath->R1_path = strdup(path);
	sprintf(path, "%s/R2.sub.fq", prefix);
	rpath->R2_path = strdup(path);
	rpath->idx_path = NULL;
	struct barcode_hash_t *bc1, *bc2;
	bc1 = g->edges[e1].barcodes + contig_level;
	bc2 = g->edges[e2].barcodes + contig_level;

	khash_t(gint) *h1 = barcode_hash_2_khash(bc1);
	khash_t(gint) *h2 = barcode_hash_2_khash(bc2);
	khash_t(gint) *h_union = get_union_bc(h1, h2);
	uint64_t *res;
	int n;
	khash_2_arr(h_union, &res, &n);
	filter_read(reads, dict, rpath, res, n);
	free(res);
	kh_destroy(gint, h1);
	kh_destroy(gint, h2);
	kh_destroy(gint, h_union);
}

/**
 * @brief: gets reads to check when finding paths
 * @param reads: path to the original read files
 * @param rpath: path to the local reads after being constructed
 * @param dict: barcodes offset in file
 * @param g: the global graph
 * @param e1, e2: the two bridge contigs
 * @param prefix: prefix of the paths to the local reads
 */
int get_reads_kmer_check(struct opt_proc_t *opt, struct asm_graph_t *g,
		int e1, int e2, struct read_path_t *local_read_path)
{
	struct read_path_t read_sorted_path;
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
        	log_error("Reads must be sorted\n");
	}
	khash_t(bcpos) *dict = kh_init(bcpos);
	construct_read_index(&read_sorted_path, dict);

	char work_dir[MAX_PATH];
	sprintf(work_dir, "%s/kmer_check_%d_%d", opt->out_dir, e1, e2);
	mkdir(work_dir, 0755);
	const int contig_level = 2;
	get_shared_reads(&read_sorted_path, local_read_path, dict, g,
			g->edges[e1].rc_id, e2, work_dir, contig_level);
	kh_destroy(bcpos, dict);
	int res = 1;
	if (check_file_empty(local_read_path->R1_path) ||
		check_file_empty(local_read_path->R2_path))
		res = 0;
	return res;
}

int mark_visit(struct asm_graph_t *g, gint_t u, gint_t t, uint8_t *flag, uint8_t mask)
{
	if (u == t)
		return 1;
	int ret = 0, next_ret;
	gint_t j, e, v;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		if (flag[e] & mask)
			continue;
		flag[e] |= mask;
		v = g->edges[e].target;
		next_ret = mark_visit(g, v, t, flag, mask);
		ret |= next_ret;
	}
	return ret;
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

struct seq_path_t {
	int len;
	gint_t *seq;
	int e0_beg;
	int en_end;
	int seq_beg;
	int seq_end;
};

struct ref_pos_t {
	gint_t id;
	uint32_t pos;
};

struct alg_seed_t {
	gint_t ref;
	uint32_t ref_pos;
	uint32_t seq_pos;
};

KHASH_MAP_INIT_INT64(graph_index, struct ref_pos_t);

void build_local_graph_index(struct asm_graph_t *g, khash_t(graph_index) *dict,  uint32_t ksize)
{
	gint_t e, e_rc;
	uint32_t k, c;
	uint64_t kmer, krev, kmask;
	int hash_ret;
	kmask = ((uint64_t)1 << (ksize << 1)) - 1;
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		// if (e > e_rc)
		// 	continue;
		kmer = krev = 0;
		for (k = 0; k < g->edges[e].seq_len; ++k) {
			c = __binseq_get(g->edges[e].seq, k);
			kmer = ((kmer << 2) & kmask) | c;
			krev = (krev >> 2) | ((uint64_t)(c ^ 3) << ((ksize - 1) << 1));
			khiter_t it;
			if (k + 1 >= ksize) {
				it = kh_put(graph_index, dict, kmer, &hash_ret);
				// if (kmer <= krev) {
				// 	it = kh_put(graph_index, dict, kmer, &hash_ret);
				// } else {
				// 	it = kh_put(graph_index, dict, krev, &hash_ret);
				// }
				if (hash_ret == 0)
					kh_value(dict, it) = (struct ref_pos_t){(gint_t)-1, 0};
				else
					kh_value(dict, it) = (struct ref_pos_t){e, k};
			}
		}
	}
}

void find_path_on_graph(struct asm_graph_t *g, khash_t(graph_index) *dict,
		int ksize, uint32_t *seq, int len, struct seq_path_t *ret)
{
	uint64_t kmask, kmer, krev;
	uint32_t c;
	int i, k, n_seed;
	khiter_t it;
	struct alg_seed_t *aseed = malloc(len * sizeof(struct alg_seed_t));
	n_seed = 0;
	kmask = ((uint64_t)1 << (ksize << 1)) - 1;
	kmer = krev = 0;
	for (k = 0; k < len; ++k) {
		c = __binseq_get(seq, k);
		kmer = ((kmer << 2) & kmask) | c;
		krev = (krev >> 2) | ((uint64_t)(c ^ 3) << ((ksize - 1) << 1));
		if (k + 1 >= ksize) {
			it = kh_get(graph_index, dict, kmer);
			// if (kmer <= krev) {
			// 	it = kh_get(graph_index, dict, kmer);
			// } else {
			// 	it = kh_get(graph_index, dict, krev);
			// }
			if (it == kh_end(dict))
				continue;
			if (kh_value(dict, it).id != -1) {
				aseed[n_seed++] = (struct alg_seed_t){kh_value(dict, it).id,
						kh_value(dict, it).pos, k};
			} else {
			}
		}
	}
	if (n_seed == 0) { /* find no seed hit */
		ret->len = 0;
		ret->seq = NULL;
		goto find_path_clean;
	}
	ret->e0_beg = aseed[0].ref_pos + 1 - ksize;
	ret->en_end = aseed[n_seed - 1].ref_pos;
	ret->seq_beg = aseed[0].seq_pos + 1 - ksize;
	ret->seq_end = aseed[n_seed - 1].seq_pos;
	ret->seq = NULL;
	ret->len = 0;
	gint_t prev_e = -1;
	int cnt_fw, cnt_rv;
	cnt_fw = cnt_rv = 0;
	for (i = 0; i < n_seed; ++i) {
		if (aseed[i].ref != prev_e) {
			ret->seq = realloc(ret->seq, (ret->len + 1) * sizeof(gint_t));
			ret->seq[ret->len++] = aseed[i].ref;
		} else {
			if (aseed[i].ref_pos > aseed[i - 1].ref_pos)
				++cnt_fw;
			else
				++cnt_rv;
		}
		prev_e = aseed[i].ref;
	}
	// if (cnt_rv > cnt_fw) {
	// 	gint_t *new_seq = malloc(ret->len * sizeof(gint_t));
	// 	for (i = 0; i < ret->len; ++i)
	// 		new_seq[i] = g->edges[ret->seq[ret->len - i - 1]].rc_id;
	// 	ret->e0_beg = g->edges[aseed[n_seed - 1].ref].seq_len - aseed[n_seed - 1].ref_pos - 1;
	// 	ret->en_end = g->edges[aseed[0].ref].seq_len - (aseed[0].ref_pos + 1 - ksize) - 1;
	// 	free(ret->seq);
	// 	ret->seq = new_seq;
	// }

find_path_clean:
	free(aseed);
}

void get_sub_edge(uint32_t *ref, int l, int r, uint32_t **seq, int *len)
{
	if (l >= r) {
		*seq = NULL;
		*len = r - l;
		return;
	}
	*len = r - l;
	*seq = calloc((*len + 15) >> 4, sizeof(uint32_t));
	int i, k;
	uint32_t c;
	for (k = l, i = 0; k < r; ++k, ++i) {
		c = __binseq_get(ref, k);
		__binseq_set(*seq, i, c);
	}
}

void print_one_path(struct asm_graph_t *g, gint_t e1, gint_t e2, uint32_t id,
	uint32_t n_e, gint_t *list_e, uint32_t *cand_len, FILE *fp)
{
	fprintf(fp, ">%u\n", id);
	char *buf = malloc(81);
	buf[80] = '\0';
	int buf_len = 0, i, k, e1_len, e2_len;
	uint32_t sum_len = 0;
	gint_t e, j;
	e1_len = __min(g->edges[e1].seq_len, 200);
	for (i = e1_len; i > 0; --i) {
		k = g->edges[e1].seq_len - i;
		buf[buf_len++] = nt4_char[__binseq_get(g->edges[e1].seq, k)];
		if (buf_len == 80) {
			fprintf(fp, "%s\n", buf);
			buf_len = 0;
		}
	}
	for (i = 0; (uint32_t) i < n_e; ++i) {
		e = list_e[i];
		sum_len += (g->edges[e].seq_len - g->ksize);
		for (k = g->ksize; (uint32_t) k < g->edges[e].seq_len; ++k) {
			buf[buf_len++] = nt4_char[__binseq_get(g->edges[e].seq, k)];
			if (buf_len == 80) {
				fprintf(fp, "%s\n", buf);
				buf_len = 0;
			}
		}
	}
	log_debug("Candidate %u, len: %u", id, sum_len + g->ksize);
	e2_len = __min(g->edges[e2].seq_len, 200);
	for (k = g->ksize; k < e2_len; ++k) {
		buf[buf_len++] = nt4_char[__binseq_get(g->edges[e2].seq, k)];
		if (buf_len == 80) {
			fprintf(fp, "%s\n", buf);
			buf_len = 0;
		}
	}
	if (buf_len) {
		buf[buf_len] = '\0';
		fprintf(fp, "%s\n", buf);
	}
	cand_len[id] = e1_len + sum_len + e2_len - g->ksize;
}

/**
 * @brief: get the local assembly graph between two contigs
 * @param opt: application options
 * @param g: the glbal graph
 * @param e1, e2: the two bridge contigs
 * @param dict: barcode posittions
 * @description:
 * 	build the local graph for a specific region between two contigs
 * 	given the original read files, the global graph, the id of the two
 * 	bridge contigs, it needs to:
 * 		+ get the local reads
 * 		+ build the local assembly graph using those reads
 *
 * @expected outcome: if the graph is built succesfully, a file named
 * 	graph_k_%d_local_level_1.bin must be placed in the working directory,
 * 	othewise such a file must not be written.
 *
 * @pre-conditions:
 * 	Original read files must be sorted
 * 
 * @post-conditions:
 * 	if read files are empty, no assembly graph file must exist in the working dir,
 * 	otherwise such a file must exist
 *
 * @potential bugs:
 * 	Local read files are empty
 */
void get_local_assembly(struct opt_proc_t *opt, struct asm_graph_t *g,
					gint_t e1, gint_t e2, khash_t(bcpos) *dict)
{
	pre_test(get_local_assembly, opt);
	struct read_path_t read_sorted_path = parse_read_path_from_opt(opt);

	struct read_path_t local_read_path;
	char work_dir_build_graph[MAX_PATH];
	sprintf(work_dir_build_graph, "%s/local_assembly_%ld_%ld", opt->out_dir,
			e1, e2);
	mkdir(work_dir_build_graph, 0755);

	char work_dir_build_cov[MAX_PATH];
	sprintf(work_dir_build_cov, "%s/build_cov_%ld_%ld", opt->out_dir, e1, e2);
	mkdir(work_dir_build_cov, 0755);

	struct read_path_t read_build_graph, read_build_cov;
	int ret_graph = get_reads_local_graph(&read_sorted_path, &read_build_graph,
			dict, g, e1, e2, work_dir_build_graph);
	int ret_cov = get_reads_build_cov(&read_sorted_path, &read_build_cov,
			dict, g, e1, e2, work_dir_build_cov);
	if (!ret_graph || !ret_cov){
		log_warn("Something supicious happends, probably the read files are empty, please check these files: %s %s",
				read_build_graph.R1_path, read_build_graph.R2_path);
	} else {
		struct asm_graph_t lg, lg1;
		build_local_assembly_graph(opt->lk, opt->n_threads, opt->mmem, 1,
			&(read_build_graph.R1_path), &(read_build_graph.R2_path),
			work_dir_build_graph, &lg, g, e1, e2);
		save_graph_info(work_dir_build_graph, &lg, "local_lvl_0");
		build_local_0_1(&lg, &lg1);
		build_local_graph_cov(opt, g, &lg1, e1, e2, work_dir_build_cov,
				read_build_cov.R1_path, read_build_cov.R2_path);
		save_graph_info(work_dir_build_graph, &lg1, "local_lvl_1");
		destroy_read_path(&read_build_graph);
		asm_graph_destroy(&lg1);

	}
	post_test(get_local_assembly, opt->lk, work_dir_build_graph, ret_graph && ret_cov);
}

KHASH_INIT(used_pair, struct pair_contig_t, char, 0, __mix_2_64, __cmp_2_64);
KHASH_MAP_INIT_INT64(cap_set, int);

struct dfs_listpath_t {
	gint_t e1;
	gint_t e2;
	uint32_t max_len;
	uint32_t max_cnt;
	FILE *fp;
	uint32_t cnt;
	uint32_t *cand_len;
	uint32_t n_e;
	gint_t *list_e;
	double uni_cov;
};

void list_all_path(struct asm_graph_t *g, gint_t u, gint_t t, uint32_t cur_len,
	khash_t(gint) *set_e, khash_t(cap_set) *cap, struct dfs_listpath_t *opt)
{
	if (opt->cnt > opt->max_cnt || cur_len > opt->max_len)
		return;
	if (u == t) {
		print_one_path(g, opt->e1, opt->e2, opt->cnt, opt->n_e,
					opt->list_e, opt->cand_len, opt->fp);
		++opt->cnt;
	}
	gint_t e, j;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		if (kh_get(gint, set_e, e) == kh_end(set_e))
			continue;
		khiter_t it = kh_get(cap_set, cap, e);
		if (kh_value(cap, it)) {
			--kh_value(cap, it);
			opt->list_e[opt->n_e++] = e;
			list_all_path(g, g->edges[e].target, t,
				cur_len + g->edges[e].seq_len - g->ksize, set_e,
				cap, opt);
			--opt->n_e;
			++kh_value(cap, it);
		}
	}
}

void list_all_path_local(struct asm_graph_t *g, gint_t u, gint_t t, uint32_t cur_len,
			khash_t(cap_set) *cap, struct dfs_listpath_t *opt)
{
	if (opt->cnt > opt->max_cnt || cur_len > opt->max_len)
		return;
	if (u == t) {
		print_one_path(g, opt->e1, opt->e2, opt->cnt, opt->n_e,
					opt->list_e, opt->cand_len, opt->fp);
		++opt->cnt;
	}
	gint_t e, j;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		khiter_t it = kh_get(cap_set, cap, e);
		if (kh_value(cap, it)) {
			--kh_value(cap, it);
			opt->list_e[opt->n_e++] = e;
			list_all_path_local(g, g->edges[e].target, t,
				cur_len + g->edges[e].seq_len - g->ksize,
				cap, opt);
			--opt->n_e;
			++kh_value(cap, it);
		}
	}
}

void list_path_slim_local(struct asm_graph_t *g, gint_t u, gint_t t, uint32_t cur_len,
			khash_t(cap_set) *cap, struct dfs_listpath_t *opt)
{
	if (cur_len > opt->max_len)
		return;
	if (u == t) {
		++opt->cnt;
	}
	if (opt->cnt > opt->max_cnt)
		return;
	gint_t e, j;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		khiter_t it = kh_get(cap_set, cap, e);
		if (kh_value(cap, it)) {
			--kh_value(cap, it);
			opt->list_e[opt->n_e++] = e;
			list_path_slim_local(g, g->edges[e].target, t,
				cur_len + g->edges[e].seq_len - g->ksize,
				cap, opt);
			if (opt->cnt > opt->max_cnt)
				return;
			--opt->n_e;
			++kh_value(cap, it);
		}
	}
}

void list_path_slim(struct asm_graph_t *g, gint_t u, gint_t t, uint32_t cur_len,
	khash_t(gint) *set_e, khash_t(cap_set) *cap, struct dfs_listpath_t *opt)
{
	if (cur_len > opt->max_len)
		return;
	if (u == t) {
		++opt->cnt;
	}
	if (opt->cnt > opt->max_cnt)
		return;
	gint_t e, j;
	for (j = 0; j < g->nodes[u].deg; ++j) {
		e = g->nodes[u].adj[j];
		if (kh_get(gint, set_e, e) == kh_end(set_e))
			continue;
		khiter_t it = kh_get(cap_set, cap, e);
		if (kh_value(cap, it)) {
			--kh_value(cap, it);
			opt->list_e[opt->n_e++] = e;
			list_path_slim(g, g->edges[e].target, t,
				cur_len + g->edges[e].seq_len - g->ksize, set_e,
				cap, opt);
			if (opt->cnt > opt->max_cnt)
				return;
			--opt->n_e;
			++kh_value(cap, it);
		}
	}
}

void reconstruct_path(struct asm_graph_t *g, gint_t *list_e, uint32_t n_e,
	struct result_local_t *sret)
{
	uint32_t len = g->ksize, i, k, c;
	gint_t e;
	for (i = 0; i < n_e; ++i) {
		e = list_e[i];
		len += (g->edges[e].seq_len - g->ksize);
	}
	sret->len = len;
	sret->seq = calloc((len + 15) >> 4, sizeof(uint32_t));
	for (k = 0; k < (uint32_t) g->ksize; ++k) {
		c = __binseq_get(g->edges[list_e[0]].seq, k);
		__binseq_set(sret->seq, k, c);
	}
	len = g->ksize;
	for (i = 0; i < n_e; ++i) {
		e = list_e[i];
		for (k = g->ksize; k < g->edges[e].seq_len; ++k) {
			c = __binseq_get(g->edges[e].seq, k);
			__binseq_set(sret->seq, len, c);
			++len;
		}
	}
}

void reconstruct_path_global(struct asm_graph_t *g, gint_t e1, gint_t e2,
	uint32_t e1_len, uint32_t e2_len, gint_t *list_e, uint32_t n_e,
	struct result_local_t *sret)
{
	__VERBOSE_LOG("DEBUG", "n_e = %u\n", n_e);
	uint32_t len, i, k, c;
	gint_t e;
	len = e1_len;
	for (i = 0; i < n_e; ++i) {
		e = list_e[i];
		len += (g->edges[e].seq_len - g->ksize);
	}
	len += (e2_len - g->ksize);
	sret->len = len;
	sret->seq = calloc((len + 15) >> 4, sizeof(uint32_t));
	for (k = 0; k < e1_len; ++k) {
		c = __binseq_get(g->edges[e1].seq, g->edges[e1].seq_len - e1_len + k);
		__binseq_set(sret->seq, k, c);
	}
	len = e1_len;
	for (i = 0; i < n_e; ++i) {
		e = list_e[i];
		for (k = g->ksize; k < g->edges[e].seq_len; ++k) {
			c = __binseq_get(g->edges[e].seq, k);
			__binseq_set(sret->seq, len, c);
			++len;
		}
	}
	for (k = g->ksize; k < e2_len; ++k) {
		c = __binseq_get(g->edges[e2].seq, k);
		__binseq_set(sret->seq, len, c);
		++len;
	}
}

void correct_local_path(struct asm_graph_t *g, struct result_local_t *sret)
{
	if (sret->len >= g->ksize * 2) {
		sret->trim_e1 = sret->trim_e2 = g->ksize;
	} else {
		sret->len -= g->ksize * 2;
		sret->trim_e1 = sret->trim_e2 = 0;
		free(sret->seq);
		sret->seq = NULL;
	}
}

int find_best_local_path(struct opt_local_t *opt, struct asm_graph_t *g,
			gint_t e1, gint_t e2, uint32_t e1_len, uint32_t e2_len,
			struct result_local_t *sret)
{
	gint_t e, max_list_e, e1_rc, e2_rc;
	int ret = -1, hash_ret, i;
	khiter_t it, k;
	char fasta_path[MAX_PATH];
	khash_t(cap_set) *cap = kh_init(cap_set);
	for (e = 0; e < g->n_e; ++e) {
		it = kh_put(cap_set, cap, e, &hash_ret);
		kh_value(cap, it) = 4;
	}
	max_list_e = g->n_e * 4;

	sprintf(fasta_path, "%s/candidates.fasta", opt->out_dir);
	__VERBOSE_LOG("DEBUG", "fasta path: %s\n", fasta_path);

	struct dfs_listpath_t dfs_opt;
	dfs_opt.e1 = e1;
	dfs_opt.e2 = e2;
	dfs_opt.max_len = 3000;
	dfs_opt.max_cnt = 500;
	dfs_opt.fp = xfopen(fasta_path, "wb");
	dfs_opt.cnt = 0;
	dfs_opt.cand_len = malloc((dfs_opt.max_cnt + 1) * sizeof(uint32_t));
	dfs_opt.n_e = 0;
	dfs_opt.list_e = malloc(max_list_e * sizeof(gint_t));
	list_all_path_local(g, g->edges[e1].target, g->edges[e2].source, 0, cap,
								&dfs_opt);
	fclose(dfs_opt.fp);
	if (dfs_opt.cnt == 0) {
		ret = -1;
		goto resolve_failed;
	}
	if (dfs_opt.cnt > dfs_opt.max_cnt) {
		ret = 2;
		goto resolve_failed;
	}
	khash_t(contig_count) *ctg_cnt = kh_init(contig_count);
	for (i = 0; (uint32_t) i < dfs_opt.cnt; ++i) {
		k = kh_put(contig_count, ctg_cnt, i, &hash_ret);
		kh_value(ctg_cnt, k) = 0;
	}
	count_readpair_path(opt->n_threads, opt->read_path, fasta_path, ctg_cnt);

	log_debug( "Testing path on local graph");
	int best_cnt, best_len, best_cand;
	best_cnt = 0;
	best_len = 0;
	best_cand = -1;
	for (k = kh_begin(ctg_cnt); k != kh_end(ctg_cnt); ++k) {
		if (!kh_exist(ctg_cnt, k))
			continue;
		if (kh_value(ctg_cnt, k) > best_cnt) {
			best_cnt = kh_value(ctg_cnt, k);
			best_len = dfs_opt.cand_len[kh_key(ctg_cnt, k)];
			best_cand = kh_key(ctg_cnt, k);
		} else if (kh_value(ctg_cnt, k) == best_cnt &&
			dfs_opt.cand_len[kh_key(ctg_cnt, k)] > (uint32_t) best_len) {
			best_len = dfs_opt.cand_len[kh_key(ctg_cnt, k)];
			best_cand = kh_key(ctg_cnt, k);
		}
		log_debug( "Candidate %ld with %d counts!", kh_key(ctg_cnt, k),
			kh_value(ctg_cnt, k));
	}
	dfs_opt.max_cnt = best_cand;
	dfs_opt.cnt = 0;
	dfs_opt.n_e = 0;

	list_path_slim_local(g, g->edges[e1].target, g->edges[e2].source, 0, cap,
								&dfs_opt);
	sret->seq = NULL;
	sret->len = sret->trim_e1 = sret->trim_e2 = 0;
	reconstruct_path_global(g, e1, e2, e1_len, e2_len, dfs_opt.list_e, dfs_opt.n_e, sret);
	kh_destroy(contig_count, ctg_cnt);
	ret = 1;

resolve_failed:
	kh_destroy(cap_set, cap);
	free(dfs_opt.cand_len);
	free(dfs_opt.list_e);
	return ret;
}

/* return:
 * -1: no path
 *  1: found good path
 *  2: have path but not pick one
 */
int fill_path_local(struct opt_local_t *opt, struct asm_graph_t *g0, struct asm_graph_t *g,
			gint_t e1, gint_t e2, struct result_local_t *sret)
{
	int ksize, ret, i;
	ret = 0;
	log_debug( "Find path local assembly %ld(%ld) <-> %ld(%ld)",
		g0->edges[e1].rc_id, e1, e2, g0->edges[e2].rc_id);
	khash_t(graph_index) *kdict = kh_init(graph_index);
	ksize = __min(g->ksize, 31);
	build_local_graph_index(g, kdict, ksize);
	gint_t e1_rc = g0->edges[e1].rc_id;
	struct seq_path_t ep1, ep2;
	log_debug("Find path for e1");
	find_path_on_graph(g, kdict, ksize, g0->edges[e1_rc].seq, g0->edges[e1_rc].seq_len, &ep1);
	if (ep1.len == 0) {
		log_debug( "e1: null path");
	} else {
		log_debug("ep1.len = %d", ep1.len);
		log_debug( "e0_beg = %d; en_end = %d; seq_beg = %d; seq_end = %d",
			ep1.e0_beg, ep1.en_end, ep1.seq_beg, ep1.seq_end);
		for (i = 0; i < ep1.len; ++i) {
			fprintf(stderr, "%ld, ", ep1.seq[i]);
		}
		fprintf(stderr, "\n");
	}
	log_debug("Find path for e2");
	find_path_on_graph(g, kdict, ksize, g0->edges[e2].seq, g0->edges[e2].seq_len, &ep2);
	if (ep2.len == 0) {
		log_debug("e2: null path");
	} else {
		log_debug("ep2.len = %d", ep2.len);
		log_debug("e0_beg = %d; en_end = %d; seq_beg = %d; seq_end = %d",
			ep2.e0_beg, ep2.en_end, ep2.seq_beg, ep2.seq_end);
		for (i = 0; i < ep2.len; ++i) {
			fprintf(stderr, "%ld, ", ep2.seq[i]);
		}
		fprintf(stderr, "\n");
	}
	if (ep1.len == 0 || ep2.len == 0) {
		ret = 0;
		goto fill_path_clean;
	}
	if (ep1.seq[ep1.len - 1] == ep2.seq[0]) {
		sret->trim_e1 = g0->edges[e1].seq_len - ep1.seq_end - 1;
		sret->trim_e2 = ep2.seq_beg;
		get_sub_edge(g->edges[ep2.seq[0]].seq, ep1.en_end + 1, ep2.e0_beg, &(sret->seq), &(sret->len));
		log_debug("trim_e1 = %d; trim_e2 = %d; len = %d",
			sret->trim_e1, sret->trim_e2, sret->len);
		log_debug("Easy peasy case %ld <-> %ld: success", e1, e2);
		ret = 1;
	} else {
		sret->trim_e1 = g0->edges[e1].seq_len - ep1.seq_end - 1;
		sret->trim_e2 = ep2.seq_beg;
		int sub_e1, sub_e2;
		sub_e1 = g->edges[ep1.seq[ep1.len - 1]].seq_len - ep1.en_end - 1;
		sub_e2 = ep2.e0_beg;
		if (sub_e1 < g->ksize) {
			sret->trim_e1 += g->ksize - sub_e1;
			sub_e1 = g->ksize;
		}
		if (sub_e2 < g->ksize) {
			sret->trim_e2 += g->ksize - sub_e2;
			sub_e2 = g->ksize;
		}
		ret = find_best_local_path(opt, g, ep1.seq[ep1.len - 1], ep2.seq[0],
			sub_e1, sub_e2, sret);
		if (ret == 1)
			log_debug("Found path: %ld <-> %ld", e1, e2);
		else if (ret == 2)
			log_debug("Too many path:  %ld <-> %ld", e1, e2);
		else if (ret == -1)
			log_debug("No path:  %ld <-> %ld", e1, e2);
	}
fill_path_clean:
	free(ep1.seq);
	free(ep2.seq);
	kh_destroy(graph_index, kdict);
	return ret;
}


int join_1_1_jungle_la(struct asm_graph_t *g, khash_t(gint) *set_e,
			khash_t(gint) *set_leg, struct opt_local_t *opt,
			khash_t(used_pair) *assemblied_pair)
{
	gint_t *legs = alloca(2 * sizeof(gint_t));
	int n_leg, hash_ret, max_list_e, i, ret;
	ret = 0;
	khiter_t k, it;
	gint_t e1, e2, e1_rc, e2_rc, e;
	n_leg = 0;
	for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
		if (!kh_exist(set_leg, k))
			continue;
		legs[n_leg++] = kh_key(set_leg, k);
	}
	e1 = legs[0];
	e2 = legs[1];
	e1_rc = g->edges[e1].rc_id;
	e2_rc = g->edges[e2].rc_id;
	char work_dir[MAX_PATH], fasta_path[MAX_PATH];
	sprintf(work_dir, "%s/construct_bridge_local_%ld_%ld", opt->out_dir, e1, e2);
	mkdir(work_dir, 0755);
	struct read_path_t local_read;
	get_local_reads(opt->read_path, &local_read, opt->dict, g, e1, e2, work_dir);

	double uni_cov = (__get_edge_cov(g->edges + e1, g->ksize) +
				__get_edge_cov(g->edges + e2, g->ksize)) / 2;
	khash_t(cap_set) *cap = kh_init(cap_set);
	max_list_e = 0;
	for (k = kh_begin(set_e); k != kh_end(set_e); ++k) {
		if (!kh_exist(set_e, k))
			continue;
		e = kh_key(set_e, k);
		struct cov_range_t rcov = get_edge_cov_range(g, e, uni_cov);
		max_list_e += rcov.hi;
		it = kh_put(cap_set, cap, e, &hash_ret);
		kh_value(cap, it) = rcov.hi;
	}

	sprintf(fasta_path, "%s/candidates.fasta", work_dir);

	struct dfs_listpath_t dfs_opt;
	dfs_opt.e1 = e1_rc;
	dfs_opt.e2 = e2;
	dfs_opt.max_len = 3000;
	dfs_opt.max_cnt = 500;
	dfs_opt.fp = xfopen(fasta_path, "wb");
	dfs_opt.cnt = 0;
	dfs_opt.cand_len = malloc((dfs_opt.max_cnt + 1) * sizeof(uint32_t));
	dfs_opt.n_e = 0;
	dfs_opt.list_e = malloc(max_list_e * sizeof(gint_t));
	dfs_opt.uni_cov = uni_cov;
	list_all_path(g, g->edges[e1_rc].target, g->edges[e2].source, 0, set_e,
		cap, &dfs_opt);
	fclose(dfs_opt.fp);
	if (dfs_opt.cnt == 0 || dfs_opt.cnt > dfs_opt.max_cnt) {
		ret = 0;
		goto resolve_failed;
	}
	khash_t(contig_count) *ctg_cnt = kh_init(contig_count);
	for (i = 0; (uint32_t) i < dfs_opt.cnt; ++i) {
		k = kh_put(contig_count, ctg_cnt, i, &hash_ret);
		kh_value(ctg_cnt, k) = 0;
	}
	count_readpair_path(opt->n_threads, &local_read, fasta_path, ctg_cnt);
	destroy_read_path(&local_read);

	log_debug( "Testing local path %ld <-> %ld", e1, e2);
	int best_cnt, best_len, best_cand;
	best_cnt = 0;
	best_len = 0;
	best_cand = -1;
	for (k = kh_begin(ctg_cnt); k != kh_end(ctg_cnt); ++k) {
		if (!kh_exist(ctg_cnt, k))
			continue;
		if (kh_value(ctg_cnt, k) > best_cnt) {
			best_cnt = kh_value(ctg_cnt, k);
			best_len = dfs_opt.cand_len[kh_key(ctg_cnt, k)];
			best_cand = kh_key(ctg_cnt, k);
		} else if (kh_value(ctg_cnt, k) == best_cnt &&
			dfs_opt.cand_len[kh_key(ctg_cnt, k)] > (uint32_t) best_len) {
			best_len = dfs_opt.cand_len[kh_key(ctg_cnt, k)];
			best_cand = kh_key(ctg_cnt, k);
		}
		log_debug( "Candidate %ld with %d counts!", kh_key(ctg_cnt, k),
			kh_value(ctg_cnt, k));
	}
	dfs_opt.max_cnt = best_cand;
	dfs_opt.cnt = 0;
	dfs_opt.n_e = 0;
	list_path_slim(g, g->edges[e1_rc].target, g->edges[e2].source, 0, set_e,
		cap, &dfs_opt);
	struct result_local_t sret;
	reconstruct_path(g, dfs_opt.list_e, dfs_opt.n_e, &sret);
	correct_local_path(g, &sret);
	log_debug( "Middle part len: %d", sret.len);
	asm_join_edge_with_fill(g, e1_rc, e1, e2, e2_rc,
		sret.seq, sret.len, sret.trim_e1, sret.trim_e2);
	free(sret.seq);
	kh_destroy(contig_count, ctg_cnt);
	ret = 1;

resolve_failed:
	kh_destroy(cap_set, cap);
	free(dfs_opt.cand_len);
	free(dfs_opt.list_e);
	return ret;
}

gint_t barcode_find_pair_alter(struct asm_graph_t *g, khash_t(gint) *set_leg, gint_t se, gint_t ae, int *stat)
{
	gint_t e, ret_e, sec_e;
	khiter_t k;
	ret_e = ae;
	sec_e = -1;
	for (k = kh_begin(set_leg); k != kh_end(set_leg); ++k) {
		if (!kh_exist(set_leg, k))
			continue;
		e = kh_key(set_leg, k);
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
	if (ret_e == -1) {
		return -1;
	}
	if (sec_e != -1 && !check_barcode_superior(g, se, ret_e, sec_e)) {
		*stat = 0;
	} else {
		*stat = 1;
	}
	return ret_e;
}

void get_rc_seq(uint32_t **seq, uint32_t *ref, uint32_t len)
{
	uint32_t i, k, c;
	*seq = calloc((len + 15) >> 4, sizeof(uint32_t));
	for (i = 0; i < len; ++i) {
		k = len - i - 1;
		c = __binseq_get(ref, k) ^ 3;
		__binseq_set(*seq, i, c);
	}
}

void do_something_local(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	struct read_path_t read_sorted_path;
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
	sprintf(work_dir, "%s/level_4_temp_dir", opt->out_dir);
	mkdir(work_dir, 0755);

	double uni_cov = get_genome_coverage(g);
	khash_t(gint) *visited, *set_e, *set_v, *set_leg, *set_self;
	visited = kh_init(gint);
	set_e = kh_init(gint);
	set_v = kh_init(gint);
	set_leg = kh_init(gint);
	set_self = kh_init(gint);
	gint_t e, ret, resolve_local;
	uint32_t n_leg, n_self;
	khash_t(used_pair) *assemblied_pair = kh_init(used_pair);
	ret = 0;
	struct opt_local_t opt_local;
	opt_local.out_dir = work_dir;
	opt_local.read_path = &read_sorted_path;
	opt_local.dict = dict;
	opt_local.ksize = g->ksize;
	opt_local.n_threads = opt->n_threads;
	opt_local.mmem = opt->mmem;
	do {
		resolve_local = 0;
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			uint32_t len = get_edge_len(g->edges + e);
			if (kh_get(gint, visited, e) != kh_end(visited) || len < MIN_CONTIG_BARCODE)
				continue;
			find_region(g, e, MIN_CONTIG_BARCODE, MAX_EDGE_COUNT, uni_cov, set_v, set_e);
			if (kh_size(set_e) < MAX_EDGE_COUNT) {
				kh_merge_set(visited, set_e);
				detect_leg(g, MIN_LONG_CONTIG, 6000,
						set_v, set_e, set_leg, set_self);
				n_leg = kh_size(set_leg);
				n_self = kh_size(set_self);
				if (n_self == 0 && n_leg == 2){
					resolve_local += join_1_1_jungle_la(g, set_e, set_leg,
						&opt_local, assemblied_pair);
				} else if (n_self + n_leg >= 2){
					/*resolve_local += join_n_m_complex_jungle_la(g, set_e, set_leg, set_self,
						&opt_local, assemblied_pair);*/
				}
			}
			kh_clear(gint, set_leg);
			kh_clear(gint, set_e);
			kh_clear(gint, set_v);
			kh_clear(gint, set_self);
		}
		kh_clear(gint, visited);
		ret += resolve_local;
	} while (resolve_local);
	kh_destroy(used_pair, assemblied_pair);
	kh_destroy(gint, set_leg);
	kh_destroy(gint, set_e);
	kh_destroy(gint, set_v);
	kh_destroy(gint, set_self);
	log_debug("Number of joined pair(s) through jungle: %ld", ret);
}

khash_t(gint) *barcode_hash_2_khash(struct barcode_hash_t *bc)
{
	khash_t(gint) *res = kh_init(gint);
	for (uint32_t i = 0; i < bc->size; ++i) {
		if (bc->keys[i] != (uint64_t) -1) {
			int ret;
			kh_put(gint, res, bc->keys[i], &ret);
		}
	}
	return res;
}

khash_t(gint) *get_shared_bc(khash_t(gint) *h1, khash_t(gint) *h2)
{
	khash_t(gint) *res = kh_init(gint);
	for (khiter_t it = kh_begin(h1); it != kh_end(h1); ++it){
		if (!kh_exist(h1, it))
			continue;
		gint_t key = kh_key(h1, it);
		khiter_t it2 = kh_get(gint, h2, key);
		if (it2 != kh_end(h2)){
			int ret;
			kh_put(gint, res, key, &ret);
		}
	}
	return res;
}

khash_t(gint) *get_union_bc(khash_t(gint) *h1, khash_t(gint) *h2)
{
	khash_t(gint) *res = kh_init(gint);
	for (khiter_t it = kh_begin(h1); it != kh_end(h1); ++it) {
		if (!kh_exist(h1, it))
			continue;
		int ret;
		kh_put(gint, res, kh_key(h1, it), &ret);
	}
	for (khiter_t it = kh_begin(h2); it != kh_end(h2); ++it) {
		if (!kh_exist(h2, it))
			continue;
		khiter_t it2 = kh_get(gint, h1, kh_key(h2, it));
		if (it2 == kh_end(h1)) {
			int ret;
			kh_put(gint, res, kh_key(h2, it), &ret);
		}
	}
	return res;
}

khash_t(gint) *get_exclude_bc(khash_t(gint) *h_in, khash_t(gint) *h_ex)
{
	khash_t(gint) *res = kh_init(gint);
	for (khiter_t it = kh_begin(h_in); it != kh_end(h_in); ++it) {
		if (!kh_exist(h_in, it))
			continue;
		khiter_t it2 = kh_get(gint, h_ex, kh_key(h_in, it));
		if (it2 == kh_end(h_ex)) {
			int ret;
			kh_put(gint, res, kh_key(h_in, it), &ret);
		}
	}
	return res;
}

void khash_2_arr(khash_t(gint) *h, uint64_t **arr, int *n)
{
	*n = 0;
	*arr = calloc(kh_size(h), sizeof(uint64_t));
	for (khiter_t it = kh_begin(h); it != kh_end(h); ++it) {
		if (!kh_exist(h, it))
			continue;
		(*arr)[(*n)++] = kh_key(h, it);
	}
}


