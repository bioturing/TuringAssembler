#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "barcode_hash.h"
#include "fastq_producer.h"
#include "io_utils.h"
#include "khash.h"
#include "kmer_build.h"
#include "kseq.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"
#include "scaffolding/global_params.h"
#include "../include/kmc_skipping.h"

KSEQ_INIT(gzFile, gzread);

__KHASH_IMPL(pair_contig_count, , struct pair_contig_t, struct contig_count_t, 1,
							__mix_2_64, __cmp_2_64);

__KHASH_IMPL(contig_count, , gint_t, int, 1, kh_int64_hash_func, kh_int64_hash_equal);

static inline int is_hole_rc(struct asm_edge_t *e1, struct asm_edge_t *e2)
{
	if (e1->n_holes != e2->n_holes)
		return 0;
	uint32_t i, len;
	len = e1->seq_len;
	for (i = 0; i < e1->n_holes; ++i) {
		if (e1->l_holes[i] != e2->l_holes[e2->n_holes - i - 1])
			return 0;
		if (e1->p_holes[i] != len - e2->p_holes[e2->n_holes - i - 1] - 2)
			return 0;
	}
	return 1;
}

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

gint_t get_longest_edge(struct asm_graph_t *g)
{
	gint_t e, ret = -1;
	uint32_t max_len = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		uint32_t len = get_edge_len(g->edges + e);
		if (len > max_len) {
			max_len = len;
			ret = e;
		}
	}
	return ret;
}

static uint32_t get_hash_edge32(struct asm_edge_t *e)
{
	uint32_t ret = 0, i;
	for (i = 0; i < ((e->seq_len + 15) >> 4); ++i)
		ret ^= e->seq[i];
	return ret;
}

double get_genome_coverage(struct asm_graph_t *g)
{
	/* Using the coverage of the longest contigs */
	gint_t e;
	double ret_cov = 0.0;
	uint32_t max_len = 0;
	uint32_t sum_len = 0;
	double sum_cov = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		if (g->edges[e].seq_len > max_len) {
			max_len = g->edges[e].seq_len;
			ret_cov = __get_edge_cov(g->edges + e, g->ksize);
		}
	}
	return ret_cov;
}

double get_genome_coverage_h(struct asm_graph_t *g)
{
	/* Using the coverage of the longest contigs */
	gint_t e;
	double ret_cov = 0.0;
	uint32_t max_len = 0;
	uint32_t sum_len = 0;
	double sum_cov = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		int len = get_edge_len(&g->edges[e]);
		float cov = __get_edge_cov(g->edges +e, g->ksize);
		if (len < 1000)
			continue;
		sum_len += g->edges[e].seq_len;
		sum_cov += g->edges[e].seq_len * cov;
		VERBOSE_FLAG(0, "sumlen %d sumcov %lf\n", sum_len, sum_cov);
	}
	VERBOSE_FLAG(0, "sumlen %d sumcov %lf\n", sum_len, sum_cov);
	return sum_cov/sum_len;
}

gint_t dump_edge_seq_h(char **seq, uint32_t *m_seq, struct asm_edge_t *e)
{
	uint32_t i, j, k, len = e->seq_len;
	for (i = 0; i < e->n_holes; ++i)
		len += e->l_holes[i];
	if (*m_seq < len + 1) {
		*m_seq = len + 1;
		*seq = realloc(*seq, *m_seq);
	}
	j = k = 0;
	for (i = 0; i < e->seq_len; ++i) {
		(*seq)[k++] = nt4_char[__binseq_get(e->seq, i)];
		/* append holes to the sequences */
		if (j < e->n_holes && e->p_holes[j] == i) {
			/* fill with 'N's */
			memset((*seq) + k, 0x4e, e->l_holes[j]);
			k += e->l_holes[j];
			++j;
		}
	}
	(*seq)[k] = '\0';
	return (gint_t)k;
}

// void asm_duplicate_edge_seq(struct asm_graph_t *g, gint_t e, int cov)
// {
// 	if (cov == 1)
// 		return;
// 	double fcov;
// 	fcov = __get_edge_cov(g->edges + e, g->ksize);
// 	// __VERBOSE("Edge %ld[len=%u][cov~%.3lf] before duplicate x%d, ", e,
// 	// 			get_edge_len(g->edges + e), fcov, cov);

// 	g->edges = realloc(g->edges, (g->n_e + 1) * sizeof(struct asm_edge_t));
// 	asm_clone_edge2(g, g->n_e, e);
// 	int n;
// 	for (n = 0; n + 1 < cov; ++n)
// 		asm_append_seq2(g, e, g->n_e);
// 	asm_clean_edge_seq(g->edges + g->n_e);

// 	fcov = __get_edge_cov(g->edges + e, g->ksize);
// 	// __VERBOSE("[len=%u][cov~%.3lf] after duplicate\n",
// 	// 		get_edge_len(g->edges + e), fcov);
// }

// void asm_duplicate_edge_seq2(struct asm_graph_t *g, gint_t e1, gint_t e2, int cov)
// {
// 	double fcov1, fcov2;
// 	fcov1 = __get_edge_cov(g->edges + e1, g->ksize);
// 	fcov2 = __get_edge_cov(g->edges + e2, g->ksize);
// 	// __VERBOSE("Edge e1=%ld[len=%u][cov~%.3lf], e2=%ld[len=%u][cov~%.3lf] before double duplicate x%d, ",
// 	// 	e1, get_edge_len(g->edges + e1), fcov1,
// 	// 	e2, get_edge_len(g->edges + e2), fcov2, cov);

// 	g->edges = realloc(g->edges, (g->n_e + 1) * sizeof(struct asm_edge_t));
// 	asm_clone_edge2(g, g->n_e, e1);
// 	int n;
// 	for (n = 0; n < cov; ++n) {
// 		asm_append_seq2(g, e1, e2);
// 		asm_append_seq2(g, e1, g->n_e);
// 	}
// 	g->edges[e1].count += g->edges[e2].count;
// 	asm_clean_edge_seq(g->edges + g->n_e);

// 	fcov1 = __get_edge_cov(g->edges + e1, g->ksize);
// 	// __VERBOSE("[len=%u][cov~%.3lf] after duplicate\n",
// 	// 		get_edge_len(g->edges + e1), fcov1);
// }

gint_t asm_create_node(struct asm_graph_t *g)
{
	g->nodes = realloc(g->nodes, (g->n_v + 2) * sizeof(struct asm_node_t));
	g->nodes[g->n_v].rc_id = g->n_v + 1;
	g->nodes[g->n_v + 1].rc_id = g->n_v;
	g->nodes[g->n_v].adj = g->nodes[g->n_v + 1].adj = NULL;
	g->nodes[g->n_v].deg = g->nodes[g->n_v + 1].deg = 0;
	g->n_v += 2;
	return g->n_v - 2;
}

void asm_clone_seq(struct asm_edge_t *dst, struct asm_edge_t *src)
{
	dst->count = src->count;
	dst->seq_len = src->seq_len;
	dst->seq = calloc((dst->seq_len + 15) >> 4, sizeof(uint32_t));
	memcpy(dst->seq, src->seq,
		((dst->seq_len + 15) >> 4) * sizeof(uint32_t));
	dst->n_holes = src->n_holes;
	if (dst->n_holes) {
		dst->p_holes = malloc(dst->n_holes * sizeof(uint32_t));
		memcpy(dst->p_holes, src->p_holes, dst->n_holes * sizeof(uint32_t));
		dst->l_holes = malloc(dst->n_holes * sizeof(uint32_t));
		memcpy(dst->l_holes, src->l_holes, dst->n_holes * sizeof(uint32_t));
	} else {
		dst->p_holes = NULL;
		dst->l_holes = NULL;
	}
}

void asm_clone_seq_reverse(struct asm_edge_t *dst, struct asm_edge_t *src)
{
	dst->count = src->count;
	dst->seq_len = src->seq_len;
	dst->seq = calloc((dst->seq_len + 15) >> 4, sizeof(uint32_t));
	uint32_t i, k;
	for (i = 0; i < dst->seq_len; ++i) {
		k = dst->seq_len - i - 1;
		dst->seq[i >> 4] |= (uint32_t)(__binseq_get(src->seq, k) ^ 3)
							<< ((i & 15) << 1);
	}
	dst->n_holes = src->n_holes;
	dst->l_holes = malloc(dst->n_holes * sizeof(uint32_t));
	dst->p_holes = malloc(dst->n_holes * sizeof(uint32_t));
	for (i = 0; i < dst->n_holes; ++i) {
		dst->l_holes[i] = src->l_holes[dst->n_holes - i - 1];
		dst->p_holes[i] = dst->seq_len - 1
			- (src->p_holes[dst->n_holes - i - 1] + 1);
	}
}

void asm_clone_edge(struct asm_graph_t *g, gint_t dst, gint_t src)
{
	/* clone the topology */
	asm_clone_seq(g->edges + dst, g->edges + src);
	g->edges[dst].source = g->edges[src].source;
	g->edges[dst].target = g->edges[src].target;
	/* clone the barcode */
	if (g->aux_flag & ASM_HAVE_BARCODE) {
		g->edges[dst].barcodes = calloc(3, sizeof(struct barcode_hash_t));
		barcode_hash_clone(g->edges[dst].barcodes, g->edges[src].barcodes);
		barcode_hash_clone(g->edges[dst].barcodes + 1, g->edges[src].barcodes + 1);
		barcode_hash_clone(g->edges[dst].barcodes + 2, g->edges[src].barcodes + 2);
	}
}

gint_t asm_create_clone_edge(struct asm_graph_t *g, gint_t src)
{
	g->edges = realloc(g->edges, (g->n_e + 2) * sizeof(struct asm_edge_t));
	memset(g->edges + g->n_e, 0, 2 * sizeof(struct asm_edge_t));
	g->n_e += 2;
	asm_clone_edge(g, g->n_e - 2, src);
	asm_clone_edge(g, g->n_e - 1, g->edges[src].rc_id);
	g->edges[g->n_e - 2].rc_id = g->n_e - 1;
	g->edges[g->n_e - 1].rc_id = g->n_e - 2;
	asm_add_node_adj(g, g->edges[g->n_e - 2].source, g->n_e - 2);
	asm_add_node_adj(g, g->edges[g->n_e - 1].source, g->n_e - 1);
	return g->n_e - 2;
}

void asm_append_seq_with_gap(struct asm_edge_t *dst,
				struct asm_edge_t *src, uint32_t gap_size)
{
	/* append the bin seq */
	uint32_t seq_len, new_m, m, i, k;
	seq_len = dst->seq_len + src->seq_len;
	new_m = (seq_len + 15) >> 4;
	m = (dst->seq_len + 15) >> 4;
	if (new_m > m) {
		dst->seq = realloc(dst->seq, new_m * sizeof(uint32_t));
		memset(dst->seq + m, 0, (new_m - m) * sizeof(uint32_t));
	}
	for (i = 0; i < src->seq_len; ++i) {
		k = i + dst->seq_len;
		dst->seq[k >> 4] |= ((src->seq[i >> 4] >> ((i & 15) << 1) & 3)
							<< ((k & 15) << 1));
	}
	uint32_t n_holes = dst->n_holes + src->n_holes + 1;
	dst->p_holes = realloc(dst->p_holes, n_holes * sizeof(uint32_t));
	dst->l_holes = realloc(dst->l_holes, n_holes * sizeof(uint32_t));
	/* new gap */
	dst->p_holes[dst->n_holes] = dst->seq_len - 1;
	dst->l_holes[dst->n_holes] = gap_size;
	for (i = 0; i < src->n_holes; ++i)
		dst->p_holes[dst->n_holes + i + 1] = src->p_holes[i] + dst->seq_len;
	if (src->n_holes)
		memcpy(dst->l_holes + dst->n_holes + 1, src->l_holes, src->n_holes * sizeof(uint32_t));
	dst->n_holes = n_holes;
	dst->seq_len = seq_len;
}

void asm_append_barcode_edge(struct asm_edge_t *dst, struct asm_edge_t *src)
{
	if (dst->seq_len < CONTIG_LEVEL_0) {
		barcode_hash_merge(dst->barcodes, src->barcodes);
	}
	if (dst->seq_len < CONTIG_LEVEL_1) {
		if (dst->seq_len + CONTIG_LEVEL_0 >= CONTIG_LEVEL_1)
			barcode_hash_merge(dst->barcodes + 1, src->barcodes);
		else
			barcode_hash_merge(dst->barcodes + 1, src->barcodes + 1);
	}
	if (dst->seq_len < CONTIG_LEVEL_2) {
		if (dst->seq_len + CONTIG_LEVEL_0 >= CONTIG_LEVEL_2)
			barcode_hash_merge(dst->barcodes + 2, src->barcodes);
		else if (dst->seq_len + CONTIG_LEVEL_1 >= CONTIG_LEVEL_2)
			barcode_hash_merge(dst->barcodes + 2, src->barcodes + 1);
		else
			barcode_hash_merge(dst->barcodes + 2, src->barcodes + 2);
	}
}

void asm_append_barcode_readpair(struct asm_graph_t *g, gint_t dst, gint_t src)
{
	if (g->aux_flag & ASM_HAVE_BARCODE) {
		if (g->edges[dst].seq_len < CONTIG_LEVEL_0) {
			barcode_hash_merge(g->edges[dst].barcodes, g->edges[src].barcodes);
		}
		if (g->edges[dst].seq_len < CONTIG_LEVEL_1) {
			if (g->edges[dst].seq_len + CONTIG_LEVEL_0 >= CONTIG_LEVEL_1)
				barcode_hash_merge(g->edges[dst].barcodes + 1,
						g->edges[src].barcodes);
			else
				barcode_hash_merge(g->edges[dst].barcodes + 1,
						g->edges[src].barcodes + 1);
		}
		if (g->edges[dst].seq_len < CONTIG_LEVEL_2) {
			if (g->edges[dst].seq_len + CONTIG_LEVEL_0 >= CONTIG_LEVEL_2)
				barcode_hash_merge(g->edges[dst].barcodes + 2,
						g->edges[src].barcodes);
			else if (g->edges[dst].seq_len + CONTIG_LEVEL_1 >= CONTIG_LEVEL_2)
				barcode_hash_merge(g->edges[dst].barcodes + 2,
						g->edges[src].barcodes + 1);
			else
				barcode_hash_merge(g->edges[dst].barcodes + 2,
						g->edges[src].barcodes + 2);
		}
	}
}

void asm_append_seq_with_fill(struct asm_edge_t *dst, struct asm_edge_t *src,
			uint32_t *seq, int len, int trim_dst, int trim_src)
{
	uint32_t seq_len, new_m, m, c;
	int i, k;
	seq_len = dst->seq_len - trim_dst + len + src->seq_len - trim_src;
	new_m = (seq_len + 15) >> 4;
	m = (dst->seq_len - trim_dst + 15) >> 4;
	dst->seq = realloc(dst->seq, new_m * sizeof(uint32_t));
	if (dst->seq == NULL)
		__ERROR("Unable to realloc");
	if (new_m > m)
		memset(dst->seq + m, 0, (new_m - m) * sizeof(uint32_t));
	uint32_t mask = ((dst->seq_len - trim_dst) & 15) == 0 ? (uint32_t)-1 :
		(((uint32_t)1 << ((((dst->seq_len - trim_dst) & 15) << 1))) - 1);
	dst->seq[m - 1] &= mask;

	if (len >= 0) {
		for (i = 0, k = dst->seq_len - trim_dst; i < len; ++i, ++k) {
			c = __binseq_get(seq, i);
			__binseq_set(dst->seq, k, c);
		}

		for (i = trim_src, k = dst->seq_len - trim_dst + len; i < src->seq_len; ++i, ++k) {
			c = __binseq_get(src->seq, i);
			__binseq_set(dst->seq, k, c);
		}
	} else {
		for (i = trim_src - len, k = dst->seq_len - trim_dst; i < src->seq_len; ++i, ++k) {
			c = __binseq_get(src->seq, i);
			__binseq_set(dst->seq, k, c);
		}
	}
	dst->seq_len = seq_len;
	/* there is no gap at this stage */
}

void asm_append_seq_with_fill_reverse(struct asm_edge_t *dst, struct asm_edge_t *src,
			uint32_t *seq, int len, int trim_dst, int trim_src)
{
	uint32_t c, seq_len, new_m, m;
	int i, k;
	seq_len = dst->seq_len - trim_dst + len + src->seq_len - trim_src;
	new_m = (seq_len + 15) >> 4;
	m = (dst->seq_len - trim_dst + 15) >> 4;
	dst->seq = realloc(dst->seq, new_m * sizeof(uint32_t));
	if (dst->seq == NULL)
		__ERROR("Unable to realloc");
	if (new_m > m)
		memset(dst->seq + m, 0, (new_m - m) * sizeof(uint32_t));
	uint32_t mask = ((dst->seq_len - trim_dst) & 15) == 0 ? (uint32_t)-1 :
		(((uint32_t)1 << ((((dst->seq_len - trim_dst) & 15) << 1))) - 1);
	dst->seq[m - 1] &= mask;

	if (len >= 0) {
		for (i = len - 1, k = dst->seq_len - trim_dst; i >= 0; --i, ++k) {
			c = __binseq_get(seq, i) ^ 3;
			__binseq_set(dst->seq, k, c);
		}

		for (i = trim_src, k = dst->seq_len - trim_dst + len; i < src->seq_len; ++i, ++k) {
			c = __binseq_get(src->seq, i);
			__binseq_set(dst->seq, k, c);
		}
	} else {
		for (i = trim_src - len, k = dst->seq_len - trim_dst; i < src->seq_len; ++i, ++k) {
			c = __binseq_get(src->seq, i);
			__binseq_set(dst->seq, k, c);
		}
	}
	dst->seq_len = seq_len;
	/* there is no gap at this stage */
}

void asm_append_seq(struct asm_edge_t *dst, struct asm_edge_t *src, uint32_t overlap)
{
	uint32_t i, k;
	/* append the bin seq */
	uint32_t seq_len, new_m, m;
	seq_len = dst->seq_len + src->seq_len - overlap;
	new_m = (seq_len + 15) >> 4;
	m = (dst->seq_len + 15) >> 4;
	if (new_m > m) {
		dst->seq = realloc(dst->seq, new_m * sizeof(uint32_t));
		if (dst->seq == NULL)
			__ERROR("Unable to realloc");
		memset(dst->seq + m, 0, (new_m - m) * sizeof(uint32_t));
	}

	for (i = overlap; i < src->seq_len; ++i) {
		k = i - overlap + dst->seq_len;
		dst->seq[k >> 4] |= ((src->seq[i >> 4] >> ((i & 15) << 1)) & 3)
							<< ((k & 15) << 1);
	}
	/* append the gaps */
	if (src->n_holes) {
		uint32_t n_holes = dst->n_holes + src->n_holes;
		dst->p_holes = realloc(dst->p_holes, n_holes * sizeof(uint32_t));
		dst->l_holes = realloc(dst->l_holes, n_holes * sizeof(uint32_t));
		for (i = 0; i < src->n_holes; ++i) {
			dst->p_holes[dst->n_holes + i] = src->p_holes[i] + dst->seq_len - overlap;
		}
		memcpy(dst->l_holes + dst->n_holes, src->l_holes, src->n_holes * sizeof(uint32_t));
		dst->n_holes = n_holes;
	}
	dst->seq_len = seq_len;
}

void asm_join_edge_with_gap(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
				gint_t e2, gint_t e_rc2, uint32_t gap_size)
{
	/*    contig 1  |  gap  |   contig2
	 * AAAAAAAAAAAAANNNNNNNNNAAAAAAAAAAAA
	 * Since the barcode + read pair is now preserve only on the 1st contig
	 * we do not need to append the barcode + read pair information */
	asm_append_seq_with_gap(g->edges + e1, g->edges + e2, gap_size);
	g->edges[e1].target = g->edges[e2].target;
	g->edges[e1].count += g->edges[e2].count;

	asm_append_seq_with_gap(g->edges + e_rc2, g->edges + e_rc1, gap_size);
	g->edges[e_rc2].target = g->edges[e_rc1].target;
	g->edges[e_rc2].count += g->edges[e_rc1].count;

	g->edges[e1].rc_id = e_rc2;
	g->edges[e_rc2].rc_id = e1;

	asm_remove_edge(g, e2);
	asm_remove_edge(g, e_rc1);
}

void asm_join_edge_with_fill(struct asm_graph_t *g, gint_t e1, gint_t e_rc1, gint_t e2, gint_t e_rc2,
	uint32_t *aseq, int alen, int trim_e1, int trim_e2)
{
	asm_append_barcode_readpair(g, e1, e2);
	asm_append_seq_with_fill(g->edges + e1, g->edges + e2, aseq, alen, trim_e1, trim_e2);
	g->edges[e1].target = g->edges[e2].target;
	g->edges[e1].count += g->edges[e2].count;

	asm_append_barcode_readpair(g, e_rc2, e_rc1);
	asm_append_seq_with_fill_reverse(g->edges + e_rc2, g->edges + e_rc1, aseq, alen, trim_e2, trim_e1);
	g->edges[e_rc2].target = g->edges[e_rc1].target;
	g->edges[e_rc2].count += g->edges[e_rc1].count;

	g->edges[e1].rc_id = e_rc2;
	g->edges[e_rc2].rc_id = e1;

	asm_remove_edge(g, e2);
	asm_remove_edge(g, e_rc1);
}

void asm_join_edge(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
					gint_t e2, gint_t e_rc2)
{
	/*    contig 1  | overlap |
	 * ATCTTCGGTTTTTCTTTAAAAAAG
	 *              CTTTAAAAAAGAATTAAAAACTT
	 *                           contig 2
	 * Since the barcode + read pair is now preserve only on the 1st contig
	 * we do not need to append the barcode + read pair information */
	asm_append_barcode_readpair(g, e1, e2);
	asm_append_seq(g->edges + e1, g->edges + e2, g->ksize);
	g->edges[e1].target = g->edges[e2].target;
	g->edges[e1].count += g->edges[e2].count;

	asm_append_barcode_readpair(g, e_rc2, e_rc1);
	asm_append_seq(g->edges + e_rc2, g->edges + e_rc1, g->ksize);
	g->edges[e_rc2].target = g->edges[e_rc1].target;
	g->edges[e_rc2].count += g->edges[e_rc1].count;

	g->edges[e1].rc_id = e_rc2;
	g->edges[e_rc2].rc_id = e1;

	asm_remove_edge(g, e2);
	asm_remove_edge(g, e_rc1);
}

void asm_unroll_loop_forward(struct asm_graph_t *g, gint_t e1, gint_t e2, int rep)
{
	g->edges = realloc(g->edges, (g->n_e + 1) * sizeof(struct asm_edge_t));
	memset(g->edges + g->n_e, 0, sizeof(struct asm_edge_t));
	++g->n_e;
	asm_clone_edge(g, g->n_e - 1, e1);
	int i;
	for (i = 0; i < rep; ++i) {
		asm_append_barcode_readpair(g, e1, e2);
		asm_append_seq(g->edges + e1, g->edges + e2, g->ksize);
		asm_append_barcode_readpair(g, e1, g->n_e - 1);
		asm_append_seq(g->edges + e1, g->edges + (g->n_e - 1), g->ksize);
		if (g->edges[e1].seq_len >= 2000)
			break;
	}
	g->edges[e1].count += g->edges[e2].count;
	asm_clean_edge(g, g->n_e - 1);
	--g->n_e;
}

void asm_join_edge3(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
	gint_t e2, gint_t e_rc2, gint_t e3, gint_t e_rc3, uint64_t e2_count)
{
	/*    contig 1  | overlap |      | overlap |  contig 3
	 * ATCTTCGGTTTTTCTTTAAAAAAG      AAACTTTTTTTGGGGGATACCC
	 *              CTTTAAAAAAGAATTAAAAACTTTTTTT
	 *                        contig 2
	 * Since e2 is usually a repetitive edges, we need to pre-estimate the
	 * count that e2 contributes to final edge
	 */
	asm_append_barcode_readpair(g, e1, e2);
	asm_append_seq(g->edges + e1, g->edges + e2, g->ksize);
	asm_append_barcode_readpair(g, e1, e3);
	asm_append_seq(g->edges + e1, g->edges + e3, g->ksize);
	g->edges[e1].target = g->edges[e3].target;
	g->edges[e1].count += g->edges[e3].count + e2_count;

	asm_append_barcode_readpair(g, e_rc3, e_rc2);
	asm_append_seq(g->edges + e_rc3, g->edges + e_rc2, g->ksize);
	asm_append_barcode_readpair(g, e_rc3, e_rc1);
	asm_append_seq(g->edges + e_rc3, g->edges + e_rc1, g->ksize);
	g->edges[e_rc3].target = g->edges[e_rc1].target;
	g->edges[e_rc3].count += g->edges[e_rc1].count + e2_count;
	
	g->edges[e1].rc_id = e_rc3;
	g->edges[e_rc3].rc_id = e1;

	asm_remove_edge(g, e3);
	asm_remove_edge(g, e_rc1);
}

// void asm_join_edge_loop_reverse(struct asm_graph_t *g, gint_t e1, gint_t e2,
// 				gint_t e_rc2, gint_t e_rc1)
// {
// 	g->edges = realloc(g->edges, (g->n_e + 2) * sizeof(struct asm_edge_t));
// 	g->n_e += 2;
// 	asm_clone_edge2(g, g->n_e - 2, e2);
// 	asm_clone_edge2(g, g->n_e - 1, e_rc2);
// 	asm_append_seq2(g, g->n_e - 2, e_rc1);
// 	asm_append_seq2(g, g->n_e - 1, e_rc1);

// 	asm_remove_edge(g, e_rc1);
// 	asm_clone_edge2(g, e_rc1, e1);

// 	asm_append_seq2(g, e1, g->n_e - 1);
// 	asm_append_seq2(g, e_rc1, g->n_e - 2);

// 	g->edges[e1].count += g->edges[e2].count;
// 	g->edges[e_rc1].count = g->edges[e1].count;
// 	g->edges[e1].target = g->nodes[g->edges[e1].source].rc_id;
// 	asm_add_node_adj(g, g->edges[e_rc1].source, e_rc1);
// 	g->edges[e_rc1].target = g->nodes[g->edges[e_rc1].source].rc_id;

// 	asm_clean_edge_seq(g->edges + g->n_e - 2);
// 	asm_clean_edge_seq(g->edges + g->n_e - 1);
// 	g->edges[g->n_e - 2].source = g->edges[g->n_e - 2].target = -1;
// 	g->edges[g->n_e - 1].source = g->edges[g->n_e - 1].target = -1;
// 	g->n_e -= 2;
// }

// void asm_join_edge_loop(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
// 			gint_t e2, gint_t e_rc2, uint64_t added_count)
// {
// 	g->edges = realloc(g->edges, (g->n_e + 2) * sizeof(struct asm_edge_t));
// 	g->n_e += 2;
// 	asm_clone_edge2(g, g->n_e - 2, e2);
// 	asm_clone_edge2(g, g->n_e - 1, e_rc2);
// 	asm_append_seq2(g, g->n_e - 2, e1);
// 	asm_append_seq2(g, e1, g->n_e - 2);
// 	asm_append_seq2(g, g->n_e - 1, e_rc1);
// 	asm_append_seq2(g, e_rc1, g->n_e - 1);
// 	g->edges[e1].count += added_count;
// 	g->edges[e_rc1].count += added_count;
// 	asm_clean_edge_seq(g->edges + g->n_e - 2);
// 	asm_clean_edge_seq(g->edges + g->n_e - 1);
// 	g->edges[g->n_e - 2].source = g->edges[g->n_e - 2].target = -1;
// 	g->edges[g->n_e - 1].source = g->edges[g->n_e - 1].target = -1;
// 	g->n_e -= 2;
// }

static inline void asm_clean_edge_seq(struct asm_edge_t *e)
{
	free(e->seq);
	free(e->l_holes);
	free(e->p_holes);
	e->seq = NULL;
	e->l_holes = NULL;
	e->p_holes = NULL;
	e->n_holes = 0;
}

void asm_clean_edge(struct asm_graph_t *g, gint_t e)
{
	asm_clean_edge_seq(g->edges + e);
	if (g->aux_flag & ASM_HAVE_BARCODE) {
		barcode_hash_clean(g->edges[e].barcodes);
		barcode_hash_clean(g->edges[e].barcodes + 1);
		barcode_hash_clean(g->edges[e].barcodes + 2);
		free(g->edges[e].barcodes);
		barcode_hash_clean(&(g->edges[e].barcodes_scaf));
		barcode_hash_clean(&(g->edges[e].barcodes_scaf2));
	}
}

void asm_remove_edge(struct asm_graph_t *g, gint_t e)
{
	gint_t u = g->edges[e].source;
	if (u == -1) /* already clean */
		return;
	asm_clean_edge(g, e);
	gint_t j = find_adj_idx(g->nodes[u].adj, g->nodes[u].deg, e);
	if (j >= 0)
		g->nodes[u].adj[j] = g->nodes[u].adj[--g->nodes[u].deg];
	g->edges[e].source = g->edges[e].target = -1;
}

uint32_t get_edge_len(struct asm_edge_t *e)
{
	uint32_t ret = e->seq_len, i;
	for (i = 0; i < e->n_holes; ++i)
		ret += e->l_holes[i];
	return ret;
}

static void asm_edge_cc(struct asm_graph_t *g, gint_t *id_edge, gint_t **ret_size)
{
	memset(id_edge, 255, g->n_e * sizeof(gint_t));
	gint_t m_cc = 0x10000;
	gint_t *size = malloc(m_cc * sizeof(gint_t));
	gint_t n_cc = 0;
	gint_t k, l, r;
	gint_t *q = malloc(g->n_e * sizeof(gint_t));

	for (k = 0; k < g->n_e; ++k) {
		if (g->edges[k].source == -1)
			continue;
		if (id_edge[k] != -1)
			continue;
		id_edge[k] = id_edge[g->edges[k].rc_id] = n_cc;
		gint_t cur_size = 0;
		l = r = 0;
		q[0] = k;
		while (l <= r) {
			gint_t e = q[l++];
			gint_t e_rc = g->edges[e].rc_id;
			cur_size += 2 * (g->edges[e].seq_len - g->ksize);

			gint_t u, c;
			u = g->edges[e].target;
			if (g->nodes[u].deg == 0)
				cur_size += g->ksize;
			for (c = 0; c < g->nodes[u].deg; ++c) {
				gint_t next_e = g->nodes[u].adj[c];
				gint_t next_e_rc = g->edges[next_e].rc_id;
				if (id_edge[next_e] == -1) {
					id_edge[next_e] = id_edge[next_e_rc] = n_cc;
					q[++r] = next_e;
				}
			}

			u = g->edges[e_rc].target;
			if (g->nodes[u].deg == 0)
				cur_size += g->ksize;
			for (c = 0; c < g->nodes[u].deg; ++c) {
				gint_t next_e = g->nodes[u].adj[c];
				gint_t next_e_rc = g->edges[next_e].rc_id;
				if (id_edge[next_e] == -1) {
					id_edge[next_e] = id_edge[next_e_rc] = n_cc;
					q[++r] = next_e;
				}
			}
		}
		if (n_cc + 1 == m_cc) {
			m_cc <<= 1;
			size = realloc(size, m_cc * sizeof(gint_t));
		}
		size[n_cc++] = cur_size / 2;
	}
	size = realloc(size, n_cc * sizeof(gint_t));
	*ret_size = size;
	free(q);
}

static inline uint64_t get_bandage_count(struct asm_edge_t *e, int ksize)
{
	float cov = __get_edge_cov(e, ksize);
	uint32_t i, len = e->seq_len;
	for (i = 0; i < e->n_holes; ++i)
		len += e->l_holes[i];
	return (uint64_t)(cov * len);
}

static void print_debug(struct asm_edge_t *e)
{
	__DEBUG("seq_len = %u\n", e->seq_len);
	__DEBUG("n_holes = %u\n", e->n_holes);
	uint32_t i;
	for (i = 0; i < e->n_holes; ++i) {
		__DEBUG("p_holes = %u; l_holes = %u\n",
				e->p_holes[i], e->l_holes[i]);
	}
}

void write_fasta(struct asm_graph_t *g, const char *path)
{
	gint_t *id_edge, *cc_size;
	id_edge = malloc(g->n_e * sizeof(gint_t));
	cc_size = NULL;
	asm_edge_cc(g, id_edge, &cc_size);

	FILE *fp = xfopen(path, "w");
	char *seq = NULL;
	uint32_t seq_len = 0;
	char *buf = alloca(81);
	gint_t e, e_rc;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		gint_t cc_id = id_edge[e];
		if (cc_size[cc_id] < MIN_CONNECT_SIZE ||
			g->edges[e].seq_len < MIN_NOTICE_LEN)
			continue;
		gint_t len = dump_edge_seq_h(&seq, &seq_len, g->edges + e);
		fprintf(fp, ">SEQ_%lld_%lld_length_%lld_cov_%.3lf\n",
			(long long)e, (long long)e_rc, (long long)len,
			__get_edge_cov(g->edges + e, g->ksize));
		gint_t k = 0;
		while (k < len) {
			gint_t l = __min(80, len - k);
			memcpy(buf, seq + k, l);
			buf[l] = '\0';
			fprintf(fp, "%s\n", buf);
			k += l;
		}
	}
	fclose(fp);
	free(seq);
	free(id_edge);
	free(cc_size);
}

void write_gfa(struct asm_graph_t *g, const char *path)
{
	gint_t *id_edge, *cc_size;
	id_edge = malloc(g->n_e * sizeof(gint_t));
	cc_size = NULL;
	asm_edge_cc(g, id_edge, &cc_size);

	FILE *fp = xfopen(path, "w");
	char *seq = NULL;
	uint32_t seq_len = 0;
	gint_t e, e_rc;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		gint_t cc_id = id_edge[e];
		if (cc_size[cc_id] < 250)
			continue;
		dump_edge_seq_h(&seq, &seq_len, g->edges + e);
		uint64_t fake_count = get_bandage_count(g->edges + e, g->ksize);
		/* print fake count for correct coverage display on Bandage */
		fprintf(fp, "S\t%lld_%lld\t%s\tKC:i:%llu\n", (long long)e,
			(long long)e_rc, seq, (long long unsigned)fake_count);
	}
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1)
			continue;
		gint_t cc_id = id_edge[e];
		if (cc_size[cc_id] < 250)
			continue;
		e_rc = g->edges[e].rc_id;
		gint_t pe, pe_rc, next_pe, next_pe_rc;
		char ce, next_ce;
		if (e > e_rc) {
			pe = e_rc;
			pe_rc = e;
			ce = '-';
		} else {
			pe = e;
			pe_rc = e_rc;
			ce = '+';
		}
		gint_t n = g->edges[e].target;
		gint_t k;
		for (k = 0; k < g->nodes[n].deg; ++k) {
			gint_t next_e, next_e_rc;
			next_e = g->nodes[n].adj[k];
			next_e_rc = g->edges[next_e].rc_id;
			if (next_e > next_e_rc) {
				next_pe = next_e_rc;
				next_pe_rc = next_e;
				next_ce = '-';
			} else {
				next_pe = next_e;
				next_pe_rc = next_e_rc;
				next_ce = '+';
			}
			fprintf(fp, "L\t%lld_%lld\t%c\t%lld_%lld\t%c\t%dM\n",
				(long long)pe, (long long)pe_rc, ce,
				(long long)next_pe, (long long)next_pe_rc,
				next_ce, g->ksize);
		}
	}
	fclose(fp);
	free(seq);
	free(id_edge);
	free(cc_size);
}

void test2_asm_graph(struct asm_graph_t *g)
{
	gint_t u, e;
	for (u = 0; u < g->n_v; ++u) {
		gint_t j;
		for (j = 0; j < g->nodes[u].deg; ++j) {
			e = g->nodes[u].adj[j];
			gint_t e_rc = g->edges[e].rc_id;
			if (!is_seq_rc(g->edges[e].seq, g->edges[e].seq_len,
					g->edges[e_rc].seq, g->edges[e_rc].seq_len)) {
				fprintf(stderr, "seq_len = %u; rc_seq_len = %u\n",
					g->edges[e].seq_len, g->edges[e_rc].seq_len);
				assert(g->edges[e].seq_len == g->edges[e_rc].seq_len);
				char *seq = NULL;
				uint32_t lseq = 0;
				dump_edge_seq_h(&seq, &lseq, g->edges + e);
				fprintf(stderr, "seq    = %s\n", seq);
				dump_edge_seq_h(&seq, &lseq, g->edges + e_rc);
				fprintf(stderr, "seq_rc = %s\n", seq);
				assert(0 && "Smart error");
			}
		}
	}
}

static void debug_dump_adj(struct asm_graph_t *g, gint_t u)
{
	gint_t k, e;
	char *seq = NULL;
	uint32_t lseq = 0;
	for (k = 0; k < g->nodes[u].deg; ++k) {
		e = g->nodes[u].adj[k];
		dump_edge_seq_h(&seq, &lseq, g->edges + e);
		__VERBOSE("e = %ld: %s\n", e, seq);
	}
	free(seq);
}

void deb_dump_seq(struct asm_graph_t *g, gint_t e)
{
	uint32_t len = 0;
	char *seq = NULL;
	dump_edge_seq_h(&seq, &len, g->edges + e);
	printf("%s\n", seq);
}

void test_asm_graph(struct asm_graph_t *g)
{
	gint_t le_idx = get_longest_edge(g);
	if (le_idx != -1)
		__VERBOSE("Longest edge %ld_%ld, length %u\n",
			le_idx, g->edges[le_idx].rc_id, get_edge_len(g->edges + le_idx));

	gint_t u, e, j, k;
	for (u = 0; u < g->n_v; ++u) {
		/* Test 1: consistency of node's adj and edge's source */
		for (j = 0; j < g->nodes[u].deg; ++j) {
			e = g->nodes[u].adj[j];
			if (g->edges[e].source != u) {
				__VERBOSE("node = %ld; edge = [%ld](%ld->%ld)\n",
					u, e,
					g->edges[e].source, g->edges[e].target);
				assert(0 && "Node's adjs are node consistent with edges's source");
			}
		}
		/* Test 2: consistency of sequence of edges from one node */
		int c;
		for (c = 0; c < g->ksize; ++c) {
			int prev_nu = -1;
			for (k = 0; k < g->nodes[u].deg; ++k) {
				e = g->nodes[u].adj[k];
				int nu = __binseq_get(g->edges[e].seq, c);
				if (nu != prev_nu && prev_nu != -1) {
					debug_dump_adj(g, u);
					assert(0 && "Edges from same node not have same k-prefix");
				}
				prev_nu = nu;
			}
		}
		/* Test 3: Edge id within [0, g->n_e) */
		for (j = 0; j < g->nodes[u].deg; ++j) {
			e = g->nodes[u].adj[j];
			if (e < 0 || e >= g->n_e) {
				__VERBOSE("node = %ld; edge = %ld\n", u, e);
				assert(0 && "Node has undefined edges");
			}
		}
		/* Test 4: Node reverse complement id within [0, g->n_v) */
		if (g->nodes[u].rc_id < 0 || g->nodes[u].rc_id >= g->n_v) {
			__VERBOSE("node = %ld; rc_id = %ld\n", u, g->nodes[u].rc_id);
			assert(0 && "Node has undefined reverse complement");
		}
		/* Test 5: Continous edges share kmer */
		gint_t u_rc, e1, e2;
		u_rc = g->nodes[u].rc_id;
		if (g->nodes[u].deg > 0 && g->nodes[u_rc].deg > 0) {
			e1 = g->nodes[u].adj[0];
			e2 = g->edges[g->nodes[u_rc].adj[0]].rc_id;
			for (k = 0; k < g->ksize; ++k) {
				j = g->edges[e2].seq_len - g->ksize + k;
				if (__binseq_get(g->edges[e2].seq, j) !=
					__binseq_get(g->edges[e1].seq, k)) {
					printf("(%ld, %ld) -> (%ld, %ld)\n",
						g->edges[e2].source,
						g->edges[e2].target,
						g->edges[e1].source,
						g->edges[e1].target);
					deb_dump_seq(g, e2);
					deb_dump_seq(g, e1);
					assert(0 && "Continuous edges not share kmer");
				}
			}
		}
	}
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].source == -1) /* edge was removed */
			continue;
		gint_t e_rc = g->edges[e].rc_id;
		/* Test 1: Check correct reverse complement edge id */
		if (e_rc < 0 || e_rc >= g->n_e) {
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e, g->edges[e].source, g->edges[e].target,
				g->edges[e].rc_id);
			assert(0 && "Edge has undefined reverse complement");
		}
		if (e != g->edges[e_rc].rc_id) {
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e, g->edges[e].source, g->edges[e].target,
				g->edges[e].rc_id);
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e_rc, g->edges[e_rc].source, g->edges[e_rc].target,
				g->edges[e_rc].rc_id);
			assert(0 && "Edge reverse complement link is not 2-way");
		}
		/* Test 2: source and target within [0, g->n_e) */
		if (g->edges[e].source < 0 || g->edges[e].target >= g->n_v ||
			g->edges[e_rc].source < 0 || g->edges[e_rc].target >= g->n_v) {
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e, g->edges[e].source, g->edges[e].target,
				g->edges[e].rc_id);
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e_rc, g->edges[e_rc].source, g->edges[e_rc].target,
				g->edges[e_rc].rc_id);
			assert(0 && "Edge source and target node are undefined");
		}
		gint_t src, dst, src_rc, dst_rc;
		src = g->edges[e].source;
		dst = g->edges[e].target;
		/* Test 3: Edge must be in source's adj */
		gint_t idx = find_adj_idx(g->nodes[src].adj,
						g->nodes[src].deg, e);
		if (idx == -1) {
			__VERBOSE("node [%ld]; edge [%ld](%ld->%ld)\n",
				src, e, src, dst);
			assert(0 && "Edge not in source's adj");
		}
		src_rc = g->edges[e_rc].source;
		dst_rc = g->edges[e_rc].target;
		/* Test 4: Edge and its reverse complement must link reverse
		 * complemented nodes
		 */
		if (src != g->nodes[dst_rc].rc_id || dst != g->nodes[src_rc].rc_id) {
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e, g->edges[e].source, g->edges[e].target,
				g->edges[e].rc_id);
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e_rc, g->edges[e_rc].source, g->edges[e_rc].target,
				g->edges[e_rc].rc_id);
			assert(0 && "Edge and reverse complement not link between reverse complemented nodes");
		}
		/* Test 5: Sequence reverse complement */
		if (!is_seq_rc(g->edges[e].seq, g->edges[e].seq_len,
				g->edges[e_rc].seq, g->edges[e_rc].seq_len)) {
			char *seq = NULL;
			uint32_t lseq = 0;
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e, g->edges[e].source, g->edges[e].target,
				g->edges[e].rc_id);
			dump_edge_seq_h(&seq, &lseq, g->edges + e);
			// __VERBOSE("%s\n", seq);
			printf("seq_len = %lu; seq = %s\n", strlen(seq), seq);
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e_rc, g->edges[e_rc].source, g->edges[e_rc].target,
				g->edges[e_rc].rc_id);
			dump_edge_seq_h(&seq, &lseq, g->edges + e_rc);
			printf("seq_len = %lu; seq = %s\n", strlen(seq), seq);
			// __VERBOSE("%s\n", seq);
			assert(0 && "Edge and rc sequence is not reverse complemented");
		}
		if (!is_hole_rc(g->edges + e, g->edges + e_rc)) {
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e, g->edges[e].source, g->edges[e].target,
				g->edges[e].rc_id);
			__VERBOSE("n_holes = %u; seq_len = %u",
				g->edges[e].n_holes, g->edges[e].seq_len);
			uint32_t j;
			for (j = 0; j < g->edges[e].n_holes; ++j)
				__VERBOSE("(p=%u, l=%u) ",
					g->edges[e].p_holes[j],
					g->edges[e].l_holes[j]);
			__VERBOSE("\n");
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e_rc, g->edges[e_rc].source, g->edges[e_rc].target,
				g->edges[e_rc].rc_id);
			__VERBOSE("n_holes = %u; seq_len = %u",
				g->edges[e_rc].n_holes, g->edges[e_rc].seq_len);
			for (j = 0; j < g->edges[e_rc].n_holes; ++j)
				__VERBOSE("(p=%u, l=%u) ",
					g->edges[e_rc].p_holes[j],
					g->edges[e_rc].l_holes[j]);
			__VERBOSE("\n");
			assert(0 && "Edge and rc holes is not symmetric");
		}
	}

	// for (e = 0; e < g->n_e; ++e) {
	// 	uint32_t len = get_edge_len(g->edges + e);
	// 	double cov = __get_edge_cov(g->edges + e, g->ksize);
	// 	if (len > 5000 && cov < 100.0)
	// 		__VERBOSE("WARNING: Edge %ld has length %u with cov ~ %.6lf\n",
	// 			e, len, cov);
	// }
}

void save_asm_graph(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "wb");
	char *sig = "asmg";
	xfwrite(sig, 4, 1, fp);
	xfwrite(&g->aux_flag, 4, 1, fp);

	/* save the graph topology and sequence information */
	gint_t u, e;
	xfwrite(&g->ksize, sizeof(int), 1, fp);
	xfwrite(&g->n_v, sizeof(gint_t), 1, fp);
	xfwrite(&g->n_e, sizeof(gint_t), 1, fp);
	for (u = 0; u < g->n_v; ++u) {
		xfwrite(&g->nodes[u].rc_id, sizeof(gint_t), 1, fp);
		xfwrite(&g->nodes[u].deg, sizeof(gint_t), 1, fp);
		if (g->nodes[u].deg)
			xfwrite(g->nodes[u].adj, sizeof(gint_t), g->nodes[u].deg, fp);
	}
	for (e = 0; e < g->n_e; ++e) {
		xfwrite(&g->edges[e].source, sizeof(gint_t), 1, fp);
		xfwrite(&g->edges[e].target, sizeof(gint_t), 1, fp);
		if (g->edges[e].source == -1)
			continue;
		xfwrite(&g->edges[e].rc_id, sizeof(gint_t), 1, fp);
		xfwrite(&g->edges[e].count, sizeof(uint64_t), 1, fp);
		xfwrite(&g->edges[e].seq_len, sizeof(gint_t), 1, fp);
		xfwrite(g->edges[e].seq, sizeof(uint32_t), (g->edges[e].seq_len + 15) >> 4, fp);
		xfwrite(&g->edges[e].n_holes, sizeof(uint32_t), 1, fp);
		if (g->edges[e].n_holes) {
			xfwrite(g->edges[e].p_holes, sizeof(uint32_t), g->edges[e].n_holes, fp);
			xfwrite(g->edges[e].l_holes, sizeof(uint32_t), g->edges[e].n_holes, fp);
		}
	}

	/* save the barcode information */
	if (g->aux_flag & ASM_HAVE_BARCODE) {
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			struct barcode_hash_t *h = g->edges[e].barcodes;
			xfwrite(&h->size, sizeof(uint32_t), 1, fp);
			xfwrite(&h->n_item, sizeof(uint32_t), 1, fp);
			xfwrite(h->keys, sizeof(uint64_t), h->size, fp);
			// xfwrite(h->cnts, sizeof(uint32_t), h->size, fp);
			h = g->edges[e].barcodes + 1;
			xfwrite(&h->size, sizeof(uint32_t), 1, fp);
			xfwrite(&h->n_item, sizeof(uint32_t), 1, fp);
			xfwrite(h->keys, sizeof(uint64_t), h->size, fp);

			h = g->edges[e].barcodes + 2;
			xfwrite(&h->size, sizeof(uint32_t), 1, fp);
			xfwrite(&h->n_item, sizeof(uint32_t), 1, fp);
			xfwrite(h->keys, sizeof(uint64_t), h->size, fp);
		}
	}
	if (g->aux_flag & ASM_HAVE_BARCODE_SCAF) {
        for (e = 0; e < g->n_e; ++e) {
            if (g->edges[e].source == -1)
                continue;
            struct barcode_hash_t *h = &g->edges[e].barcodes_scaf;
            xfwrite(&h->size, sizeof(uint32_t), 1, fp);
            xfwrite(&h->n_item, sizeof(uint32_t), 1, fp);
            xfwrite(h->keys, sizeof(uint64_t), h->size, fp);
            h = &g->edges[e].barcodes_scaf2;
            xfwrite(&h->size, sizeof(uint32_t), 1, fp);
            xfwrite(&h->n_item, sizeof(uint32_t), 1, fp);
            xfwrite(h->keys, sizeof(uint64_t), h->size, fp);
        }
	}
	fclose(fp);
}

void load_asm_graph(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "rb");
	char sig[4];
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, "asmg", 4))
		__ERROR("Not assembly graph format file");
	xfread(&g->aux_flag, 4, 1, fp);
	__VERBOSE("aux_flag = %u\n", g->aux_flag);

	/* load the graph topology and sequence information */
	gint_t u, e;
	xfread(&g->ksize, sizeof(int), 1, fp);
	xfread(&g->n_v, sizeof(gint_t), 1, fp);
	xfread(&g->n_e, sizeof(gint_t), 1, fp);
	g->nodes = calloc(g->n_v, sizeof(struct asm_node_t));
	g->edges = calloc(g->n_e, sizeof(struct asm_edge_t));
	for (u = 0; u < g->n_v; ++u) {
		xfread(&g->nodes[u].rc_id, sizeof(gint_t), 1, fp);
		xfread(&g->nodes[u].deg, sizeof(gint_t), 1, fp);
		if (g->nodes[u].deg) {
			g->nodes[u].adj = malloc(g->nodes[u].deg * sizeof(gint_t));
			xfread(g->nodes[u].adj, sizeof(gint_t), g->nodes[u].deg, fp);
		} else {
			g->nodes[u].adj = NULL;
		}
	}
	for (e = 0; e < g->n_e; ++e) {
		xfread(&g->edges[e].source, sizeof(gint_t), 1, fp);
		xfread(&g->edges[e].target, sizeof(gint_t), 1, fp);
		if (g->edges[e].source == -1)
			continue;
		xfread(&g->edges[e].rc_id, sizeof(gint_t), 1, fp);
		xfread(&g->edges[e].count, sizeof(uint64_t), 1, fp);
		xfread(&g->edges[e].seq_len, sizeof(gint_t), 1, fp);
		g->edges[e].seq = calloc((g->edges[e].seq_len + 15) >> 4, sizeof(uint32_t));
		xfread(g->edges[e].seq, sizeof(uint32_t), (g->edges[e].seq_len + 15) >> 4, fp);
		xfread(&g->edges[e].n_holes, sizeof(uint32_t), 1, fp);
		if (g->edges[e].n_holes) {
			g->edges[e].p_holes = malloc(g->edges[e].n_holes * sizeof(uint32_t));
			g->edges[e].l_holes = malloc(g->edges[e].n_holes * sizeof(uint32_t));
			xfread(g->edges[e].p_holes, sizeof(uint32_t), g->edges[e].n_holes, fp);
			xfread(g->edges[e].l_holes, sizeof(uint32_t), g->edges[e].n_holes, fp);
		} else {
			g->edges[e].p_holes = NULL;
			g->edges[e].l_holes = NULL;
		}
	}

	/* load the barcode information */
	if (g->aux_flag & ASM_HAVE_BARCODE) {
		for (e = 0; e < g->n_e; ++e) {
			if (g->edges[e].source == -1)
				continue;
			g->edges[e].barcodes = calloc(3, sizeof(struct barcode_hash_t));
			struct barcode_hash_t *h = g->edges[e].barcodes;
			xfread(&h->size, sizeof(uint32_t), 1, fp);
			xfread(&h->n_item, sizeof(uint32_t), 1, fp);
			h->keys = malloc(h->size * sizeof(uint64_t));
			xfread(h->keys, sizeof(uint64_t), h->size, fp);
			h->cnts = NULL;

			h = g->edges[e].barcodes + 1;
			xfread(&h->size, sizeof(uint32_t), 1, fp);
			xfread(&h->n_item, sizeof(uint32_t), 1, fp);
			h->keys = malloc(h->size * sizeof(uint64_t));
			xfread(h->keys, sizeof(uint64_t), h->size, fp);
			h->cnts = NULL;

			h = g->edges[e].barcodes + 2;
			xfread(&h->size, sizeof(uint32_t), 1, fp);
			xfread(&h->n_item, sizeof(uint32_t), 1, fp);
			h->keys = malloc(h->size * sizeof(uint64_t));
			xfread(h->keys, sizeof(uint64_t), h->size, fp);
			h->cnts = NULL;
		}
	}

	if (g->aux_flag & ASM_HAVE_BARCODE_SCAF) {
        for (e = 0; e < g->n_e; ++e) {
            if (g->edges[e].source == -1)
                continue;
            struct barcode_hash_t *h = &g->edges[e].barcodes_scaf;
            xfread(&h->size, sizeof(uint32_t), 1, fp);
            xfread(&h->n_item, sizeof(uint32_t), 1, fp);
            h->keys = malloc(h->size * sizeof(uint64_t));
            xfread(h->keys, sizeof(uint64_t), h->size, fp);
            h->cnts = NULL;

            h = &g->edges[e].barcodes_scaf2;
            xfread(&h->size, sizeof(uint32_t), 1, fp);
            xfread(&h->n_item, sizeof(uint32_t), 1, fp);
            h->keys = malloc(h->size * sizeof(uint64_t));
            xfread(h->keys, sizeof(uint64_t), h->size, fp);
            h->cnts = NULL;
        }
	}
	fclose(fp);
}

int asm_fasta_edge_convert(struct asm_graph_t *g, gint_t e, kseq_t *seq)
{
	uint32_t *p_holes, *l_holes, *bseq;
	uint32_t n_holes, i, j, k, c, last_c, m_seq;
	p_holes = NULL;
	l_holes = NULL;
	m_seq = 0x100;
	bseq = calloc(m_seq, sizeof(uint32_t));
	last_c = 0;
	n_holes = 0;
	for (i = k = 0; i < seq->seq.l; ++i) {
		c = nt4_table[(int)seq->seq.s[i]];
		if (c >= 4) {
			if (last_c >= 4) {
				++l_holes[j];
			} else {
				if (k == 0) {
					free(bseq);
					return 0;
				}
				p_holes = realloc(p_holes, (n_holes + 1) * sizeof(uint32_t));
				l_holes = realloc(l_holes, (n_holes + 1) * sizeof(uint32_t));
				j = n_holes;
				++n_holes;
				p_holes[j] = k - 1;
				l_holes[j] = 1;
			}
		} else {
			/* append new char */
			if (((k + 15) >> 4) >= m_seq) {
				uint32_t new_m = m_seq << 1;
				bseq = realloc(bseq, new_m * sizeof(uint32_t));
				memset(bseq + m_seq, 0, m_seq * sizeof(uint32_t));
				m_seq = new_m;
			}
			__binseq_set(bseq, k, c);
			++k;
		}
		last_c = c;
	}
	g->edges[e].seq_len = k;
	g->edges[e].seq = realloc(bseq, ((k + 15) >> 4) * sizeof(uint32_t));
	g->edges[e].count = 0;
	g->edges[e].n_holes = n_holes;
	g->edges[e].l_holes = l_holes;
	g->edges[e].p_holes = p_holes;
	return 1;
}

void load_asm_graph_fasta(struct asm_graph_t *g, const char *path, int ksize)
{
	g->ksize = ksize;
	g->aux_flag = 0;
	g->edges = NULL;
	g->nodes = NULL;
	g->n_v = g->n_e = 0;
	gzFile fp = gzopen(path, "r");
	if (!fp)
		__ERROR("Unable to open file [%s] to read", path);
	kseq_t *seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		/* add new edge */
		if ((int)seq->seq.l < ksize)
			continue;
		g->edges = realloc(g->edges, (g->n_e + 2) * sizeof(struct asm_edge_t));
		g->nodes = realloc(g->nodes, (g->n_v + 4) * sizeof(struct asm_node_t));

		if (!asm_fasta_edge_convert(g, g->n_e, seq))
			continue;
		asm_clone_seq_reverse(g->edges + (g->n_e + 1), g->edges + g->n_e);
		// asm_clone_reverse(g->edges + (g->n_e + 1), g->edges + g->n_e);
		g->edges[g->n_e].rc_id = g->n_e + 1;
		g->edges[g->n_e + 1].rc_id = g->n_e;

		g->edges[g->n_e].source = g->n_v;
		g->edges[g->n_e].target = g->n_v + 1;
		g->nodes[g->n_v].adj = malloc(sizeof(gint_t));
		g->nodes[g->n_v].adj[0] = g->n_e;
		g->nodes[g->n_v].deg = 1;
		g->nodes[g->n_v + 1].adj = NULL;
		g->nodes[g->n_v + 1].deg = 0;

		g->edges[g->n_e + 1].source = g->n_v + 2;
		g->edges[g->n_e + 1].target = g->n_v + 3;
		g->nodes[g->n_v + 2].adj = malloc(sizeof(gint_t));
		g->nodes[g->n_v + 2].adj[0] = g->n_e + 1;
		g->nodes[g->n_v + 2].deg = 1;
		g->nodes[g->n_v + 3].adj = NULL;
		g->nodes[g->n_v + 3].deg = 0;

		g->nodes[g->n_v].rc_id = g->n_v + 3;
		g->nodes[g->n_v + 3].rc_id = g->n_v;
		g->nodes[g->n_v + 1].rc_id = g->n_v + 2;
		g->nodes[g->n_v + 2].rc_id = g->n_v + 1;
		g->n_e += 2;
		g->n_v += 4;
	}
	kseq_destroy(seq);
	gzclose(fp);
}

void asm_graph_destroy(struct asm_graph_t *g)
{
	gint_t u, e;
	for (e = 0; e < g->n_e; ++e){
		int u = g->edges[e].source;
		if (u == -1)
			continue;
		asm_clean_edge(g, e);
	}
	free(g->edges);
	g->edges = NULL;
	for (u = 0; u < g->n_v; ++u) {
		free(g->nodes[u].adj);
		g->nodes[u].adj = NULL;
	}
	free(g->nodes);
	g->nodes = NULL;
	g->n_e = g->n_v = 0;
	g->ksize = 0;
	g->bin_size = 0;
}

void asm_clone_graph(struct asm_graph_t *g0, struct asm_graph_t *g1)
{
	save_asm_graph(g0, "tmp_graph.bin");
	load_asm_graph(g1, "tmp_graph.bin");
	return;
	g1->ksize = g0->ksize;
	g1->bin_size = g0->bin_size;
	g1->aux_flag = g0->aux_flag;
	g1->n_v = g0->n_v;
	g1->n_e = g0->n_e;
	if (g0->candidates != NULL){
		g1->candidates = kh_init(pair_contig_count);
		for (khiter_t it = kh_begin(g0->candidates); it != kh_end(g0->candidates);
				++it){
			if (!kh_exist(g0->candidates, it))
				continue;
			struct pair_contig_t key = kh_key(g0->candidates, it);
			struct contig_count_t val = kh_val(g0->candidates, it);
			int ret;
			khiter_t it2 = kh_put(pair_contig_count, g1->candidates, key, &ret);
			kh_val(g1->candidates, it2) = val;
		}
	} else {
		g1->candidates = NULL;
	}
	g1->nodes = (struct asm_node_t *) calloc(g0->n_v, sizeof(struct asm_node_t));
	for (int i = 0; i < g1->n_v; ++i){
		g1->nodes[i].rc_id = g0->nodes[i].rc_id;
		g1->nodes[i].deg = g0->nodes[i].deg;
		g1->nodes[i].adj = (gint_t *) calloc(g1->nodes[i].deg, sizeof(gint_t));
		memcpy(g1->nodes[i].adj, g0->nodes[i].adj, sizeof(gint_t) * g1->nodes[i].deg);
	}

	g1->edges = (struct asm_edge_t *) calloc(g0->n_e, sizeof(struct asm_edge_t));
	for (int i = 0; i < g1->n_e; ++i){
		g1->edges[i].count = g0->edges[i].count;
		g1->edges[i].seq_len = g0->edges[i].seq_len;
		g1->edges[i].n_holes = g0->edges[i].n_holes;
		g1->edges[i].source = g0->edges[i].source;
		g1->edges[i].target = g0->edges[i].target;
		g1->edges[i].rc_id = g0->edges[i].rc_id;
		pthread_mutex_init(&g1->edges[i].lock, NULL);
		__VERBOSE("seq len %d\n", g1->edges[i].seq_len);
		__VERBOSE("%d\n", g0->edges[i].seq == NULL);
		g1->edges[i].seq = (uint32_t *) calloc((g1->edges[i].seq_len + 3) / 4,
				sizeof(uint32_t));
		for (int j = 0; j < g0->edges[i].seq_len; ++j)
			__VERBOSE("%d", __binseq_get(g0->edges[i].seq, j));
		__VERBOSE("\n");
		__VERBOSE("copy mem\n");
		memcpy(g1->edges[i].seq, g0->edges[i].seq,
				sizeof(uint32_t) * ((g1->edges[i].seq_len + 3) / 4));
		__VERBOSE("done seq\n");
		if (g1->edges[i].n_holes > 0){
			g1->edges[i].p_holes = (uint32_t *) calloc(g1->edges[i].n_holes,
					sizeof(uint32_t));
			memcpy(g1->edges[i].p_holes, g0->edges[i].p_holes,
					sizeof(uint32_t) * g1->edges[i].n_holes);
			g1->edges[i].l_holes = (uint32_t *) calloc(g1->edges[i].n_holes,
					sizeof(uint32_t));
			memcpy(g1->edges[i].l_holes, g0->edges[i].l_holes,
					sizeof(uint32_t) * g1->edges[i].n_holes);
		}
	}
	__VERBOSE("DONE\n");
}
