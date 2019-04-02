#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "fastq_producer.h"
#include "io_utils.h"
#include "k31hash.h"
#include "k63hash.h"
#include "k31_count.h"
#include "k63_count.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

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


void k63_build0(struct opt_count_t *opt, int ksize, struct asm_graph_t *g0)
{
	struct k63hash_t *kmer_hash;
	kmer_hash = calloc(1, sizeof(struct k63hash_t));
	__VERBOSE("Estimating kmer\n");
	build_k63_table_lazy(opt, kmer_hash, ksize);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Estimating kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nBuilding assembly graph\n");
	build_asm_graph_from_k63(opt, ksize, kmer_hash, g0);
	__VERBOSE_LOG("kmer_%d_graph_#0", "Number of nodes: %lld\n", ksize,
							(long long)g0->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#0", "Number of edges: %lld\n", ksize,
							(long long)g0->n_e);
}

void k31_build0(struct opt_count_t *opt, int ksize, struct asm_graph_t *g0)
{
	struct k31hash_t *kmer_hash;
	kmer_hash = calloc(1, sizeof(struct k31hash_t));
	__VERBOSE("Estimating kmer\n");
	build_k31_table_lazy(opt, kmer_hash, ksize);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Estimating kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nBuilding assembly graph\n");
	build_asm_graph_from_k31(opt, ksize, kmer_hash, g0);
	__VERBOSE_LOG("kmer_%d_graph_#0", "Number of nodes: %lld\n", ksize,
							(long long)g0->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#0", "Number of edges: %lld\n", ksize,
							(long long)g0->n_e);
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
		uint32_t len = get_edge_len(g->edges + e);
		if (len > max_len) {
			max_len = len;
			ret = e;
		}
	}
	return ret;
}

float get_genome_coverage(struct asm_graph_t *g)
{
	/* Using the coverage of the longest contigs */
	gint_t e;
	float ret_cov = 0.0;
	uint32_t max_len = 0;
	for (e = 0; e < g->n_e; ++e) {
		if (g->edges[e].seq_len > max_len) {
			max_len = g->edges[e].seq_len;
			ret_cov = __get_edge_cov(g->edges + e, g->ksize);
		}
	}
	return ret_cov;
}

gint_t dump_edge_seq(char **seq, uint32_t *m_seq, struct asm_edge_t *e)
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

void asm_clone_edge2(struct asm_graph_t *g, gint_t dst, gint_t src)
{
	asm_clone_edge(g->edges + dst, g->edges + src);
	g->edges[dst].source = g->edges[src].source;
	g->edges[dst].target = g->edges[src].target;
	/* clone the bucks */
	gint_t slen, nbin;
	slen = get_edge_len(g->edges + dst);
	nbin = (slen + g->bin_size - 1) / g->bin_size;
	g->edges[dst].bucks = malloc(nbin * sizeof(struct barcode_hash_t));
	gint_t i;
	for (i = 0; i < nbin; ++i)
		barcode_hash_clone(g->edges[dst].bucks + i,
						g->edges[src].bucks + i);
}

void asm_clone_edge(struct asm_edge_t *dst, struct asm_edge_t *src)
{
	dst->count = src->count;
	dst->seq_len = src->seq_len;
	dst->seq = calloc((dst->seq_len + 15) >> 4, sizeof(uint32_t));
	memcpy(dst->seq, src->seq,
		((dst->seq_len + 15) >> 4) * sizeof(uint32_t));
	dst->n_holes = src->n_holes;
	dst->p_holes = malloc(dst->n_holes * sizeof(uint32_t));
	memcpy(dst->p_holes, src->p_holes, dst->n_holes * sizeof(uint32_t));
	dst->l_holes = malloc(dst->n_holes * sizeof(uint32_t));
	memcpy(dst->l_holes, src->l_holes, dst->n_holes * sizeof(uint32_t));
}

void asm_clone_reverse(struct asm_edge_t *dst, struct asm_edge_t *src)
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

void asm_append_seq_with_gap2(struct asm_graph_t *g, gint_t e1, gint_t e2,
							uint32_t gap_size)
{
	gint_t slen, slen1, slen2, nbin, nbin1, nbin2;
	/* append the bucket */
	slen1 = get_edge_len(g->edges + e1);
	slen2 = get_edge_len(g->edges + e2);
	nbin1 = (slen1 + g->bin_size - 1) / g->bin_size;
	nbin2 = (slen2 + g->bin_size - 1) / g->bin_size;
	asm_append_seq_with_gap(g->edges + e1, g->edges + e2, gap_size);
	if (g->edges[e1].bucks == NULL)
		return;
	slen = get_edge_len(g->edges + e1);
	nbin = (slen + g->bin_size - 1) / g->bin_size;
	g->edges[e1].bucks = realloc(g->edges[e1].bucks,
					nbin * sizeof(struct barcode_hash_t));
	gint_t i;
	if (nbin == nbin1 + nbin2) {
		/* just concat 2 bucks */
		for (i = 0; i < nbin2; ++i)
			barcode_hash_clone(g->edges[e1].bucks + nbin1 + i,
				g->edges[e2].bucks + i);
	} else if (nbin + 1 == nbin1 + nbin2) {
		/* merge the first bucket of e2 to last bucket of e1 */
		barcode_hash_merge(g->edges[e1].bucks + (nbin1 - 1),
							g->edges[e2].bucks);
		for (i = 1; i < nbin2; ++i)
			barcode_hash_clone(g->edges[e1].bucks + nbin1 + i - 1,
				g->edges[e2].bucks + i);
	} else if (nbin + 2 == nbin1 + nbin2) {
		barcode_hash_merge(g->edges[e1].bucks + (nbin1 - 1),
							g->edges[e2].bucks);
		for (i = 1; i + 1 < nbin2; ++i)
			barcode_hash_clone(g->edges[e1].bucks + nbin1 + i - 1,
				g->edges[e2].bucks + i);
	} else if (nbin > nbin1 + nbin2) {
		/* append some empty bin in the middle */
		for (i = 0; i < nbin - nbin1 - nbin2; ++i)
			barcode_hash_init(g->edges[e1].bucks + nbin1 + i, 4);
		for (i = 0; i < nbin2; ++i)
			barcode_hash_clone(g->edges[e1].bucks + (nbin - nbin2 + i),
				g->edges[e2].bucks + i);
	} else {
		assert(0 && "wrong barcode list merge");
	}
}

void asm_append_edge_seq2(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	gint_t slen, slen1, slen2, nbin, nbin1, nbin2;
	/* append the bucket */
	slen1 = get_edge_len(g->edges + e1);
	slen2 = get_edge_len(g->edges + e2);
	nbin1 = (slen1 + g->bin_size - 1) / g->bin_size;
	nbin2 = (slen2 + g->bin_size - 1) / g->bin_size;
	asm_append_edge_seq(g->edges + e1, g->edges + e2, g->ksize);
	if (g->edges[e1].bucks == NULL)
		return;
	slen = get_edge_len(g->edges + e1);
	nbin = (slen + g->bin_size - 1) / g->bin_size;
	g->edges[e1].bucks = realloc(g->edges[e1].bucks,
					nbin * sizeof(struct barcode_hash_t));
	gint_t i;
	if (nbin == nbin1 + nbin2) {
		/* just concat 2 bucks */
		for (i = 0; i < nbin2; ++i)
			barcode_hash_clone(g->edges[e1].bucks + nbin1 + i,
				g->edges[e2].bucks + i);
		// memcpy(g->edges[e1].bucks + nbin1, g->edges[e2].bucks,
		// 			nbin2 * sizeof(struct barcode_hash_t));
	} else if (nbin + 1 == nbin1 + nbin2) {
		/* merge the first bucket of e2 to last bucket of e1 */
		barcode_hash_merge(g->edges[e1].bucks + (nbin1 - 1),
							g->edges[e2].bucks);
		for (i = 1; i < nbin2; ++i)
			barcode_hash_clone(g->edges[e1].bucks + nbin1 + i - 1,
				g->edges[e2].bucks + i);
		// memcpy(g->edges[e1].bucks + nbin1, g->edges[e2].bucks + 1,
		// 		(nbin2 - 1) * sizeof(struct barcode_hash_t));
	} else if (nbin + 2 == nbin1 + nbin2) {
		barcode_hash_merge(g->edges[e1].bucks + (nbin1 - 1),
							g->edges[e2].bucks);
		for (i = 1; i + 1 < nbin2; ++i)
			barcode_hash_clone(g->edges[e1].bucks + nbin1 + i - 1,
				g->edges[e2].bucks + i);
		// memcpy(g->edges[e1].bucks + nbin1, g->edges[e2].bucks + 1,
		// 		(nbin2 - 2) * sizeof(struct barcode_hash_t));
	} else {
		assert(0 && "wrong barcode list merge");
	}
}

void asm_join_edge_with_gap(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
		gint_t e2, gint_t e_rc2, uint32_t gap_size, uint64_t gap_count)
{
	if (!is_hole_rc(g->edges + e1, g->edges + e_rc1))
		__ERROR("Error from start e1");
	if (!is_hole_rc(g->edges + e2, g->edges + e_rc2))
		__ERROR("Error from start e2");
	uint32_t j;
	fprintf(stderr, "e = %ld; n_holes = %u\n",
		e1, g->edges[e1].n_holes);
	if (g->edges[e1].n_holes == 0) {
		assert(g->edges[e1].p_holes == NULL);
		assert(g->edges[e1].l_holes == NULL);
	}
	for (j = 0; j < g->edges[e1].n_holes; ++j)
		fprintf(stderr, "p=%u; l=%u\n", g->edges[e1].p_holes[j],
			g->edges[e1].l_holes[j]);
	if (g->edges[e_rc2].n_holes == 0) {
		assert(g->edges[e_rc2].p_holes == NULL);
		assert(g->edges[e_rc2].l_holes == NULL);
	}
	fprintf(stderr, "e = %ld; n_holes = %u\n",
		e_rc2, g->edges[e_rc2].n_holes);
	for (j = 0; j < g->edges[e_rc2].n_holes; ++j)
		fprintf(stderr, "p=%u; l=%u\n", g->edges[e_rc2].p_holes[j],
			g->edges[e_rc2].l_holes[j]);

	asm_append_seq_with_gap2(g, e1, e2, gap_size);
	g->edges[e1].target = g->edges[e2].target;
	g->edges[e1].count += g->edges[e2].count + gap_count;

	asm_append_seq_with_gap2(g, e_rc2, e_rc1, gap_size);
	g->edges[e_rc2].target = g->edges[e_rc1].target;
	g->edges[e_rc2].count += g->edges[e_rc1].count + gap_count;

	g->edges[e1].rc_id = e_rc2;
	g->edges[e_rc2].rc_id = e1;
	if (!is_hole_rc(g->edges + e1, g->edges + e_rc2)) {
		uint32_t j;
		fprintf(stderr, "e = %ld; n_holes = %u; seq_len = %u\n",
			e1, g->edges[e1].n_holes, g->edges[e1].seq_len);
		for (j = 0; j < g->edges[e1].n_holes; ++j)
			fprintf(stderr, "p=%u; l=%u\n", g->edges[e1].p_holes[j],
				g->edges[e1].l_holes[j]);
		fprintf(stderr, "e = %ld; n_holes = %u; seq_len = %u\n",
			e_rc2, g->edges[e_rc2].n_holes, g->edges[e_rc2].seq_len);
		for (j = 0; j < g->edges[e_rc2].n_holes; ++j)
			fprintf(stderr, "p=%u; l=%u\n", g->edges[e_rc2].p_holes[j],
				g->edges[e_rc2].l_holes[j]);
		__ERROR("Join edge with gap failed %ld_%ld %ld_%ld\n",
			e1, e_rc1, e2, e_rc2);
	}

	asm_remove_edge(g, e2);
	asm_remove_edge(g, e_rc1);
}

void asm_join_edge3(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
	gint_t e2, gint_t e_rc2, gint_t e3, gint_t e_rc3, uint64_t added_count)
{
	asm_append_edge_seq2(g, e1, e2);
	asm_append_edge_seq2(g, e1, e3);
	g->edges[e1].target = g->edges[e3].target;
	g->edges[e1].count += g->edges[e3].count + added_count;

	asm_append_edge_seq2(g, e_rc3, e_rc2);
	asm_append_edge_seq2(g, e_rc3, e_rc1);
	g->edges[e_rc3].target = g->edges[e_rc1].target;
	g->edges[e_rc3].count += g->edges[e_rc1].count + added_count;
	
	g->edges[e1].rc_id = e_rc3;
	g->edges[e_rc3].rc_id = e1;

	asm_remove_edge(g, e3);
	asm_remove_edge(g, e_rc1);
}

void asm_append_seq_with_gap(struct asm_edge_t *dst, struct asm_edge_t *src,
							uint32_t gap_size)
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
	memcpy(dst->l_holes + dst->n_holes + 1, src->l_holes, src->n_holes * sizeof(uint32_t));
	dst->n_holes = n_holes;
	dst->seq_len = seq_len;
}

void asm_append_edge_seq(struct asm_edge_t *dst, struct asm_edge_t *src,
							uint32_t overlap)
{
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

	uint32_t i, k;
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

void asm_append_edge(struct asm_edge_t *dst, struct asm_edge_t *src,
							uint32_t overlap)
{
	if (dst->target != src->source)
		__VERBOSE_INFO("WARNING", "Append edge not consecutive\n");
	asm_append_edge_seq(dst, src, overlap);
	dst->count += src->count;
	dst->target = src->target;
}

void asm_append_edge2(struct asm_graph_t *g, gint_t dst, gint_t src)
{
	if (g->edges[dst].target != g->edges[src].source)
		__VERBOSE_INFO("WARNING", "Append edge not consecutive\n");
	asm_append_edge_seq2(g, dst, src);
	g->edges[dst].count += g->edges[src].count;
	g->edges[dst].target = g->edges[src].target;
}
void asm_clean_edge_seq(struct asm_edge_t *e)
{
	free(e->seq);
	free(e->l_holes);
	free(e->p_holes);
	e->seq = NULL;
	e->l_holes = NULL;
	e->p_holes = NULL;
}

void asm_remove_edge(struct asm_graph_t *g, gint_t e)
{
	assert(e < g->n_e);
	asm_clean_edge_seq(g->edges + e);
	gint_t u = g->edges[e].source;
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
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		gint_t cc_id = id_edge[e];
		if (cc_size[cc_id] < MIN_CONNECT_SIZE ||
			g->edges[e].seq_len < MIN_CONTIG_LEN)
			continue;
		if (e == 1300586)
		// if (e == 74)
			print_debug(g->edges + e);
		gint_t len = dump_edge_seq(&seq, &seq_len, g->edges + e);
		fprintf(fp, ">SEQ_%lld_length_%lld_count_%llu\n", (long long)e,
			(long long)len, (long long unsigned)g->edges[e].count);
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
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		gint_t cc_id = id_edge[e];
		if (cc_size[cc_id] < 250)
			continue;
		dump_edge_seq(&seq, &seq_len, g->edges + e);
		uint64_t fake_count = get_bandage_count(g->edges + e, g->ksize);
		/* print fake count for correct coverage display on Bandage */
		fprintf(fp, "S\t%lld_%lld\t%s\tKC:i:%llu\n", (long long)e,
			(long long)e_rc, seq, (long long unsigned)fake_count);
	}
	for (e = 0; e < g->n_e; ++e) {
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
				dump_edge_seq(&seq, &lseq, g->edges + e);
				fprintf(stderr, "seq    = %s\n", seq);
				dump_edge_seq(&seq, &lseq, g->edges + e_rc);
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
		dump_edge_seq(&seq, &lseq, g->edges + e);
		__VERBOSE("e = %ld: %s\n", e, seq);
	}
	free(seq);
}

void test_asm_graph(struct asm_graph_t *g)
{
	gint_t le_idx = get_longest_edge(g);
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
			dump_edge_seq(&seq, &lseq, g->edges + e);
			__VERBOSE("%s\n", seq);
			__VERBOSE("edge [%ld](%ld->%ld); rc_id = %ld\n",
				e_rc, g->edges[e_rc].source, g->edges[e_rc].target,
				g->edges[e_rc].rc_id);
			dump_edge_seq(&seq, &lseq, g->edges + e_rc);
			__VERBOSE("%s\n", seq);
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

	for (e = 0; e < g->n_e; ++e) {
		uint32_t len = get_edge_len(g->edges + e);
		double cov = __get_edge_cov(g->edges + e, g->ksize);
		if (len > 5000 && cov < 100.0)
			__VERBOSE("WARNING: Edge %ld has length %u with cov ~ %.6lf\n",
				e, len, cov);
	}
}

void save_graph(struct asm_graph_t *g, FILE *fp)
{
	xfwrite(&g->ksize, sizeof(int), 1, fp);
	xfwrite(&g->n_v, sizeof(gint_t), 1, fp);
	xfwrite(&g->n_e, sizeof(gint_t), 1, fp);
	gint_t u, e;
	for (u = 0; u < g->n_v; ++u) {
		xfwrite(&g->nodes[u].rc_id, sizeof(gint_t), 1, fp);
		xfwrite(&g->nodes[u].deg, sizeof(gint_t), 1, fp);
		if (g->nodes[u].deg)
			xfwrite(g->nodes[u].adj, sizeof(gint_t), g->nodes[u].deg, fp);
	}
	for (e = 0; e < g->n_e; ++e) {
		xfwrite(&g->edges[e].source, sizeof(gint_t), 1, fp);
		xfwrite(&g->edges[e].target, sizeof(gint_t), 1, fp);
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
}

void save_asm_graph_simple(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "wb");
	char *sig = "asmg";
	xfwrite(sig, 4, 1, fp);
	save_graph(g, fp);
	fclose(fp);
}

void save_asm_graph_barcode(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "wb");
	char *sig = "asmb";
	xfwrite(sig, 4, 1, fp);
	save_graph(g, fp);
	xfwrite(&g->bin_size, sizeof(int), 1, fp);
	gint_t e;
	for (e = 0; e < g->n_e; ++e) {
		gint_t n, k, len;
		len = get_edge_len(g->edges + e);
		n = (len + g->bin_size - 1) / g->bin_size;
		for (k = 0; k < n; ++k) {
			struct barcode_hash_t *h = g->edges[e].bucks + k;
			xfwrite(&h->size, sizeof(uint32_t), 1, fp);
			xfwrite(&h->n_item, sizeof(uint32_t), 1, fp);
			xfwrite(h->keys, sizeof(uint64_t), h->size, fp);
			xfwrite(h->cnts, sizeof(uint32_t), h->size, fp);
		}
	}
	fclose(fp);
}

void load_graph(struct asm_graph_t *g, FILE *fp)
{
	xfread(&g->ksize, sizeof(int), 1, fp);
	xfread(&g->n_v, sizeof(gint_t), 1, fp);
	xfread(&g->n_e, sizeof(gint_t), 1, fp);
	g->nodes = malloc(g->n_v * sizeof(struct asm_node_t));
	g->edges = malloc(g->n_e * sizeof(struct asm_edge_t));
	gint_t u, e;
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
}

void load_barcode(struct asm_graph_t *g, FILE *fp)
{
	xfread(&g->bin_size, sizeof(int), 1, fp);
	gint_t e;
	for (e = 0; e < g->n_e; ++e) {
		gint_t n, k;
		gint_t len = get_edge_len(g->edges + e);
		n = (len + g->bin_size - 1) / g->bin_size;
		g->edges[e].bucks = calloc(n, sizeof(struct barcode_hash_t));
		for (k = 0; k < n; ++k) {
			struct barcode_hash_t *h = g->edges[e].bucks + k;
			xfread(&h->size, sizeof(uint32_t), 1, fp);
			xfread(&h->n_item, sizeof(uint32_t), 1, fp);
			h->keys = malloc(h->size * sizeof(uint64_t));
			h->cnts = malloc(h->size * sizeof(uint32_t));
			xfread(h->keys, sizeof(uint64_t), h->size, fp);
			xfread(h->cnts, sizeof(uint32_t), h->size, fp);
		}
	}
}

void load_asm_graph(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "rb");
	char sig[4];
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, "asmg", 4) && strncmp(sig, "asmb", 4))
		__ERROR("Not assembly graph format file");
	load_graph(g, fp);
	if (strncmp(sig, "asmb", 4) == 0)
		load_barcode(g, fp);
	fclose(fp);
}

void load_asm_graph_simple(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "rb");
	char sig[4];
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, "asmg", 4) && strncmp(sig, "asmb", 4))
		__ERROR("Not assembly graph format file");
	load_graph(g, fp);
	fclose(fp);
}

void load_asm_graph_complex(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "rb");
	char sig[4];
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, "asmb", 4))
		__ERROR("Not assembly graph with barcode format file");
	load_graph(g, fp);
	load_barcode(g, fp);
	fclose(fp);
}
