#include <stdlib.h>

#include "fastq_producer.h"
#include "graph_assembly.h"
#include "io_utils.h"
#include "kmer_count.h"
#include "time_utils.h"
#include "utils.h"
#include "verbose.h"

#define BUCK_SORT_SIZE		64

#define __get_revc_num(y, x, l, mask)					       \
(	(x) = (y) << (64 - ((l) << 1)),					       \
	(x) = (((x) & 0xffffffff00000000ull) >> 32) | (((x) & 0x00000000ffffffffull) << 32), \
	(x) = (((x) & 0xffff0000ffff0000ull) >> 16) | (((x) & 0x0000ffff0000ffffull) << 16), \
	(x) = (((x) & 0xff00ff00ff00ff00ull) >>  8) | (((x) & 0x00ff00ff00ff00ffull) <<  8), \
	(x) = (((x) & 0xf0f0f0f0f0f0f0f0ull) >>  4) | (((x) & 0x0f0f0f0f0f0f0f0full) <<  4), \
	(x) = (((x) & 0xccccccccccccccccull) >>  2) | (((x) & 0x3333333333333333ull) <<  2), \
	(x) ^= 0xffffffffffffffffull, (x) &= (mask))

#define __off_bit(a, i) ((a)[(i) >> 5] &= (~((uint32_t)1 << ((i) & 31))))
#define __on_bit(a, i) ((a)[(i) >> 5] |= ((uint32_t)1 << ((i) & 31)))
#define __get_bit(a, i) (1 & ((a)[(i) >> 5] >> ((i) & 31)))

#define __degree(e) (((e) & 1) + (((e) >> 1) & 1) + (((e) >> 2) & 1) + (((e) >> 3) & 1))

/* e must have degree 1 */
#define __only_edge(e) ((((e) >> 1) & 1) * 1 + (((e) >> 2) & 1) * 2 + (((e) >> 3) & 1) * 3)

#define __set_degree(bin, k, deg) ((bin)[(k) >> 4] ^= (((bin)[(k) >> 4] & (3u << (((k) & 15) << 1))) ^ ((uint32_t)(deg) << (((k) & 15) << 1))))

#define __get_degree(bin, k) (((bin)[(k) >> 4] >> (((k) & 15) << 1)) & 3)

struct raw_node_t {
	kmkey_t kmer;
	kmval_t cnt;
	uint8_t adj;
};

struct raw_graph_t {
	kmint_t n;
	struct raw_node_t *nodes;
};

struct edgecount_bundle_t {
	struct dqueue_t *q;
	struct raw_graph_t *g;
	int ksize;
	int64_t *n_reads;
};

struct raw_graph_t *extract_kmer(struct opt_count_t *opt, struct kmhash_t *h);

static kmint_t bin_search_id(struct raw_graph_t *g, kmkey_t x);

static void dump_bin_seq(uint32_t *bin, char *seq, int len);

static uint32_t *init_seq_bin(kmkey_t key, int l);

static void rev_bin_seq(uint32_t *seq, int beg, int len);

static uint32_t *append_bin_seq_rv(uint32_t *dest, int cur_len, uint32_t *src, int beg, int len);

static uint32_t *append_bin_seq_fw(uint32_t *dest, int cur_len, uint32_t *src, int beg, int len);

void sort_kmer(struct raw_graph_t *g);

void test_sort_kmer(struct raw_graph_t *g);

void get_edges(struct opt_count_t *opt, struct raw_graph_t *g);

void get_edge_stat(struct raw_graph_t *g);

struct scrap_graph_t *sketch_graph(struct raw_graph_t *pre_graph, int ksize);

void dump_scrap_graph(struct scrap_graph_t *g, struct opt_count_t *opt);

struct scrap_graph_t *remove_tips_round_1(struct scrap_graph_t *pre_g);

void assembly_process(struct opt_count_t *opt)
{
	init_clock();
	__VERBOSE("Counting kmer\n");
	struct kmhash_t *kmer_hash;
	kmer_hash = count_kmer(opt);

	struct raw_graph_t *pre_graph;
	pre_graph = extract_kmer(opt, kmer_hash);
	kmhash_destroy(kmer_hash);

	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Counting kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nBuilding kmer graph\n");
	sort_kmer(pre_graph);
	// test_sort_kmer(pre_graph);

	get_edges(opt, pre_graph);
	// get_edge_stat(pre_graph);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Building kmer graph time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nGluing non-branching path #1\n");
	struct scrap_graph_t *scratch_graph;
	scratch_graph = sketch_graph(pre_graph, opt->kmer_master);
	// dump_scrap_graph(scratch_graph, opt);
	
	__VERBOSE_LOG("TIMER", "Gluing non-branching path #1 time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nRemoving tips #1 and gluing non-branching path #2\n");
	struct scrap_graph_t *bad_graph;
	bad_graph = remove_tips_round_1(scratch_graph);

	__VERBOSE_LOG("TIMER", "Removing tips #1 and gluing non-branching path #2 time: %.3f\n", sec_from_prev_time());

	__VERBOSE("Printing graph in gfa format\n");
	dump_scrap_graph(bad_graph, opt);

	__VERBOSE_LOG("TIMER", "Total time: %.3f\n", sec_from_initial_time());
}

void test_graph(struct scrap_graph_t *g)
{
	gint_t k, v;
	int deg, c;
	for (k = 0; k < g->n_v; ++k) {
		// forward edges
		deg = __get_degree(g->bin_fdeg, k);
		if (deg == 0)
			deg += (int)(g->fadj[k] != NULL);
		else
			assert(g->fadj[k] != NULL);
		assert(deg >= 0 && deg <= 4);
		for (c = 0; c < deg; ++c) {
			v = g->fadj[k][c];
			if (v > 0)
				--v;
			else
				v = -v - 1;
			assert(v >= 0 && v < g->n_v);
		}

		// reverse edges
		deg = __get_degree(g->bin_rdeg, k);
		if (deg == 0)
			deg += (int)(g->radj[k] != NULL);
		else
			assert(g->radj[k] != NULL);
		assert(deg >= 0 && deg <= 4);
		for (c = 0; c < deg; ++c) {
			v = g->radj[k][c];
			if (v > 0)
				--v;
			else
				v = -v - 1;
			assert(v >= 0 && v < g->n_v);
		}
	}
	fprintf(stderr, "test graph ok\n");
}

void test_graph_2(struct scrap_graph_t *g)
{
	// Expecting graph does not have 1-1 edge
	gint_t k, v;
	int deg, rdeg;
	for (k = 0; k < g->n_v; ++k) {
		deg = __get_degree(g->bin_fdeg, k) + (int)(g->fadj[k] != NULL);
		if (deg == 1) {
			v = g->fadj[k][0];
			if (v > 0)
				rdeg = __get_degree(g->bin_rdeg, v - 1) + (int)(g->radj[v - 1] != NULL);
			else
				rdeg = __get_degree(g->bin_fdeg, -v - 1) + (int)(g->fadj[-v - 1] != NULL);
			if (rdeg == 1 && v != k + 1 && v != -k - 1) {
				__ERROR("1-1 edge exist (%lld) -> (%lld)\n", (long long)k + 1, (long long)v);
			}
		}

		deg = __get_degree(g->bin_rdeg, k) + (int)(g->radj[k] != NULL);
		if (deg == 1) {
			v = g->radj[k][0];
			if (v > 0)
				rdeg = __get_degree(g->bin_rdeg, v - 1) + (int)(g->radj[v - 1] != NULL);
			else
				rdeg = __get_degree(g->bin_fdeg, -v - 1) + (int)(g->fadj[-v - 1] != NULL);
			if (rdeg == 1 && v != k + 1 && v != -k - 1) {
				__ERROR("1-1 edge exist (%lld) -> (%lld)\n", (long long)(-k - 1), (long long)v);
			}
		}
	}
	__VERBOSE("test graph 2 ok\n");
}

uint32_t *remove_tips(struct scrap_graph_t *g)
{
	gint_t k, uk, v;
	gint_t *uk_adj;
	int deg, old_deg, new_deg, fdeg, rdeg, uk_deg, c;
	float cov, max_cov;

	uint32_t *removed;
	removed = calloc((g->n_v + 31) / 32, sizeof(uint32_t));

	// mark remove nodes
	for (k = 0; k < g->n_v; ++k) {
		fdeg = __get_degree(g->bin_fdeg, k) + (int)(g->fadj[k] != NULL);
		rdeg = __get_degree(g->bin_rdeg, k) + (int)(g->radj[k] != NULL);
		deg = fdeg + rdeg;
		if (deg <= 1) {
			if (deg == 0) {
				__on_bit(removed, k);
				continue;
			}
			if (fdeg == 1)
				uk = g->fadj[k][0];
			else
				uk = g->radj[k][0];
			if (uk > 0) {
				uk_adj = g->radj[uk - 1];
				uk_deg = __get_degree(g->bin_rdeg, uk - 1) + (int)(g->radj[uk - 1] != NULL);
			} else {
				uk_adj = g->fadj[-uk - 1];
				uk_deg = __get_degree(g->bin_fdeg, -uk - 1) + (int)(g->fadj[-uk - 1] != NULL);
			}
			if (uk_deg > 1) {
				max_cov = 0.0;
				for (c = 0; c < uk_deg; ++c) {
					v = uk_adj[c];
					if (v > 0)
						cov = 1.0 * g->kmer_count[v - 1] / g->seq_len[v - 1];
					else
						cov = 1.0 * g->kmer_count[-v - 1] / g->seq_len[-v - 1];
					if (cov > max_cov)
						max_cov = cov;
				}
				cov = 1.0 * g->kmer_count[k] / g->seq_len[k];
				if (cov / max_cov < 0.1)
					__on_bit(removed, k);
			}
		}
	}

	// remove edge
	for (k = 0; k < g->n_v; ++k) {
		if (__get_bit(removed, k))
			continue;
		// forward edge
		old_deg = __get_degree(g->bin_fdeg, k) + (int)(g->fadj[k] != NULL);
		new_deg = 0;
		if (old_deg) {
			for (c = 0; c < old_deg; ++c) {
				uk = g->fadj[k][c];
				if (uk > 0)
					--uk;
				else
					uk = -uk - 1;
				if (__get_bit(removed, uk) == 0)
					g->fadj[k][new_deg++] = g->fadj[k][c];
			}
			if (new_deg != old_deg) {
				if (new_deg == 0) {
					__set_degree(g->bin_fdeg, k, 0);
					free(g->fadj[k]);
					g->fadj[k] = NULL;
				} else {
					__set_degree(g->bin_fdeg, k, new_deg - 1);
					g->fadj[k] = realloc(g->fadj[k], new_deg * sizeof(gint_t));
				}
			}
		}
		// reverse edge
		old_deg = __get_degree(g->bin_rdeg, k) + (int)(g->radj[k] != NULL);
		new_deg = 0;
		if (old_deg) {
			for (c = 0; c < old_deg; ++c) {
				uk = g->radj[k][c];
				if (uk > 0)
					--uk;
				else
					uk = -uk - 1;
				if (__get_bit(removed, uk) == 0)
					g->radj[k][new_deg++] = g->radj[k][c];
			}
			if (new_deg != old_deg) {
				if (new_deg == 0) {
					__set_degree(g->bin_rdeg, k, 0);
					free(g->radj[k]);
					g->radj[k] = NULL;
				} else {
					__set_degree(g->bin_rdeg, k, new_deg - 1);
					g->radj[k] = realloc(g->radj[k], new_deg * sizeof(gint_t));
				}
			}
		}
	}

	return removed;
}

struct scrap_graph_t *condense_graph(struct scrap_graph_t *g, uint32_t *removed)
{
	struct scrap_graph_t *ret_g;
	ret_g = calloc(1, sizeof(struct scrap_graph_t));

	gint_t *new_id, *node_beg, *node_end, *v_radj, *adj, *new_adj;
	gint_t n_v, k, m_v, u, v, id_u, id_v, v_node;

	uint32_t *node_seq;
	uint32_t deg, new_deg, v_rdeg, vis, j;

	int seq_len, old_seq_len, ksize;
	ksize = g->ksize;

	new_id = malloc(g->n_v * sizeof(gint_t));

	m_v = 0x1000;
	node_beg = malloc(m_v * sizeof(gint_t));
	node_end = malloc(m_v * sizeof(gint_t));
	ret_g->node_seq = malloc(m_v * sizeof(uint32_t *));
	ret_g->seq_len = malloc(m_v * sizeof(int));
	ret_g->kmer_count = malloc(m_v * sizeof(gint_t));

	// condense path
	n_v = 0;
	for (k = 0; k < g->n_v; ++k) {
		if (__get_bit(removed, k))
			continue;

		if (n_v + 1 > m_v) {
			m_v <<= 1;
			node_beg = realloc(node_beg, m_v * sizeof(gint_t));
			node_end = realloc(node_end, m_v * sizeof(gint_t));
			ret_g->node_seq = realloc(ret_g->node_seq, m_v * sizeof(uint32_t *));
			ret_g->seq_len = realloc(ret_g->seq_len, m_v * sizeof(int));
			ret_g->kmer_count = realloc(ret_g->kmer_count, m_v * sizeof(gint_t));
		}

		node_beg[n_v] = node_end[n_v] = k + 1;
		u = k;
		adj = g->fadj[k];
		deg = __get_degree(g->bin_fdeg, k) + (int)(adj != NULL);
		__on_bit(removed, k);
		new_id[k] = n_v;
		seq_len = g->seq_len[u];
		node_seq = append_bin_seq_fw(NULL, 0, g->node_seq[u], 0, g->seq_len[u]);
		ret_g->kmer_count[n_v] = g->kmer_count[k];

		while (deg == 1) {
			v = adj[0];
			if (v > 0) {
				v_radj = g->radj[v - 1];
				v_rdeg = __get_degree(g->bin_rdeg, v - 1) + (int)(v_radj != NULL);
				vis = __get_bit(removed, v - 1);
			} else {
				v_radj = g->fadj[-v - 1];
				v_rdeg = __get_degree(g->bin_fdeg, -v - 1) + (int)(v_radj != NULL);
				vis = __get_bit(removed, -v - 1);
			}
			assert(v_rdeg > 0);
			if (v_rdeg > 1 || vis == 1)
				break;
			assert(u == v_radj[0] - 1 || u == -v_radj[0] - 1);
			if (v > 0) {
				u = v - 1;
				adj = g->fadj[u];
				deg = __get_degree(g->bin_fdeg, u) + (int)(adj != NULL);
				__on_bit(removed, u);
				node_seq = append_bin_seq_fw(node_seq, seq_len, g->node_seq[u], ksize - 1, g->seq_len[u] - ksize + 1);
				seq_len += g->seq_len[u] - ksize + 1;
				new_id[u] = n_v;
				ret_g->kmer_count[n_v] += g->kmer_count[u];
			} else {
				u = -v - 1;
				adj = g->radj[u];
				deg = __get_degree(g->bin_rdeg, u) + (int)(adj != NULL);
				__on_bit(removed, u);
				node_seq = append_bin_seq_rv(node_seq, seq_len, g->node_seq[u], 0, g->seq_len[u] - ksize + 1);
				seq_len += g->seq_len[u] - ksize + 1;
				new_id[u] = n_v;
				ret_g->kmer_count[n_v] += g->kmer_count[u];
			}
			node_end[n_v] = v;
		}

		old_seq_len = seq_len;
		u = k;
		adj = g->radj[k];
		deg = __get_degree(g->bin_rdeg, k) + (int)(adj != NULL);

		while (deg == 1) {
			v = adj[0];
			if (v > 0) {
				v_radj = g->radj[v - 1];
				v_rdeg = __get_degree(g->bin_rdeg, v - 1) + (int)(v_radj != NULL);
				vis = __get_bit(removed, v - 1);
			} else {
				v_radj = g->fadj[-v - 1];
				v_rdeg = __get_degree(g->bin_fdeg, -v - 1) + (int)(v_radj != NULL);
				vis = __get_bit(removed, -v - 1);
			}
			assert(v_rdeg > 0);
			if (v_rdeg > 1 || vis == 1)
				break;
			assert(u == v_radj[0] - 1 || u == -v_radj[0] - 1);
			if (v > 0) {
				u = v - 1;
				adj = g->fadj[u];
				deg = __get_degree(g->bin_fdeg, u) + (int)(adj != NULL);
				__on_bit(removed, u);
				node_seq = append_bin_seq_fw(node_seq, seq_len, g->node_seq[u], ksize - 1, g->seq_len[u] - ksize + 1);
				seq_len += g->seq_len[u] - ksize + 1;
				new_id[u] = n_v;
				ret_g->kmer_count[n_v] += g->kmer_count[u];
			} else {
				u = -v - 1;
				adj = g->radj[u];
				deg = __get_degree(g->bin_rdeg, u) + (int)(adj != NULL);
				__on_bit(removed, u);
				node_seq = append_bin_seq_rv(node_seq, seq_len, g->node_seq[u], ksize - 1, g->seq_len[u] - ksize + 1);
				seq_len += g->seq_len[u] - ksize + 1;
				new_id[u] = n_v;
				ret_g->kmer_count[n_v] += g->kmer_count[u];
			}
			node_beg[n_v] = v;
		}

		if (node_beg[n_v] != k + 1) {
			node_beg[n_v] = -node_beg[n_v];
			rev_bin_seq(node_seq, 0, old_seq_len);
			rev_bin_seq(node_seq, 0, seq_len);
		}

		ret_g->seq_len[n_v] = seq_len;
		ret_g->node_seq[n_v] = node_seq;

		++n_v;
	}

	node_beg = realloc(node_beg, n_v * sizeof(gint_t));
	node_end = realloc(node_end, n_v * sizeof(gint_t));
	ret_g->seq_len = realloc(ret_g->seq_len, n_v * sizeof(int));
	ret_g->node_seq = realloc(ret_g->node_seq, n_v * sizeof(uint32_t *));
	ret_g->kmer_count = realloc(ret_g->kmer_count, n_v * sizeof(gint_t));
	ret_g->fadj = malloc(n_v * sizeof(gint_t *));
	ret_g->radj = malloc(n_v * sizeof(gint_t *));
	ret_g->bin_fdeg = calloc((n_v + 15) >> 4, sizeof(uint32_t));
	ret_g->bin_rdeg = calloc((n_v + 15) >> 4, sizeof(uint32_t));

	for (k = 0; k < n_v; ++k) {
		// forward edge
		u = node_end[k];
		if (u > 0) {
			id_u = u - 1;
			adj = g->fadj[id_u];
			deg = __get_degree(g->bin_fdeg, id_u) + (int)(adj != NULL);
		} else {
			id_u = -u - 1;
			adj = g->radj[id_u];
			deg = __get_degree(g->bin_rdeg, id_u) + (int)(adj != NULL);
		}
		if (deg) {
			new_adj = ret_g->fadj[k] = malloc(deg * sizeof(gint_t));
			new_deg = 0;
			for (j = 0; j < deg; ++j) {
				v = adj[j];
				assert(v >= -g->n_v && v <= g->n_v);
				if (v > 0) {
					id_v = v - 1;
					v_node = new_id[id_v];
				} else {
					id_v = -v - 1;
					v_node = new_id[id_v];
				}
				if (v == node_beg[v_node]) {
					new_adj[new_deg++] = v_node + 1;
				} else if (v == -node_end[v_node]) {
					new_adj[new_deg++] = -v_node - 1;
				} else {
					fprintf(stderr, "n_v = %llu\n", (long long)g->n_v);
					fprintf(stderr, "k = %lld\nu = %lld\nv = %lld\nv_node = %lld\nbeg[v_node]=%lld\nend[v_node]=%lld\n",
						(long long)k, (long long)u, (long long)v,
						(long long)v_node, (long long)node_beg[v_node],
						(long long)node_end[v_node]);
					assert(0);
				}
			}
			assert(deg == new_deg);
			__set_degree(ret_g->bin_fdeg, k, new_deg - 1);
		} else {
			ret_g->fadj[k] = NULL;
			__set_degree(ret_g->bin_fdeg, k, 0);
		}

		// reverse edge
		u = node_beg[k];
		if (u > 0) {
			id_u = u - 1;
			adj = g->radj[id_u];
			deg = __get_degree(g->bin_rdeg, id_u) + (int)(adj != NULL);
		} else {
			id_u = -u - 1;
			adj = g->fadj[id_u];
			deg = __get_degree(g->bin_fdeg, id_u) + (int)(adj != NULL);
		}
		if (deg) {
			new_adj = ret_g->radj[k] = malloc(deg * sizeof(gint_t));
			new_deg = 0;
			for (j = 0; j < deg; ++j) {
				v = adj[j];
				if (v > 0) {
					id_v = v - 1;
					v_node = new_id[id_v];
				} else {
					id_v = -v - 1;
					v_node = new_id[id_v];
				}
				if (v == node_beg[v_node]) {
					new_adj[new_deg++] = v_node + 1;
				} else if (v == -node_end[v_node]) {
					new_adj[new_deg++] = -v_node - 1;
				} else {
					assert(0);
				}
			}
			assert(deg == new_deg);
			__set_degree(ret_g->bin_rdeg, k, new_deg - 1);
		} else {
			ret_g->radj[k] = NULL;
			__set_degree(ret_g->bin_rdeg, k, 0);
		}
	}

	free(new_id);
	free(node_beg);
	free(node_end);

	ret_g->n_v = n_v;

	return ret_g;
}

struct scrap_graph_t *remove_tips_round_1(struct scrap_graph_t *g)
{
	struct scrap_graph_t *ret_g;
	uint32_t *removed;
	
	test_graph(g);

	__VERBOSE("||-- removing tips\n");
	removed = remove_tips(g);

	test_graph(g);

	__VERBOSE("||-- condesing graph\n");
	ret_g = condense_graph(g, removed);
	test_graph(ret_g);
	test_graph_2(ret_g);

	free(removed);

	return ret_g;
}

static uint32_t *append_bin_seq_fw(uint32_t *dest, int cur_len, uint32_t *src, int beg, int len)
{
	int m, i, k, cur_m;
	uint32_t *new_ptr;
	cur_m = (cur_len + 15) >> 4;
	m = (cur_len + len + 15) >> 4;
	if (m > cur_m) {
		new_ptr = realloc(dest, m * sizeof(uint32_t));
		memset(new_ptr + cur_m, 0, (m - cur_m) * sizeof(uint32_t));
	} else {
		new_ptr = dest;
	}
	for (i = beg; i < beg + len; ++i) {
		k = cur_len + (i - beg);
		new_ptr[k >> 4] |= ((src[i >> 4] >> ((i & 15) << 1)) & (uint32_t)3) << ((k & 15) << 1);
	}
	return new_ptr;
}

static uint32_t *append_bin_seq_rv(uint32_t *dest, int cur_len, uint32_t *src, int beg, int len)
{
	int m, i, k, cur_m;
	uint32_t *new_ptr;
	cur_m = (cur_len + 15) >> 4;
	m = (cur_len + len + 15) >> 4;
	if (m > cur_m) {
		new_ptr = realloc(dest, m * sizeof(uint32_t));
		memset(new_ptr + cur_m, 0, (m - cur_m) * sizeof(uint32_t));
	} else {
		new_ptr = dest;
	}
	for (i = beg; i < beg + len; ++i) {
		k = cur_len + (len - (i - beg) - 1);
		new_ptr[k >> 4] |= (((src[i >> 4] >> ((i & 15) << 1)) & (uint32_t)3) ^ (uint32_t)3) << ((k & 15) << 1);
	}
	return new_ptr;
}

struct scrap_graph_t *sketch_graph(struct raw_graph_t *pre_g, int ksize)
{
	struct scrap_graph_t *g;
	g = calloc(1, sizeof(struct scrap_graph_t));

	g->ksize = ksize;

	gint_t m_v, n_v, node_id, u_node, v_node, k, uk;
	gint_t *kmer_chain_id;
	kmint_t i, n_k;
	kmkey_t node_kmer, node_rkmer, u_kmer, u_rkmer, v_kmer, v_rkmer, kmask;
	kmkey_t *kmer_beg, *kmer_end;
	uint8_t u_adj, v_radj, c, u_deg, deg;
	int lmc, seq_len, old_seq_len;
	struct raw_node_t *nodes;
	nodes = pre_g->nodes;
	n_k = pre_g->n;
	kmask = ((kmkey_t)1 << (ksize << 1)) - 1;
	lmc = (ksize - 1) << 1;

	uint32_t *visited, *node_seq;
	visited = calloc((n_k + 31) / 32, sizeof(uint32_t));

	// g->kmer_chain_id = malloc(n_k * sizeof(gint_t));
	kmer_chain_id = malloc(n_k * sizeof(gint_t));
	n_v = 0;
	m_v = 0x1000;
	g->kmer_count = malloc(m_v * sizeof(gint_t));
	// g->kmer_beg = malloc(m_v * sizeof(kmkey_t));
	// g->kmer_end = malloc(m_v * sizeof(kmkey_t));
	g->seq_len = malloc(m_v * sizeof(int));
	g->node_seq = malloc(m_v * sizeof(uint32_t *));

	kmer_beg = malloc(m_v * sizeof(kmkey_t));
	kmer_end = malloc(m_v * sizeof(kmkey_t));

	for (i = 0; i < n_k; ++i) {
		node_id = (gint_t)i;
		node_kmer = nodes[i].kmer;
		if (__get_bit(visited, node_id))
			continue;
		__get_revc_num(node_kmer, node_rkmer, ksize, kmask);
		assert(node_kmer < node_rkmer);

		if (n_v + 1 > m_v) {
			m_v <<= 1;
			g->kmer_count = realloc(g->kmer_count, m_v * sizeof(gint_t));
			// g->kmer_beg = realloc(g->kmer_beg, m_v * sizeof(kmkey_t));
			// g->kmer_end = realloc(g->kmer_end, m_v * sizeof(kmkey_t));
			g->seq_len = realloc(g->seq_len, m_v * sizeof(int));
			g->node_seq = realloc(g->node_seq, m_v * sizeof(uint32_t *));
			kmer_beg = realloc(kmer_beg, m_v * sizeof(kmkey_t));
			kmer_end = realloc(kmer_end, m_v * sizeof(kmkey_t));
		}

		// chain = malloc(sizeof(kmkey_t));
		// chain[0] = node_kmer;
		// lchain = 1;
		g->kmer_count[n_v] = nodes[node_id].cnt;
		__on_bit(visited, node_id);
		// g->kmer_beg[n_v] = g->kmer_end[n_v] = node_kmer;
		kmer_beg[n_v] = kmer_end[n_v] = node_kmer;
		kmer_chain_id[node_id] = n_v;
		node_seq = init_seq_bin(node_kmer, ksize);
		seq_len = ksize;

		// forward
		u_node = node_id;
		u_kmer = node_kmer;
		u_rkmer = node_rkmer;
		u_adj = nodes[u_node].adj & (uint8_t)0xf;

		while (__degree(u_adj) == 1) {
			c = __only_edge(u_adj);
			v_kmer = ((u_kmer << 2) & kmask) | c;
			v_rkmer = (u_rkmer >> 2) | ((kmkey_t)(c ^ 3) << lmc);
			if (v_kmer < v_rkmer) {
				v_node = (gint_t)bin_search_id(pre_g, v_kmer);
				v_radj = nodes[v_node].adj >> 4;
			} else {
				v_node = (gint_t)bin_search_id(pre_g, v_rkmer);
				v_radj = nodes[v_node].adj & (uint8_t)0xf;
			}
			// Check if node v is on another chain
			if (__get_bit(visited, v_node))
				break;
			// Check if 1-1 edge
			if (__degree(v_radj) != 1)
				break;
			assert(u_rkmer == ((((v_rkmer << 2) & kmask) | __only_edge(v_radj))));

			// chain = realloc(chain, lchain + 1);
			// chain[lchain++] = v_kmer;
			// g->kmer_end[n_v] = v_kmer;
			kmer_end[n_v] = v_kmer;
			g->kmer_count[n_v] += nodes[v_node].cnt;
			kmer_chain_id[v_node] = n_v;
			if ((seq_len & 15) == 0) {
				node_seq = realloc(node_seq, ((seq_len >> 4) + 1) * sizeof(uint32_t));
				node_seq[seq_len >> 4] = 0;
			}
			node_seq[seq_len >> 4] |= ((uint32_t)(c) << ((seq_len & 15) << 1));
			++seq_len;
			// ++g->seq_len[n_v];
			__on_bit(visited, v_node);
			u_node = v_node;
			u_kmer = v_kmer;
			u_rkmer = v_rkmer;
			u_adj = nodes[u_node].adj >> (4 * (u_kmer > u_rkmer)) & (uint8_t)0xf;
		}

		// old_lchain = lchain;
		old_seq_len = seq_len;
		u_node = node_id;
		u_kmer = node_rkmer;
		u_rkmer = node_kmer;
		u_adj = nodes[u_node].adj >> 4;

		// reverse
		while (__degree(u_adj) == 1) {
			c = __only_edge(u_adj);
			v_kmer = ((u_kmer << 2) & kmask) | c;
			v_rkmer = (u_rkmer >> 2) | ((kmkey_t)(c ^ 3) << lmc);
			if (v_kmer < v_rkmer) {
				v_node = (gint_t)bin_search_id(pre_g, v_kmer);
				v_radj = nodes[v_node].adj >> 4;
			} else {
				v_node = (gint_t)bin_search_id(pre_g, v_rkmer);
				v_radj = nodes[v_node].adj & (uint8_t)0xf;
			}
			// Check if node v is on another chain
			if (__get_bit(visited, v_node))
				break;
			// Check if 1-1 edge
			if (__degree(v_radj) != 1)
				break;
			assert(u_rkmer == ((((v_rkmer << 2) & kmask) | __only_edge(v_radj))));

			// chain = realloc(chain, lchain + 1);
			// chain[lchain++] = v_kmer;
			// g->kmer_beg[n_v] = v_kmer;
			kmer_beg[n_v] = v_kmer;
			g->kmer_count[n_v] += nodes[v_node].cnt;
			kmer_chain_id[v_node] = n_v;
			if ((seq_len & 15) == 0) {
				node_seq = realloc(node_seq, ((seq_len >> 4) + 1) * sizeof(uint32_t));
				node_seq[seq_len >> 4] = 0;
			}
			node_seq[seq_len >> 4] |= ((uint32_t)(c) << ((seq_len & 15) << 1));
			++seq_len;
			// ++g->seq_len[n_v];
			__on_bit(visited, v_node);
			u_node = v_node;
			u_kmer = v_kmer;
			u_rkmer = v_rkmer;
			u_adj = nodes[u_node].adj >> (4 * (u_kmer > u_rkmer)) & (uint8_t)0xf;
		}

		// correct the chain
		if (kmer_beg[n_v] != node_kmer) {
			// __get_revc_num(g->kmer_beg[n_v], g->kmer_beg[n_v], ksize, kmask);
			__get_revc_num(kmer_beg[n_v], kmer_beg[n_v], ksize, kmask);
			rev_bin_seq(node_seq, 0, old_seq_len);
			rev_bin_seq(node_seq, 0, seq_len);
		}
		g->seq_len[n_v] = seq_len;
		g->node_seq[n_v] = node_seq;
		++n_v;
	}

	g->kmer_count = realloc(g->kmer_count, n_v * sizeof(gint_t));
	// g->kmer_beg = realloc(g->kmer_beg, n_v * sizeof(kmkey_t));
	// g->kmer_end = realloc(g->kmer_end, n_v * sizeof(kmkey_t));
	kmer_beg = realloc(kmer_beg, n_v * sizeof(kmkey_t));
	kmer_end = realloc(kmer_end, n_v * sizeof(kmkey_t));
	g->seq_len = realloc(g->seq_len, n_v * sizeof(int));
	g->node_seq = realloc(g->node_seq, n_v * sizeof(uint32_t *));

	g->fadj = malloc(n_v * sizeof(gint_t *));
	g->radj = malloc(n_v * sizeof(gint_t *));
	g->bin_fdeg = calloc((n_v + 15) / 16, sizeof(uint32_t));
	g->bin_rdeg = calloc((n_v + 15) / 16, sizeof(uint32_t));

	// count edge
	for (k = 0; k < n_v; ++k) {
		// forward
		// u_kmer = g->kmer_end[k];
		u_kmer = kmer_end[k];
		__get_revc_num(u_kmer, u_rkmer, ksize, kmask);
		if (u_kmer < u_rkmer) {
			u_node = (gint_t)bin_search_id(pre_g, u_kmer);
			u_adj = nodes[u_node].adj & (uint8_t)0xf;
		} else {
			u_node = (gint_t)bin_search_id(pre_g, u_rkmer);
			u_adj = nodes[u_node].adj >> 4;
		}
		u_deg = __degree(u_adj);
		if (u_deg) {
			deg = 0;
			g->fadj[k] = malloc(u_deg * sizeof(gint_t));
			for (c = 0; c < 4; ++c) {
				if (!((u_adj >> c) & 1))
					continue;
				v_kmer = ((u_kmer << 2) & kmask) | c;
				v_rkmer = (u_rkmer >> 2) | ((kmkey_t)(c ^ 3) << lmc);
				if (v_kmer < v_rkmer)
					v_node = (gint_t)bin_search_id(pre_g, v_kmer);
				else
					v_node = (gint_t)bin_search_id(pre_g, v_rkmer);
				uk = kmer_chain_id[v_node];
				assert(uk >= 0 && uk < n_v);
				if (v_kmer == kmer_beg[uk])
					g->fadj[k][deg++] = uk + 1;
				else if (v_rkmer == kmer_end[uk])
					g->fadj[k][deg++] = -(uk + 1);
				else
					assert(0);
			}
			if (deg != u_deg)
				fprintf(stderr, "[DEBUG] deg = %d; u_deg = %d\n", (int)deg, (int)u_deg);
			assert(deg == u_deg);
			__set_degree(g->bin_fdeg, k, u_deg - 1);
			assert(__get_degree(g->bin_fdeg, k) == u_deg - 1);
		} else {
			g->fadj[k] = NULL;
			__set_degree(g->bin_fdeg, k, 0);
			assert(__get_degree(g->bin_fdeg, k) == 0);
		}

		// reverse
		// u_rkmer = g->kmer_beg[k];
		u_rkmer = kmer_beg[k];
		__get_revc_num(u_rkmer, u_kmer, ksize, kmask);
		if (u_kmer < u_rkmer) {
			u_node = (gint_t)bin_search_id(pre_g, u_kmer);
			u_adj = nodes[u_node].adj & (uint8_t)0xf;
		} else {
			u_node = (gint_t)bin_search_id(pre_g, u_rkmer);
			u_adj = nodes[u_node].adj >> 4;
		}
		u_deg = __degree(u_adj);
		if (u_deg) {
			deg = 0;
			g->radj[k] = malloc(u_deg * sizeof(gint_t));
			for (c = 0; c < 4; ++c) {
				if (!((u_adj >> c) & 1))
					continue;
				v_kmer = ((u_kmer << 2) & kmask) | c;
				v_rkmer = (u_rkmer >> 2) | ((kmkey_t)(c ^ 3) << lmc);
				if (v_kmer < v_rkmer)
					v_node = (gint_t)bin_search_id(pre_g, v_kmer);
				else
					v_node = (gint_t)bin_search_id(pre_g, v_rkmer);
				uk = kmer_chain_id[v_node];
				assert(uk >= 0 && uk < n_v);
				if (v_rkmer == kmer_end[uk])
					g->radj[k][deg++] = -(uk + 1);
				else if (v_kmer == kmer_beg[uk])
					g->radj[k][deg++] = uk + 1;
				else
					assert(0);
			}
			assert(deg == u_deg);
			__set_degree(g->bin_rdeg, k, u_deg - 1);
			assert(__get_degree(g->bin_rdeg, k) == u_deg - 1);
		} else {
			g->radj[k] = NULL;
			__set_degree(g->bin_rdeg, k, 0);
			assert(__get_degree(g->bin_rdeg, k) == 0);
		}
	}

	g->n_v = n_v;
	__VERBOSE("[INFO] Number of node on kmer glued graph: %llu\n", (long long unsigned)n_v);

	free(visited);
	free(kmer_beg);
	free(kmer_end);
	free(kmer_chain_id);

	return g;
}

void dump_scrap_graph(struct scrap_graph_t *g, struct opt_count_t *opt)
{
	__VERBOSE_LOG("RESULT", "Number of current graph node: %lld\n", (long long)g->n_v);
	char path[1024];
	strcpy(path, opt->out_dir);
	strcat(path, "/graph.gfa");
	FILE *fp = xfopen(path, "w");

	gint_t k, uk;
	int m_seq, ksize;
	char *seq;
	uint8_t c, deg;
	m_seq = 0x100;
	seq = malloc(m_seq);

	ksize = opt->kmer_master;

	for (k = 0; k < g->n_v; ++k) {
		if (g->seq_len[k] + 1 > m_seq) {
			m_seq = g->seq_len[k] + 1;
			seq = realloc(seq, m_seq);
		}
		dump_bin_seq(g->node_seq[k], seq, g->seq_len[k]);
		fprintf(fp, "S\t%lld\t%s\tKC:i:%lld\n", (long long)k + 1, seq,
						(long long)g->kmer_count[k]);
	}

	for (k = 0; k < g->n_v; ++k) {
		// forward edge
		if (g->fadj[k] != NULL) {
			deg = __get_degree(g->bin_fdeg, k) + 1;
			for (c = 0; c < deg; ++c) {
				uk = g->fadj[k][c];
				if (uk > 0) {
					fprintf(fp, "L\t%lld\t+\t%lld\t+\t%lldM\n",
						(long long)k + 1, (long long)uk,
						(long long)ksize - 1);
				} else {
					fprintf(fp, "L\t%lld\t+\t%lld\t-\t%lldM\n",
						(long long)k + 1, (long long)(-uk),
						(long long)ksize - 1);
				}
			}
		}
		// reverse edge
		if (g->radj[k] != NULL) {
			deg = __get_degree(g->bin_rdeg, k) + 1;
			for (c = 0; c < deg; ++c) {
				uk = g->radj[k][c];
				if (uk > 0) {
					fprintf(fp, "L\t%lld\t-\t%lld\t+\t%lldM\n",
						(long long)k + 1, (long long)uk,
						(long long)ksize - 1);
				} else {
					fprintf(fp, "L\t%lld\t-\t%lld\t-\t%lldM\n",
						(long long)k + 1, (long long)(-uk),
						(long long)ksize - 1);
				}
			}
		}
	}

	free(seq);

	fclose(fp);
}

static inline void atomic_on_bit_uint8(uint8_t *ptr, int pos)
{
	uint8_t old_bin, cur_bin, new_bin;
	cur_bin = *(volatile uint8_t *)ptr;
	do {
		old_bin = cur_bin;
		new_bin = cur_bin | (1 << pos);
		cur_bin = __sync_val_compare_and_swap(ptr, old_bin, new_bin);
	} while (cur_bin != old_bin);
}

static kmint_t bin_search_id(struct raw_graph_t *g, kmkey_t x)
{
	struct raw_node_t *nodes;
	kmint_t l, r, mid;
	l = 0;
	r = g->n;
	nodes = g->nodes;
	while (l < r) {
		mid = l + ((r - l) >> 1);
		if (nodes[mid].kmer < x)
			l = mid + 1;
		else
			r = mid;
	}
	if (l < g->n && nodes[l].kmer == x)
		return l;
	return g->n;
}

void count_edge(struct read_t *r, struct raw_graph_t *g, int ksize)
{
	int i, last, ci, ck, len, lmc, kedge;
	char *seq;
	len = r->len;
	seq = r->seq;

	kmkey_t knum, krev, pknum, pkrev, kmask;
	kmint_t ki, kk;
	kmask = ((kmkey_t)1 << (ksize << 1)) - 1;
	knum = krev = 0;
	last = 0;
	lmc = (ksize - 1) << 1;
	kedge = ksize + 1;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		knum = (knum << 2) & kmask;
		krev = krev >> 2;
		if (ci < 4) {
			knum |= ci;
			krev |= (kmkey_t)(ci ^ 3) << lmc;
			++last;
		} else {
			last = 0;
		}
		if (last >= kedge) {
			ck = nt4_table[(int)seq[i - ksize]] ^ 3;
			if (pknum < pkrev) {
				ki = bin_search_id(g, pknum);
			} else {
				ki = bin_search_id(g, pkrev);
				ci += 4;
			}

			if (knum < krev) {
				kk = bin_search_id(g, knum);
				ck += 4;
			} else {
				kk = bin_search_id(g, krev);
			}
			if (ki != g->n && kk != g->n) {
				// fprintf(stderr, "adding edge %llu %llu\n", (long long unsigned)ki, (long long unsigned)kk);
				atomic_on_bit_uint8(&(g->nodes[ki].adj), ci);
				atomic_on_bit_uint8(&(g->nodes[kk].adj), ck);
			}
		}
		pknum = knum;
		pkrev = krev;
	}
}

void *PE_edge_constructer(void *data)
{
	struct edgecount_bundle_t *bundle = (struct edgecount_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct raw_graph_t *g = bundle->g;
	int ksize = bundle->ksize;

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	int64_t n_reads;
	int64_t *gcnt_reads = bundle->n_reads;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf)
			break;
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		buf1 = ext_buf->buf1;
		buf2 = ext_buf->buf2;
		input_format = ext_buf->input_format;

		n_reads = 0;
		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, buf1, &pos1) :
				get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, buf2, &pos2) :
				get_read_from_fa(&read2, buf2, &pos2);


			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file\n");

			++n_reads;
			count_edge(&read1, g, ksize);
			count_edge(&read2, g, ksize);

			if (rc1 == READ_END)
				break;
		}
		n_reads = __sync_add_and_fetch(gcnt_reads, n_reads);
		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

void get_edges(struct opt_count_t *opt, struct raw_graph_t *g)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	struct producer_bundle_t *producer_bundles;
	producer_bundles = init_fastq_PE(opt);

	struct edgecount_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct edgecount_bundle_t));

	int64_t n_reads;
	n_reads = 0;

	int i;
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].g = g;
		worker_bundles[i].ksize = opt->kmer_master;
		worker_bundles[i].n_reads = &n_reads;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, PE_edge_constructer,
				worker_bundles + i);

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	free_fastq_PE(producer_bundles, opt->n_files);
	free(worker_bundles);

	free(producer_threads);
	free(worker_threads);
}

struct raw_graph_t *extract_kmer(struct opt_count_t *opt, struct kmhash_t *h)
{
	struct raw_graph_t *g;
	kmint_t n, k;
	kmval_t threshold;
	g = calloc(1, sizeof(struct raw_graph_t));
	g->nodes = malloc(h->n_items * sizeof(struct raw_node_t));
	n = 0;
	threshold = opt->filter_thres;
	for (k = 0; k < h->size; ++k) {
		if (h->keys[k] == TOMB_STONE || h->vals[k] <= threshold)
			continue;
		g->nodes[n].kmer = h->keys[k];
		g->nodes[n].cnt = h->vals[k];
		g->nodes[n].adj = 0;
		++n;
	}
	g->nodes = realloc(g->nodes, n * sizeof(struct raw_node_t));
	g->n = n;
	__VERBOSE("\n");
	__VERBOSE_LOG("KMER COUNT", "Number of %d-mer: %llu\n", opt->kmer_master,
		(long long unsigned)h->n_items);

	__VERBOSE_LOG("KMER_COUNT", "Number of %d-mer with count greater than (%d): %llu\n",
		opt->kmer_master, opt->filter_thres, (long long unsigned)n);
	return g;
}

static inline void kmer_insertion_sort(struct raw_node_t *b, struct raw_node_t *e)
{
	struct raw_node_t *i, *j, tmp;
	for (i = b + 1; i <  e; ++i) {
		if (i->kmer < (i - 1)->kmer) {
			tmp = *i;
			for (j = i; j > b && tmp.kmer < (j - 1)->kmer; j--)
				*j = *(j - 1);
			*j = tmp;
		}
	}
}

static void kmer_merge_sort(struct raw_node_t *a, struct raw_node_t *tmp,
				kmint_t l, kmint_t r, kmint_t m)
{
	struct raw_node_t *a1, *a2;
	kmint_t len1, len2, i1, i2, k;
	len1 = m - l;
	len2 = r - m;
	memcpy(tmp, a + l, len1 * sizeof(struct raw_node_t));

	a1 = tmp;
	a2 = a + m;
	a = a + l;

	i1 = i2 = k = 0;
	while (i1 < len1 && i2 < len2) {
		if (a1[i1].kmer < a2[i2].kmer)
			a[k++] = a1[i1++];
		else
			a[k++] = a2[i2++];
	}

	if (i1 < len1)
		memcpy(a + k, a1 + i1, (len1 - i1) * sizeof(struct raw_node_t));

	if (i2 < len2)
		memcpy(a + k, a2 + i2, (len2 - i2) * sizeof(struct raw_node_t));
}

void sort_kmer(struct raw_graph_t *g)
{
	struct raw_node_t *tmp, *nodes;
	kmint_t n, m, i, l, r, mid;
	m = n = g->n;
	__round_up_kmint(m);
	tmp = malloc((m >> 1) * sizeof(struct raw_node_t));
	nodes = g->nodes;

	for (i = 0; i < n; i += BUCK_SORT_SIZE)
		kmer_insertion_sort(nodes + i, nodes + __min(i + BUCK_SORT_SIZE, n));

	for (i = BUCK_SORT_SIZE; i < n; i <<= 1)
	{
		m = i << 1;
		for (l = 0; l < n; l += m) {
			r = __min(l + m, n);
			mid = __min(l + i, r);
			if (r - mid)
				kmer_merge_sort(nodes, tmp, l, r, mid);
		}
	}
	free(tmp);
}

void test_sort_kmer(struct raw_graph_t *g)
{
	kmint_t i;
	for (i = 1; i < g->n; ++i) {
		if (g->nodes[i].kmer < g->nodes[i - 1].kmer) {
			__ERROR("Test kmer sort fail");
		}
	}
	__VERBOSE("[DEBUG] Kmer sort success\n");
}

void get_edge_stat(struct raw_graph_t *g)
{
	uint64_t cnt;
	uint8_t v;
	kmint_t k;
	cnt = 0;
	for (k = 0; k < g->n; ++k) {
		v = g->nodes[k].adj;
		cnt += (v & 1) + ((v >> 1) & 1) + ((v >> 2) & 1) + ((v >> 3) & 1)
			+ ((v >> 4) & 1) + ((v >> 5) & 1) + ((v >> 6) & 1) + ((v >> 7) & 1);
	}
	__VERBOSE("[DEBUG] Number of edges: %llu\n", (long long unsigned)cnt);
}

static void dump_bin_seq(uint32_t *bin, char *seq, int len)
{
	int i;
	for (i = 0; i < len; ++i)
		seq[i] = nt4_char[(bin[i >> 4] >> ((i & 15) << 1)) & 3];
	seq[len] = '\0';
}

static uint32_t *init_seq_bin(kmkey_t key, int l)
{
	int i, m;
	uint32_t *ret;
	kmkey_t mask;
	mask = ((kmkey_t)1 << 32) - 1;
	m = (l + 15) / 16;
	ret = calloc(m, sizeof(uint32_t));
	i = 0;
	while (key) {
		ret[i++] = key & mask;
		key >>= 32;
	}
	return ret;
}

static void rev_bin_seq(uint32_t *seq, int beg, int len)
{
	int i, k, v;
	uint32_t ck, cv;
	for (i = 0; i < len >> 1; ++i) {
		k = beg + i;
		v = beg + len - i - 1;
		ck = (seq[k >> 4] >> ((k & 15) << 1)) & 3;
		cv = (seq[v >> 4] >> ((v & 15) << 1)) & 3;
		// seq[k >> 4] &= (~((uint32_t)3 << ((k & 15) << 1)));
		// seq[v >> 4] &= (~((uint32_t)3 << ((v & 15) << 1)));
		seq[k >> 4] ^= ((cv ^ ck ^ 3) << ((k & 15) << 1));
		seq[v >> 4] ^= ((ck ^ cv ^ 3) << ((v & 15) << 1));
	}
	if (len & 1) {
		k = beg + ((len + 1) >> 1);
		seq[k >> 4] ^= (uint32_t)3 << ((k & 15) << 1);
	}
}
