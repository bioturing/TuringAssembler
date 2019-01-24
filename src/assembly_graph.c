#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "atomic.h"
#include "attribute.h"
#include "fastq_producer.h"
#include "io_utils.h"
#include "kmer_count.h"
#include "kmhash.h"
#include "utils.h"
#include "verbose.h"
#include "time_utils.h"

/*
 * strand   :                <------------------
 * seq      :                X X X X X X X X X X
 * bit order:   10987654321098765432109876543210 -> easy for insert and extract
 *                                                           multibyte sequence
 * strand   :            ---------------------->
 * seq      :            X X X X X X X X X X X X
 * bit order:   10987654321098765432109876543210 -> easy for lexicographically
 *                                                      singlebyte comparation
 */

#define __get_rev_num(y, x, l)					       \
(	(x) = (y) << (64 - ((l) << 1)),					       \
	(x) = (((x) & 0xffffffff00000000ull) >> 32) | (((x) & 0x00000000ffffffffull) << 32), \
	(x) = (((x) & 0xffff0000ffff0000ull) >> 16) | (((x) & 0x0000ffff0000ffffull) << 16), \
	(x) = (((x) & 0xff00ff00ff00ff00ull) >>  8) | (((x) & 0x00ff00ff00ff00ffull) <<  8), \
	(x) = (((x) & 0xf0f0f0f0f0f0f0f0ull) >>  4) | (((x) & 0x0f0f0f0f0f0f0f0full) <<  4), \
	(x) = (((x) & 0xccccccccccccccccull) >>  2) | (((x) & 0x3333333333333333ull) <<  2))

static inline void deb_dump_bin_seq(const char *label, uint32_t *bin, uint32_t len)
{
	char *seq;
	int i;
	seq = malloc(len + 1);
	for (i = 0; i < len; ++i) {
		seq[i] = nt4_char[(bin[i >> 4] >> ((i & 15) << 1)) & (uint32_t)0x3];
	}
	seq[len] = '\0';
	fprintf(stderr, "%s%s\n", label, seq);
	free(seq);
}

static inline void deb_dump_seq(const char *label, kmkey_t kmer, int size)
{
	char *seq = malloc(size + 1);
	int i;
	for (i = 0; i < size; ++i) {
		seq[size - i - 1] = nt4_char[kmer & 3];
		kmer >>= 2;
	}
	seq[size] = '\0';
	fprintf(stderr, "%s%s\n", label, seq);
	free(seq);
}

static inline uint32_t *init_seq_bin(kmkey_t key, int l)
{
	int i, m;
	uint32_t *ret;
	kmkey_t mask, tmp;
	__get_rev_num(key, tmp, l);
	mask = ((kmkey_t)1 << 32) - 1;
	m = (l + 15) / 16;
	ret = calloc(m, sizeof(uint32_t));
	i = 0;
	while (tmp) {
		ret[i++] = tmp & mask;
		tmp >>= 32;
	}
	return ret;
}

static void dump_bin_seq(char *seq, uint32_t *bin, int len)
{
	int i;
	for (i = 0; i < len; ++i)
		seq[i] = nt4_char[(bin[i >> 4] >> ((i & 15) << 1)) & 3];
	seq[len] = '\0';
}

#define __only_edge4(e) ((((e) >> 1) & 1) * 1 + (((e) >> 2) & 1) * 2 + (((e) >> 3) & 1) * 3)

#define __get_revc_num(y, x, l, mask)					       \
(	(x) = (y) << (64 - ((l) << 1)),					       \
	(x) = (((x) & 0xffffffff00000000ull) >> 32) | (((x) & 0x00000000ffffffffull) << 32), \
	(x) = (((x) & 0xffff0000ffff0000ull) >> 16) | (((x) & 0x0000ffff0000ffffull) << 16), \
	(x) = (((x) & 0xff00ff00ff00ff00ull) >>  8) | (((x) & 0x00ff00ff00ff00ffull) <<  8), \
	(x) = (((x) & 0xf0f0f0f0f0f0f0f0ull) >>  4) | (((x) & 0x0f0f0f0f0f0f0f0full) <<  4), \
	(x) = (((x) & 0xccccccccccccccccull) >>  2) | (((x) & 0x3333333333333333ull) <<  2), \
	(x) ^= 0xffffffffffffffffull, (x) &= (mask))

#define __bin_seq_get_char(seq, l) (((seq)[(l) >> 4] >> (((l) & 15) << 1)) & (uint32_t)0x3)

#define __degree4(e) (((e) & 1) + (((e) >> 1) & 1) + (((e) >> 2) & 1) + (((e) >> 3) & 1))

#define __get_only_edge4(adj) (((adj)[1] != -1) * 1 + ((adj)[2] != -1) * 2 + ((adj)[3] != -1) * 3)

#define __get_degree4(adj) (((adj)[0] != -1) + ((adj)[1] != -1) + ((adj)[2] != -1) + ((adj)[3] != -1))

#define table_exist(h, k) ((k) != (h)->size)

struct edge_table_t {
	/* map non-branching kmer to edge */
	kmint_t size;
	kmint_t n_probe;
	kmint_t n_item;

	kmkey_t *kmer;
	gint_t  *edge_id;
	uint8_t *adj;
};

struct node_table_t {
	/* map branching kmer to node */
	kmint_t size;
	kmint_t n_probe;
	kmint_t n_item;

	kmkey_t *kmer;
	gint_t  *node_id;
};

struct edgecount_bundle_t {
	struct dqueue_t *q;
	struct asm_graph_t *graph;
	struct edge_table_t *edict;
	struct node_table_t *ndict;
	int ksize;
	int64_t *n_reads;
};

void convert_table(struct edge_table_t *edge_dict, struct kmhash_t *kmer_hash);

struct node_table_t *init_node_table(kmint_t size);

void build_graph_from_kmer(struct asm_graph_t *graph, int ksize,
				struct edge_table_t *edge_dict,
				struct node_table_t *node_dict);

void count_edge(struct opt_count_t *opt, struct asm_graph_t *graph,
		struct edge_table_t *edict, struct node_table_t *ndict);

void remove_dead_end(struct asm_graph_t *graph);

void write_gfa(struct asm_graph_t *g, const char *path);

kmint_t hash_edge_table_get(struct edge_table_t *h, kmkey_t key)
{
	kmint_t step, i, n_probe, mask;
	kmkey_t k;
	mask = h->size - 1;
	k = __hash_int2(key);
	i = k & mask;
	n_probe = h->n_probe;
	if (h->kmer[i] == key)
		return i;
	if (h->kmer[i] == TOMB_STONE)
		return h->size;
	step = 0;
	do {
		++step;
		i = (i + step * (step + 1) / 2) & mask;
		if (h->kmer[i] == key)
			return i;
	} while (step < n_probe && h->kmer[i] != TOMB_STONE);
	return h->size;
}

kmint_t hash_node_table_get(struct node_table_t *h, kmkey_t key)
{
	kmint_t step, i, n_probe, mask;
	kmkey_t k;
	mask = h->size - 1;
	k = __hash_int2(key);
	i = k & mask;
	n_probe = h->n_probe;
	if (h->kmer[i] == key)
		return i;
	if (h->kmer[i] == TOMB_STONE)
		return h->size;
	step = 0;
	do {
		++step;
		i = (i + step * (step + 1) / 2) & mask;
		if (h->kmer[i] == key)
			return i;
	} while (step < n_probe && h->kmer[i] != TOMB_STONE);
	return h->size;
}

#define edge_table_get hash_edge_table_get
#define node_table_get hash_node_table_get

void check_hash_table(struct kmhash_t *h, int ksize)
{
	kmint_t i, k;
	kmkey_t mask, kmer, kmer_rc, nkmer, nkmer_rc;
	int lmc, c, rc;
	mask = ((kmkey_t)1 << (ksize << 1)) - 1;
	lmc = (ksize - 1) << 1;

	for (i = 0; i < h->size; ++i) {
		if (h->keys[i] == TOMB_STONE)
			continue;
		kmer = h->keys[i];
		__get_revc_num(kmer, kmer_rc, ksize, mask);
		for (c = 0; c < 4; ++c) {
			if (!((h->adjs[i] >> c) & 1))
				continue;
			nkmer = ((kmer << 2) & mask) | c;
			__get_revc_num(nkmer, nkmer_rc, ksize, mask);
			if (nkmer < nkmer_rc)
				k = hash_get(h, nkmer);
			else
				k = hash_get(h, nkmer_rc);
			assert(k != KMHASH_MAX_SIZE);
			rc = kmer >> lmc;
			if (nkmer < nkmer_rc)
				assert((h->adjs[k] >> ((rc ^ 3) + 4)) & 1);
			else
				assert((h->adjs[k] >> (rc ^ 3)) & 1);
		}

		for (c = 0; c < 4; ++c) {
			if (!((h->adjs[i] >> (c + 4)) & 1))
				continue;
			nkmer = ((kmer_rc << 2) & mask) | c;
			__get_revc_num(nkmer, nkmer_rc, ksize, mask);
			if (nkmer < nkmer_rc)
				k = hash_get(h, nkmer);
			else
				k = hash_get(h, nkmer_rc);
			assert(k != KMHASH_MAX_SIZE);
			rc = kmer_rc >> lmc;
			if (nkmer < nkmer_rc) {
				assert((h->adjs[k] >> ((rc ^ 3) + 4)) & 1);
			} else {
				assert((h->adjs[k] >> (rc ^ 3)) & 1);
			}
		}
	}
	__VERBOSE("\nTest hash table done\n");
}

void assembly_process(struct opt_count_t *opt)
{
	char path[1024];
	init_clock();

	struct kmhash_t *kmer_hash;
	__VERBOSE("Estimating kmer\n");
	kmer_hash = build_kmer_table_lazy(opt);
	check_hash_table(kmer_hash, opt->kmer_master);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Estimating kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("\nBuilding assembly graph\n");
	__VERBOSE("|--- Condensing non-branching path\n");
	struct edge_table_t *edge_dict;
	struct node_table_t *node_dict;
	struct asm_graph_t *asm_graph;

	edge_dict = calloc(1, sizeof(struct edge_table_t));
	convert_table(edge_dict, kmer_hash);
	kmhash_destroy(kmer_hash);

	node_dict = init_node_table(opt->hash_size - 1);
	asm_graph = calloc(1, sizeof(struct asm_graph_t));
	build_graph_from_kmer(asm_graph, opt->kmer_master, edge_dict, node_dict);
	test_graph_build(asm_graph);
	__VERBOSE_LOG("Graph #1", "Number of nodes: %lld\n", (long long)asm_graph->n_v);
	__VERBOSE_LOG("Graph #1", "Number of edges: %lld\n", (long long)asm_graph->n_e);

	__VERBOSE("|--- Estimating edge coverage\n");
	count_edge(opt, asm_graph, edge_dict, node_dict);
	clear_node_dict(node_dict);
	clear_edge_dict(edge_dict);

	strcpy(path, opt->out_dir);
	strcat(path, "/graph_0.gfa");
	write_gfa(asm_graph, path);

	__VERBOSE("\nRemoving dead-end\n");
	remove_dead_end(asm_graph);
	test_graph_build(asm_graph);

	strcpy(path, opt->out_dir);
	strcat(path, "/graph_1.gfa");
	write_gfa(asm_graph, path);
}

char sign_char[2] = {'+', '-'};

void write_gfa(struct asm_graph_t *g, const char *path)
{
	FILE *fp = xfopen(path, "w");
	gint_t k, rc_id, e_i, e_k, rc_i, rc_k;
	gint_t *adj_i, *adj_k;
	int s_i, s_k, ci, ck;
	char *seq;

	seq = NULL;

	for (k = 0; k < g->n_e; ++k) {
		rc_id = g->edges[k].rc_id;
		if (k > rc_id)
			continue;
		seq = realloc(seq, g->edges[k].seq_len + 1);
		dump_bin_seq(seq, g->edges[k].seq, g->edges[k].seq_len);
		fprintf(fp, "S\t%lld\t%s\tKC:i:%llu\n", (long long)k + 1, seq,
					(long long unsigned)g->edges[k].count);
	}

	for (k = 0; k < g->n_v; ++k) {
		adj_i = g->nodes[k].reverse_adj;
		adj_k = g->nodes[k].forward_adj;
		for (ci = 0; ci < 4; ++ci) {
			e_i = adj_i[ci];
			if (e_i == -1)
				continue;
			rc_i = g->edges[e_i].rc_id;
			if (e_i < rc_i) {
				e_i = e_i + 1;
				s_i = 1;
			} else {
				e_i = rc_i + 1;
				s_i = 0;
			}
			for (ck = 0; ck < 4; ++ck) {
				e_k = adj_k[ck];
				if (e_k == -1)
					continue;
				rc_k = g->edges[e_k].rc_id;
				if (e_k < rc_k) {
					e_k = e_k + 1;
					s_k = 0;
				} else {
					e_k = rc_k + 1;
					s_k = 1;
				}
				fprintf(fp, "L\t%lld\t%c\t%lld\t%c\t%dM\n",
					(long long)e_i, sign_char[s_i],
					(long long)e_k, sign_char[s_k],
					g->ksize);
				fprintf(fp, "L\t%lld\t%c\t%lld\t%c\t%dM\n",
					(long long)e_k, sign_char[s_k ^ 1],
					(long long)e_i, sign_char[s_i ^ 1],
					g->ksize);
			}
		}
	}
	fclose(fp);
}

uint32_t *append_bin_seq_forward(uint32_t *dst, gint_t dlen, uint32_t *src,
					gint_t slen, int skip)
{
	uint32_t *new_ptr;
	gint_t cur_m, m, i, k;
	cur_m = (dlen + 15) >> 4;
	m = (dlen + slen + 15) >> 4;
	if (m > cur_m) {
		new_ptr = realloc(dst, m * sizeof(uint32_t));
		memset(new_ptr + cur_m, 0, (m - cur_m) * sizeof(uint32_t));
	} else {
		new_ptr = dst;
	}
	for (i = skip; i < skip + slen; ++i) {
		k = dlen + (i - skip);
		new_ptr[k >> 4] |= ((src[i >> 4] >> ((i & 15) << 1)) & (uint32_t)3) << ((k & 15) << 1);
	}
	return new_ptr;
}

void clean_graph(struct asm_graph_t *graph)
{
	gint_t i;
	free(graph->nodes);
	graph->nodes = NULL;
	for (i = 0; i < graph->n_e; ++i)
		free(graph->edges[i].seq);
	free(graph->edges);
	graph->edges = NULL;
}

void remove_dead_end(struct asm_graph_t *graph)
{
	gint_t *new_id, *adj;
	uint32_t *edge_seq;
	uint64_t count;
	struct asm_node_t *nodes;
	struct asm_edge_t *edges;
	gint_t u, v, e_id, next_node, n_id, current_eid, current_nid;
	gint_t n_v, n_e, edge_seq_len, n_edge, rc_id;
	float sum_cov, cov;
	int c, current_c, deg_fw, deg_rv, ksize;
	ksize = graph->ksize;
	/* disconnect edge */
	for (u = 0; u < graph->n_v; ++u) {
		sum_cov = 0.0;
		// max_cov = 0;
		n_edge = 0;
		for (c = 0; c < 4; ++c) {
			e_id = graph->nodes[u].forward_adj[c];
			/* kedge = graph->ksize + 1 */
			if (e_id != -1) {
				cov = graph->edges[e_id].count * 1.0 /
					(graph->edges[e_id].seq_len - graph->ksize);
				sum_cov += cov;
				++n_edge;
				// max_cov = __max(max_cov, cov);
			}

			e_id = graph->nodes[u].reverse_adj[c];
			if (e_id != -1) {
				cov = graph->edges[e_id].count * 1.0 /
					(graph->edges[e_id].seq_len - graph->ksize);
				sum_cov += cov;
				++n_edge;
				// max_cov = __max(max_cov, cov);
			}
		}
		for (c = 0; c < 4; ++c) {
			e_id = graph->nodes[u].forward_adj[c];
			/* kedge = graph->ksize + 1 */
			if (e_id != -1) {
				cov = graph->edges[e_id].count * 1.0 /
					(graph->edges[e_id].seq_len - graph->ksize);
				if (cov / sum_cov < 0.1 && n_edge > 2) {
					next_node = graph->edges[e_id].target;
					if (next_node < 0)
						next_node = -next_node - 1;
					else
						next_node = next_node - 1;
					if (__get_degree4(graph->nodes[next_node].forward_adj)
						+ __get_degree4(graph->nodes[next_node].reverse_adj) == 1) {
						/* dead end */
						/* disconnect edge */
						memset(graph->nodes[next_node].forward_adj,
							255, 4 * sizeof(gint_t));
						memset(graph->nodes[next_node].reverse_adj,
							255, 4 * sizeof(gint_t));
						graph->nodes[u].forward_adj[c] = -1;
						--n_edge;
					}
				}
			}

			e_id = graph->nodes[u].reverse_adj[c];
			if (e_id != -1) {
				cov = graph->edges[e_id].count * 1.0 /
					(graph->edges[e_id].seq_len - graph->ksize);
				if (cov / sum_cov < 0.1 && n_edge > 2) {
					next_node = graph->edges[e_id].target;
					if (next_node < 0)
						next_node = -next_node - 1;
					else
						next_node = next_node - 1;
					if (__get_degree4(graph->nodes[next_node].forward_adj)
						+ __get_degree4(graph->nodes[next_node].reverse_adj) == 1) {
						/* dead end */
						/* disconnect edge */
						memset(graph->nodes[next_node].forward_adj,
							255, 4 * sizeof(gint_t));
						memset(graph->nodes[next_node].reverse_adj,
							255, 4 * sizeof(gint_t));
						graph->nodes[u].reverse_adj[c] = -1;
						--n_edge;
					}
				}
			}
		}
	}

	n_v = n_e = 0;
	new_id = malloc(graph->n_v * sizeof(gint_t));
	memset(new_id, 255, graph->n_v * sizeof(gint_t));
	/* Assign new node id */
	for (u = 0; u < graph->n_v; ++u) {
		deg_fw = __get_degree4(graph->nodes[u].forward_adj);
		deg_rv = __get_degree4(graph->nodes[u].reverse_adj);
		if ((deg_fw == 1 && deg_rv == 1) || deg_fw + deg_rv == 0)
			continue;
		new_id[u] = n_v++;
		n_e += (deg_fw + deg_rv);
	}

	/* Get list of new node */
	nodes = calloc(n_v, sizeof(struct asm_node_t));
	edges = calloc(n_e, sizeof(struct asm_edge_t));
	n_e = 0;
	for (u = 0; u < graph->n_v; ++u) {
		if (new_id[u] == -1)
			continue;
		n_id = new_id[u];
		nodes[n_id].seq = graph->nodes[u].seq;
		/* forward edges */
		for (c = 0; c < 4; ++c) {
			current_eid = graph->nodes[u].forward_adj[c];
			if (current_eid == -1) {
				nodes[n_id].forward_adj[c] = -1;
				continue;
			}
			edge_seq = graph->edges[current_eid].seq;
			graph->edges[current_eid].seq = NULL;
			edge_seq_len = graph->edges[current_eid].seq_len;
			count = graph->edges[current_eid].count;
			do {
				current_nid = graph->edges[current_eid].target;
				if (current_nid > 0) {
					current_nid = current_nid - 1;
					if (new_id[current_nid] == -1) {
						current_c = __get_only_edge4(
							graph->nodes[current_nid].forward_adj);
						current_eid = graph->nodes[current_nid].forward_adj[current_c];
						edge_seq = append_bin_seq_forward(edge_seq, edge_seq_len,
							graph->edges[current_eid].seq,
							graph->edges[current_eid].seq_len - ksize,
							ksize);
						edge_seq_len += (graph->edges[current_eid].seq_len - ksize);
						count += graph->edges[current_eid].count;
					}
				} else {
					current_nid = -current_nid - 1;
					if (new_id[current_nid] == -1) {
						current_c = __get_only_edge4(
							graph->nodes[current_nid].reverse_adj);
						current_eid = graph->nodes[current_nid].reverse_adj[current_c];
						edge_seq = append_bin_seq_forward(edge_seq, edge_seq_len,
							graph->edges[current_eid].seq,
							graph->edges[current_eid].seq_len - ksize,
							ksize);
						edge_seq_len += (graph->edges[current_eid].seq_len - ksize);
						count += graph->edges[current_eid].count;
					}
				}
			} while (new_id[current_nid] == -1);

			edges[n_e].seq = edge_seq;
			edges[n_e].seq_len = edge_seq_len;
			edges[n_e].source = n_id + 1;
			edges[n_e].count = count;
			if (graph->edges[current_eid].target < 0)
				edges[n_e].target = -new_id[current_nid] - 1;
			else
				edges[n_e].target = new_id[current_nid] + 1;
			nodes[n_id].forward_adj[c] = n_e;
			++n_e;
		}
		/* reverse edges */
		for (c = 0; c < 4; ++c) {
			current_eid = graph->nodes[u].reverse_adj[c];
			if (current_eid == -1) {
				nodes[n_id].reverse_adj[c] = -1;
				continue;
			}
			edge_seq = graph->edges[current_eid].seq;
			graph->edges[current_eid].seq = NULL;
			edge_seq_len = graph->edges[current_eid].seq_len;
			count = graph->edges[current_eid].count;
			do {
				current_nid = graph->edges[current_eid].target;
				if (current_nid > 0) {
					current_nid = current_nid - 1;
					if (new_id[current_nid] == -1) {
						current_c = __get_only_edge4(
							graph->nodes[current_nid].forward_adj);
						current_eid = graph->nodes[current_nid].forward_adj[current_c];
						edge_seq = append_bin_seq_forward(edge_seq, edge_seq_len,
							graph->edges[current_eid].seq,
							graph->edges[current_eid].seq_len - ksize,
							ksize);
						edge_seq_len += (graph->edges[current_eid].seq_len - ksize);
						count += graph->edges[current_eid].count;
					}
				} else {
					current_nid = -current_nid - 1;
					if (new_id[current_nid] == -1) {
						current_c = __get_only_edge4(
							graph->nodes[current_nid].reverse_adj);
						current_eid = graph->nodes[current_nid].reverse_adj[current_c];
						edge_seq = append_bin_seq_forward(edge_seq, edge_seq_len,
							graph->edges[current_eid].seq,
							graph->edges[current_eid].seq_len - ksize,
							ksize);
						edge_seq_len += (graph->edges[current_eid].seq_len - ksize);
						count += graph->edges[current_eid].count;
					}
				}
			} while (new_id[current_nid] == -1);

			edges[n_e].seq = edge_seq;
			edges[n_e].seq_len = edge_seq_len;
			edges[n_e].source = -n_id - 1;
			edges[n_e].count = count;
			if (graph->edges[current_eid].target < 0)
				edges[n_e].target = -new_id[current_nid] - 1;
			else
				edges[n_e].target = new_id[current_nid] + 1;
			nodes[n_id].reverse_adj[c] = n_e;
			++n_e;
		}
	}

	/* link reverse complemented edges */
	for (u = 0; u < n_e; ++u) {
		c = __bin_seq_get_char(edges[u].seq, edges[u].seq_len - ksize - 1) ^ 3;
		v = edges[u].target;
		if (v > 0)
			adj = nodes[v - 1].reverse_adj;
		else
			adj = nodes[-v - 1].forward_adj;
		assert(adj[c] != -1);
		edges[u].rc_id = adj[c];
	}

	/* double check to test whether reverse complemented edges are correctly
	 * labeled
	 */
	for (u = 0; u < n_e; ++u) {
		rc_id = edges[u].rc_id;
		assert(edges[rc_id].rc_id == u);
	}
	clean_graph(graph);
	graph->n_v = n_v;
	graph->n_e = n_e;
	graph->nodes = nodes;
	graph->edges = edges;
	free(new_id);
}

int is_reverse_complement(uint32_t *a, uint32_t *b, gint_t l)
{
	gint_t i;
	uint32_t ca, cb;
	for (i = 0; i < l; ++i) {
		ca = __bin_seq_get_char(a, i);
		cb = __bin_seq_get_char(b, l - i - 1);
		if (ca != (cb ^ 3)) {
			fprintf(stderr, "ca = %u; cb = %u\n", ca, cb);
			return 0;
		}
	}
	return 1;
}

void test_graph_build(struct asm_graph_t *graph)
{
	gint_t u, rc_id;
	for (u = 0; u < graph->n_e; ++u) {
		rc_id = graph->edges[u].rc_id;
		assert(graph->edges[rc_id].rc_id == u);
		assert(graph->edges[u].seq_len == graph->edges[rc_id].seq_len);
		if (!(is_reverse_complement(graph->edges[u].seq, graph->edges[rc_id].seq,
							graph->edges[u].seq_len))) {
			deb_dump_bin_seq("forward edge: ", graph->edges[u].seq, graph->edges[u].seq_len);
			deb_dump_bin_seq("reverse edge: ", graph->edges[rc_id].seq, graph->edges[rc_id].seq_len);
			assert(0);
		}
		assert(is_reverse_complement(graph->edges[u].seq, graph->edges[rc_id].seq,
							graph->edges[u].seq_len));
	}
	__VERBOSE("Test graph structure ok\n");
}

void clear_node_dict(struct node_table_t *h)
{
	free(h->kmer);
	free(h->node_id);
	h->n_item = 0;
}

void clear_edge_dict(struct edge_table_t *h)
{
	free(h->kmer);
	free(h->edge_id);
	free(h->adj);
	h->n_item = 0;
}

void reninit_node_table(struct node_table_t *h, kmint_t size)
{
	h->size = size - 1;
	__round_up_kmint(h->size);
	h->n_probe = estimate_probe_3(h->size);
	h->n_item = 0;
	h->kmer = malloc(h->size * sizeof(kmkey_t));
	h->node_id = calloc(h->size, sizeof(gint_t));
	kmint_t i;
	for (i = 0; i < h->size; ++i)
		h->kmer[i] = TOMB_STONE;
}

void reinit_edge_table(struct edge_table_t *h, kmint_t size)
{
	h->size = size - 1;
	__round_up_kmint(h->size);
	h->n_probe = estimate_probe_3(h->size);
	h->n_item = 0;
	h->kmer = malloc(h->size * sizeof(kmkey_t));
	h->edge_id = calloc(h->size, sizeof(gint_t));
	kmint_t i;
	for (i = 0; i < h->size; ++i)
		h->kmer[i] = TOMB_STONE;
}

struct node_table_t *init_node_table(kmint_t size)
{
	struct node_table_t *ret;
	ret = calloc(1, sizeof(struct node_table_t));
	ret->size = size;
	__round_up_kmint(ret->size);
	ret->n_probe = estimate_probe_3(ret->size);
	ret->n_item = 0;
	ret->kmer = malloc(ret->size * sizeof(kmkey_t));
	ret->node_id = calloc(ret->size, sizeof(gint_t));
	kmint_t i;
	for (i = 0; i < ret->size; ++i)
		ret->kmer[i] = TOMB_STONE;
	return ret;
}

void convert_table(struct edge_table_t *edge_dict, struct kmhash_t *kmer_hash)
{
	edge_dict->size = kmer_hash->size;
	edge_dict->n_probe = kmer_hash->n_probe;
	edge_dict->n_item = kmer_hash->n_item;

	edge_dict->kmer = kmer_hash->keys;
	edge_dict->adj = kmer_hash->adjs;
	kmer_hash->keys = NULL;
	kmer_hash->adjs = NULL;
	edge_dict->edge_id = calloc(edge_dict->size, sizeof(gint_t));
}

static kmint_t node_table_internal_put(struct node_table_t *h, kmkey_t key)
{
	kmint_t mask, step, i, n_probe;
	kmkey_t cur_key, k;

	mask = h->size - 1;
	k = __hash_int2(key);
	n_probe = h->n_probe;
	i = k & mask;

	if (h->kmer[i] == TOMB_STONE) {
		h->kmer[i] = key;
		++h->n_item;
		return i;
	}

	step = 0;
	do {
		++step;
		i = (i + (step * (step + 1) / 2)) & mask;
		if (h->kmer[i] == TOMB_STONE) {
			h->kmer[i] = key;
			++h->n_item;
			return i;
		}
	} while (step < n_probe);
	return h->size;
}

#define KMMASK_EMPTY			0
#define KMMASK_OLD			1
#define KMMASK_NEW			2
#define KMMASK_LOADING			3

#define rs_set_old(x, i) ((x)[(i) >> 4] = ((x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1)))) |       \
				((uint32_t)KMMASK_OLD << (((i) & 15) << 1)))

#define rs_set_new(x, i) ((x)[(i) >> 4] = ((x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1)))) |       \
				((uint32_t)KMMASK_NEW << (((i) & 15) << 1)))

#define rs_set_empty(x, i) ((x)[(i) >> 4] = (x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1))))

#define rs_get_flag(x, i) (((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)

#define rs_is_old(x, i) ((((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)  \
							== (uint32_t)KMMASK_OLD)

static inline kmint_t recount_table(struct node_table_t *h)
{
	kmint_t s, i;
	s = 0;
	for (i = 0; i < h->size; ++i)
		if (h->kmer[i] != TOMB_STONE)
			++s;
	return s;
}

static inline kmkey_t calhash_table(struct node_table_t *h)
{
	kmkey_t sum;
	kmint_t i;
	sum = 0;
	for (i = 0; i < h->size; ++i)
		sum += (h->kmer[i] ^ (kmkey_t)h->node_id[i]);
	return sum;
}

void node_table_resize(struct node_table_t *h)
{
	kmkey_t old_hash = calhash_table(h);
	kmint_t i, old_size, j, step, n_probe, mask;
	kmkey_t k, x, xt;
	gint_t y, yt;
	uint32_t *flag;
	uint32_t current_flag;
	old_size = h->size;
	h->size <<= 1;
	h->n_probe = estimate_probe_3(h->size);
	h->kmer = realloc(h->kmer, h->size * sizeof(kmkey_t));
	h->node_id = realloc(h->node_id, h->size * sizeof(gint_t));
	flag = calloc(h->size >> 4, sizeof(uint32_t));

	mask = h->size - 1;
	n_probe = h->n_probe;

	for (i = old_size; i < h->size; ++i) {
		h->kmer[i] = TOMB_STONE;
		h->node_id[i] = 0;
	}

	for (i = 0; i < old_size; ++i) {
		if (h->kmer[i] != TOMB_STONE)
			rs_set_old(flag, i);
	}

	for (i = 0; i < old_size; ++i) {
		if (rs_is_old(flag, i)) {
			x = h->kmer[i];
			y = h->node_id[i];
			rs_set_new(flag, i);
			h->kmer[i] = TOMB_STONE;
			h->node_id[i] = 0;
			while (1) {
				k = __hash_int2(x);
				j = k & mask;
				step = 0;
				current_flag = KMMASK_NEW;

				while (step <= n_probe) {
					j = (j + step * (step + 1) / 2) & mask;
					current_flag = rs_get_flag(flag, j);
					if (current_flag == KMMASK_EMPTY) {
						rs_set_new(flag, j);
						h->kmer[j] = x;
						h->node_id[j] = y;
						break;
					} else if (current_flag == KMMASK_OLD) {
						rs_set_new(flag, j);
						xt = h->kmer[j];
						yt = h->node_id[j];
						h->kmer[j] = x;
						h->node_id[j] = y;
						x = xt;
						y = yt;
						break;
					}
					++step;
				}
				if (current_flag == KMMASK_EMPTY)
					break;
				else if (current_flag == KMMASK_NEW)
					__ERROR("Resizing node table error");
			}
		}
	}

	assert(recount_table(h) == h->n_item);
	assert(calhash_table(h) == old_hash);

	free(flag);
}

void node_table_put(struct node_table_t *h, kmkey_t key, gint_t node_id)
{
	kmint_t k;
	k = node_table_internal_put(h, key);
	if (k != h->size)
		h->node_id[k] = node_id;
	while (k == h->size) {
		node_table_resize(h);
		k = node_table_internal_put(h, key);
		if (k != h->size)
			h->node_id[k] = node_id;
	}
}

void build_graph_from_kmer(struct asm_graph_t *graph, int ksize,
				struct edge_table_t *edge_dict,
				struct node_table_t *node_dict)
{
	struct asm_node_t *nodes;
	struct asm_edge_t *edges;

	uint32_t *edge_seq;
	gint_t *adj;

	gint_t n_e, n_v, u, v, rc_id;
	kmint_t i, j, k, node_id;
	kmkey_t node_kmer, node_kmer_rc, cur_kmer, cur_kmer_rc, kmask, hash_sum;
	uint32_t edge_seq_len;
	int lmc;
	uint8_t adj_forward, adj_reverse, deg_forward, deg_reverse, cur_c, c;

	kmask = ((kmkey_t)1 << (ksize << 1)) - 1;
	lmc = (ksize - 1) << 1;

	__VERBOSE("Counting node and edge\n");
	/* Count node and edge */
	n_e = n_v = 0;
	hash_sum = 0;
	for (i = 0; i < edge_dict->size; ++i) {
		if (edge_dict->kmer[i] == TOMB_STONE)
			continue;
		adj_forward = edge_dict->adj[i] & 0xf;
		adj_reverse = edge_dict->adj[i] >> 4;
		deg_forward = __degree4(adj_forward);
		deg_reverse = __degree4(adj_reverse);
		if (deg_forward == 1 && deg_reverse == 1)
			continue;
		node_table_put(node_dict, edge_dict->kmer[i], n_v);
		++n_v;
		n_e += (deg_forward + deg_reverse);
		hash_sum += edge_dict->kmer[i];
	}

	fprintf(stderr, "node hash sum = %llu\n", hash_sum);

	hash_sum = 0;
	gint_t *check = malloc(n_v * sizeof(gint_t));
	for (i = 0; i < n_v; ++i)
		check[i] = -1;

	for (i = 0; i < node_dict->size; ++i) {
		if (node_dict->kmer[i] == TOMB_STONE)
			continue;
		hash_sum += node_dict->kmer[i];
		if (check[node_dict->node_id[i]] != -1) {
			fprintf(stderr, "i = %llu; collide = %lld\n", i, check[node_dict->node_id[i]]);
			fprintf(stderr, "node_id_1 = %lld; node_id_2 = %lld\n",
				node_dict->node_id[i], node_dict->node_id[check[node_dict->node_id[i]]]);
			assert(0);
		}
		check[node_dict->node_id[i]] = i;
	}

	fprintf(stderr, "node hash sum recal = %llu\n", hash_sum);

	fprintf(stderr, "n_v = %llu\n", (long long unsigned)n_v);
	fprintf(stderr, "n_e = %llu\n", (long long unsigned)n_e);

	assert(recount_table(node_dict) == node_dict->n_item);
	assert(n_v == node_dict->n_item);

	nodes = calloc(n_v, sizeof(struct asm_node_t));
	edges = calloc(n_e, sizeof(struct asm_edge_t));

	__VERBOSE("Assembling node and edge\n");
	n_e = n_v = 0;
	/* Assemble node and edge */
	for (i = 0; i < node_dict->size; ++i) {
		if (node_dict->kmer[i] == TOMB_STONE)
			continue;
		node_kmer = node_dict->kmer[i];
		__get_revc_num(node_kmer, node_kmer_rc, ksize, kmask);
		node_id = node_dict->node_id[i];
		nodes[node_id].seq = node_kmer;
		k = edge_table_get(edge_dict, node_kmer);
		adj_forward = edge_dict->adj[k] & 0xf;
		adj_reverse = edge_dict->adj[k] >> 4;

		/* forward kmer's edges */
		for (c = 0; c < 4; ++c) {
			if (!((adj_forward >> c) & 1)) {
				nodes[node_id].forward_adj[c] = -1;
				continue;
			}
			edge_seq = init_seq_bin(node_kmer, ksize);
			edge_seq_len = ksize;
			cur_kmer = node_kmer;
			cur_kmer_rc = node_kmer_rc;
			cur_c = c;
			do {
				if (!(edge_seq_len & 15)) {
					edge_seq = realloc(edge_seq,
							((edge_seq_len >> 4) + 1) * sizeof(uint32_t));
					edge_seq[edge_seq_len >> 4] = 0;
				}
				edge_seq[edge_seq_len >> 4] |=
					(uint32_t)cur_c << ((edge_seq_len & 15) << 1);
				++edge_seq_len;

				cur_kmer = ((cur_kmer << 2) & kmask) | cur_c;
				cur_kmer_rc = (cur_kmer_rc >> 2) |
						((kmkey_t)(cur_c ^ 3) << lmc);
				if (cur_kmer < cur_kmer_rc) {
					k = node_table_get(node_dict, cur_kmer);
					if (!table_exist(node_dict, k)) {
						j = edge_table_get(edge_dict, cur_kmer);
						if (!table_exist(edge_dict, j))
							fprintf(stderr, "%llu %u\n",
								cur_kmer, edge_seq_len);
						assert(table_exist(edge_dict, j));
						edge_dict->edge_id[j] = n_e;
						assert(__degree4(edge_dict->adj[j] & 0xf) == 1);
						cur_c = __only_edge4(edge_dict->adj[j] & 0xf);
					}
				} else {
					k = node_table_get(node_dict, cur_kmer_rc);
					if (!table_exist(node_dict, k)) {
						j = edge_table_get(edge_dict, cur_kmer_rc);
						if (!table_exist(edge_dict, j))
							fprintf(stderr, "%llu %u\n",
								cur_kmer, edge_seq_len);
						assert(table_exist(edge_dict, j));
						assert(__degree4(edge_dict->adj[j] >> 4) == 1);
						cur_c = __only_edge4(edge_dict->adj[j] >> 4);
					}
				}
			} while (!table_exist(node_dict, k));

			edges[n_e].seq = edge_seq;
			edges[n_e].seq_len = edge_seq_len;
			edges[n_e].source = node_id + 1;
			if (cur_kmer < cur_kmer_rc)
				edges[n_e].target = node_dict->node_id[k] + 1;
			else
				edges[n_e].target = -node_dict->node_id[k] - 1;
			nodes[node_id].forward_adj[c] = n_e;
			++n_e;
		}

		/* reverse kmer's edges */
		for (c = 0; c < 4; ++c) {
			if (!((adj_reverse >> c) & 1)) {
				nodes[node_id].reverse_adj[c] = -1;
				continue;
			}
			edge_seq = init_seq_bin(node_kmer_rc, ksize);
			edge_seq_len = ksize;
			cur_kmer = node_kmer_rc;
			cur_kmer_rc = node_kmer;
			cur_c = c;
			do {
				if (!(edge_seq_len & 15)) {
					edge_seq = realloc(edge_seq,
							((edge_seq_len >> 4) + 1) * sizeof(uint32_t));
					edge_seq[edge_seq_len >> 4] = 0;
				}
				edge_seq[edge_seq_len >> 4] |= (uint32_t)cur_c << ((edge_seq_len & 15) << 1);
				++edge_seq_len;

				cur_kmer = ((cur_kmer << 2) & kmask) | cur_c;
				cur_kmer_rc = (cur_kmer_rc >> 2) |
						((kmkey_t)(cur_c ^ 3) << lmc);
				if (cur_kmer < cur_kmer_rc) {
					k = node_table_get(node_dict, cur_kmer);
					if (!table_exist(node_dict, k)) {
						j = edge_table_get(edge_dict, cur_kmer);
						if (!table_exist(edge_dict, j))
							fprintf(stderr, "%llu %u\n",
								cur_kmer, edge_seq_len);
						assert(table_exist(edge_dict, j));
						edge_dict->edge_id[j] = n_e;
						assert(__degree4(edge_dict->adj[j] & 0xf) == 1);
						cur_c = __only_edge4(edge_dict->adj[j] & 0xf);
					}
				} else {
					k = node_table_get(node_dict, cur_kmer_rc);
					if (!table_exist(node_dict, k)) {
						j = edge_table_get(edge_dict, cur_kmer_rc);
						if (!table_exist(edge_dict, j))
							fprintf(stderr, "%llu %u\n",
								cur_kmer, edge_seq_len);
						assert(table_exist(edge_dict, j));
						assert(__degree4(edge_dict->adj[j] >> 4) == 1);
						cur_c = __only_edge4(edge_dict->adj[j] >> 4);
					}
				}
			} while (!table_exist(node_dict, k));
			edges[n_e].seq = edge_seq;
			edges[n_e].seq_len = edge_seq_len;

			edges[n_e].source = -node_id - 1;
			if (cur_kmer < cur_kmer_rc)
				edges[n_e].target = node_dict->node_id[k] + 1;
			else
				edges[n_e].target = -node_dict->node_id[k] - 1;
			nodes[node_id].reverse_adj[c] = n_e;
			++n_e;
		}
		++n_v;
	}

	for (u = 0; u < n_v; ++u) {
		k = node_table_get(node_dict, nodes[u].seq);
		assert(table_exist(node_dict, k));
		assert(node_dict->node_id[k] == u);
	}

	fprintf(stderr, "n_v = %llu\n", (long long unsigned)n_v);
	fprintf(stderr, "n_e = %llu\n", (long long unsigned)n_e);

	assert(n_v == node_dict->n_item);

	/* link reverse complemented edges */
	for (u = 0; u < n_e; ++u) {
		c = __bin_seq_get_char(edges[u].seq, edges[u].seq_len - ksize - 1) ^ 3;
		v = edges[u].target;
		if (v > 0)
			adj = nodes[v - 1].reverse_adj;
		else
			adj = nodes[-v - 1].forward_adj;
		if (adj[c] == -1) {
			deb_dump_bin_seq("edge seq   : ", edges[u].seq, edges[u].seq_len);
			if (edges[u].target > 0) {
				fprintf(stderr, "target id : %lld; target seq: %llu\n",
					edges[u].target - 1, nodes[edges[u].target - 1].seq);
				deb_dump_seq("target seq : ",
					nodes[edges[u].target - 1].seq, ksize);
			} else {
				fprintf(stderr, "target id : %lld; target seq: %llu\n",
					-edges[u].target - 1, nodes[-edges[u].target - 1].seq);
				deb_dump_seq("target seq : ",
					nodes[-edges[u].target - 1].seq, ksize);
			}
			if (edges[u].source > 0) {
				fprintf(stderr, "source id : %lld; source seq: %llu\n",
					edges[u].source - 1, nodes[edges[u].source - 1].seq);
				deb_dump_seq("source seq : ",
					nodes[edges[u].source - 1].seq, ksize);
			} else {
				fprintf(stderr, "source id : %lld; source seq: %llu\n",
					-edges[u].source - 1, nodes[-edges[u].source - 1].seq);
				deb_dump_seq("source seq : ",
					nodes[-edges[u].source - 1].seq, ksize);
			}
			assert(0);
		}
		edges[u].rc_id = adj[c];
	}

	/* double check to test whether reverse complemented edges are correctly
	 * labeled
	 */
	for (u = 0; u < n_e; ++u) {
		rc_id = edges[u].rc_id;
		assert(edges[rc_id].rc_id == u);
	}

	graph->n_e = n_e;
	graph->n_v = n_v;
	graph->nodes = nodes;
	graph->edges = edges;
	graph->ksize = ksize;
}

void inc_count_edge(kmkey_t kmer, kmkey_t kmer_rc, int c, struct asm_graph_t *g,
		struct edge_table_t *edict, struct node_table_t *ndict)
{
	kmint_t e_pos, n_pos;
	gint_t e_id, e_id_rc, n_id;
	if (kmer < kmer_rc) {
		n_pos = node_table_get(ndict, kmer);
		if (table_exist(ndict, n_pos)) {
			n_id = ndict->node_id[n_pos];
			e_id = g->nodes[n_id].forward_adj[c];
			if (e_id != -1) {
				atomic_add_and_fetch64(&(g->edges[e_id].count), 1);
				e_id_rc = g->edges[e_id].rc_id;
				atomic_add_and_fetch64(&(g->edges[e_id_rc].count), 1);
			}
		} else {
			e_pos = edge_table_get(edict, kmer);
			if (!table_exist(edict, e_pos))
				return;
			e_id = edict->edge_id[e_pos];
			atomic_add_and_fetch64(&(g->edges[e_id].count), 1);
			e_id_rc = g->edges[e_id].rc_id;
			atomic_add_and_fetch64(&(g->edges[e_id_rc].count), 1);
		}
	} else {
		n_pos = node_table_get(ndict, kmer_rc);
		if (table_exist(ndict, n_pos)) {
			n_id = ndict->node_id[n_pos];
			e_id = g->nodes[n_id].reverse_adj[c];
			if (e_id != -1) {
				atomic_add_and_fetch64(&(g->edges[e_id].count), 1);
				e_id_rc = g->edges[e_id].rc_id;
				atomic_add_and_fetch64(&(g->edges[e_id_rc].count), 1);
			}
		} else {
			e_pos = edge_table_get(edict, kmer_rc);
			if (!table_exist(edict, e_pos))
				return;
			e_id = edict->edge_id[e_pos];
			atomic_add_and_fetch64(&(g->edges[e_id].count), 1);
			e_id_rc = g->edges[e_id].rc_id;
			atomic_add_and_fetch64(&(g->edges[e_id_rc].count), 1);
		}
	}
}

void count_edge_read(struct read_t *r, int ksize, struct asm_graph_t *graph,
			struct edge_table_t *edict, struct node_table_t *ndict)
{
	int i, last, last_i, ci, ck, len, lmc, kedge;
	char *seq;
	len = r->len;
	seq = r->seq;

	kmkey_t knum, krev, pknum, pkrev, kmask;
	kmint_t ki, kk;
	kmask = ((kmkey_t)1 << (ksize << 1)) - 1;
	knum = krev = 0;
	last = 0;
	last_i = -1;
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
			/*
			 * firstly, we find the edge
			 * increase count of edge
			 * and increase count of the reverse complemented edge
			 */
			inc_count_edge(pknum, pkrev, ci, graph, edict, ndict);
		}
		pknum = knum;
		pkrev = krev;
	}
}

void *count_edge_worker(void *data)
{

	struct edgecount_bundle_t *bundle = (struct edgecount_bundle_t *)data;

	struct dqueue_t *q = bundle->q;
	struct asm_graph_t *graph = bundle->graph;
	struct edge_table_t *edict = bundle->edict;
	struct node_table_t *ndict = bundle->ndict;

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	int64_t n_reads;
	int64_t *gcnt_reads;
	gcnt_reads = bundle->n_reads;

	int ksize = bundle->ksize;

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
			count_edge_read(&read1, ksize, graph, edict, ndict);
			count_edge_read(&read2, ksize, graph, edict, ndict);

			if (rc1 == READ_END)
				break;
		}
		n_reads = atomic_add_and_fetch64(gcnt_reads, n_reads);
		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

void count_edge(struct opt_count_t *opt, struct asm_graph_t *graph,
		struct edge_table_t *edict, struct node_table_t *ndict)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	int i;

	struct producer_bundle_t *producer_bundles;
	producer_bundles = init_fastq_PE(opt);

	struct edgecount_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct edgecount_bundle_t));

	int64_t n_reads;
	n_reads = 0;

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].graph = graph;
		worker_bundles[i].edict = edict;
		worker_bundles[i].ndict = ndict;
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
		pthread_create(worker_threads + i, &attr, count_edge_worker,
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

