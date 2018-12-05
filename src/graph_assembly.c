#include <stdlib.h>

#include "fastq_producer.h"
#include "graph_assembly.h"
#include "kmer_count.h"
#include "utils.h"
#include "verbose.h"

#define BUCK_SORT_SIZE		64

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

void sort_kmer(struct raw_graph_t *g);

void test_sort_kmer(struct raw_graph_t *g);

void get_edges(struct opt_count_t *opt, struct raw_graph_t *g);

void get_edge_stat(struct raw_graph_t *g);

void assembly_process(struct opt_count_t *opt)
{
	__VERBOSE("2-step kmer counting\n");
	struct kmhash_t *kmer_hash;
	kmer_hash = count_kmer(opt);

	struct raw_graph_t *pre_graph;
	pre_graph = extract_kmer(opt, kmer_hash);
	kmhash_destroy(kmer_hash);

	__VERBOSE("kmer graph building\n");
	sort_kmer(pre_graph);
	// test_sort_kmer(pre_graph);

	get_edges(opt, pre_graph);
	__VERBOSE("\n");
	// get_edge_stat(pre_graph);

	struct scrap_graph_t *scratch_graph;
	scratch_graph = sketch_graph(pre_graph);
}

struct scrap_graph_t *sketch_graph(struct raw_graph_t *pre_graph)
{
	struct scrap_graph_t *g;
	g = calloc(1, sizeof(struct scrap_graph_t));

	gint_t m_v, m_e, mchain, lchain;
	kmint_t i, n_k;
	struct raw_node_t *nodes;
	nodes = pre_graph->nodes;
	n_k = pre_graph->n;

	uint32_t *visited;
	visited = calloc((n_k + 31) / 32, sizeof(uint32_t));

	for (i = 0; i < n_k; ++i) {
		node_id = (gint_t)i;
		node_kmer = nodes[i].kmer;
		if (__get_bit(visited, node_id))
			continue;
		__get_revc_num(node_kmer, node_rkmer, ksize, kmask);
		assert(node_kmer < node_rkmer);

		if (n_v + 2 > m_v) {
			m_v <<= 1;
			g->kmer_count = realloc(g->kmer_count, m_v * sizeof(gint_t));
			g->kmer_chains = realloc(g->kmer_chains, m_v * sizeof(kmkey_t *));
		}

		chain = malloc(sizeof(kmkey_t));
		chain[0] = node_kmer;
		lchain = 1;
		g->kmer_count[n_v] = nodes[node_id].cnt;
		__on_bit(visited, node_id);

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
				v_node = bin_search_id(g, v_kmer);
				v_radj = nodes[v_node].adj >> 4;
			} else {
				v_node = bin_search_id(g, v_rkmer);
				v_radj = nodes[v_node].adj & (uint8_t)0xf;
			}
			// Check if node v is on another chain
			if (__get_bit(visited, v_node))
				break;
			// Check if 1-1 edge
			if (__degree(v_radj) != 1)
				break;
			assert(u_rkmer == (((v_rkmer << 2) & kmask | __only_edge(v_radj))));

			chain = realloc(chain, lchain + 1);
			chain[lchain++] = v_kmer;
			g->kmer_count[n_v] += nodes[v_node].cnt;
			__on_bit(visited, v_node);
			u_node = v_node;
			u_kmer = v_kmer;
			u_rkmer = v_rkmer;
			u_adj = nodes[u_node].adj >> (4 * (u_kmer > u_rkmer)) & (uint8_t)0xf;
		}

		old_lchain = lchain;
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
				v_node = bin_search_id(g, v_kmer);
				v_radj = nodes[v_node].adj >> 4;
			} else {
				v_node = bin_search_id(g, v_rkmer);
				v_radj = nodes[v_node].adj & (uint8_t)0xf;
			}
			// Check if node v is on another chain
			if (__get_bit(visited, v_node))
				break;
			// Check if 1-1 edge
			if (__degree(v_radj) != 1)
				break;
			assert(u_rkmer == (((v_rkmer << 2) & kmask | __only_edge(v_radj))));

			chain = realloc(chain, lchain + 1);
			chain[lchain++] = v_kmer;
			g->kmer_count[n_v] += nodes[v_node].cnt;
			__on_bit(visited, v_node);
			u_node = v_node;
			u_kmer = v_kmer;
			u_rkmer = v_rkmer;
			u_adj = nodes[u_node].adj >> (4 * (u_kmer > u_rkmer)) & (uint8_t)0xf;
		}

		// correct the chain
		if (lchain > old_lchain) {
			for (k = 0; k < (old_lchain >> 1); ++k) {
				__get_revc_num(chain[k], tmp, ksize, kmask);
				__get_revc_num(chain[old_lchain - k - 1], chain[k], ksize, kmask);
				chain[old_lchain - k - 1] = tmp;
			}
			if (old_lchain & 1) {
				__get_revc_num(chain[old_lchain >> 1], tmp, ksize, kmask);
				chain[old_lchain >> 1] = tmp;
			}
		}

		g->kmer_chains[n_v] = chain;
		+++n_v;
	}

	return g;
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
