#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "fastq_producer.h"
#include "k31_count.h"
#include "k31hash.h"
#include "utils.h"
#include "verbose.h"

struct edgecount_bundle_t {
	struct dqueue_t *q;
	struct asm_graph_t *graph;
	struct k63_idhash_t *edict;
	struct k63_idhash_t *ndict;
	int ksize;
	int64_t *n_reads;
};
#define __bin_degree4(e) (((e) & 1) + (((e) >> 1) & 1) + (((e) >> 2) & 1) + (((e) >> 3) & 1))

#define __bin_only4(e) ((((e) >> 1) & 1) * 1 + (((e) >> 2) & 1) * 2 + (((e) >> 3) & 1) * 3)

#define __bin_seq_get_char(seq, l) (((seq)[(l) >> 4] >> (((l) & 15) << 1)) & (uint32_t)0x3)

#define __k63_lt(x, y) ((x).bin[1] < (y).bin[1] || ((x).bin[1] == (y).bin[1] && (x).bin[0] < (y).bin[0]))

#define __k63_lshift2(k) (((k).bin[1] = ((k).bin[1] << 2) | ((k).bin[0] >> 62)), \
				((k).bin[0] <<= 2))
#define __k63_rshift2(k) (((k).bin[0] = ((k).bin[0] >> 2) | (((k).bin[1] & 0x3ull) << 62)), \
				((k).bin[1] >>= 2))

#define __k63_lshift(k, l) (((k).bin[1] = ((k).bin[1] << (l)) | ((k).bin[0] >> (64 - (l)))), \
				((k).bin[0] <<= (l)))
#define __k63_rshift(k, l) (((k).bin[0] = ((k).bin[0] >> (l)) | (((k).bin[1] & ((1ull << (l)) - 1)) << (64 - (l))), \
				((k).bin[1] >>= (l))))

#define __k63_and(k, v) ((k).bin[0] &= (v).bin[0], (k).bin[1] &= (v).bin[1])

#define __reverse_bit_order64(x)					       \
(									       \
	(x) = (((x) & 0xffffffff00000000ull) >> 32) | (((x) & 0x00000000ffffffffull) << 32), \
	(x) = (((x) & 0xffff0000ffff0000ull) >> 16) | (((x) & 0x0000ffff0000ffffull) << 16), \
	(x) = (((x) & 0xff00ff00ff00ff00ull) >>  8) | (((x) & 0x00ff00ff00ff00ffull) <<  8), \
	(x) = (((x) & 0xf0f0f0f0f0f0f0f0ull) >>  4) | (((x) & 0x0f0f0f0f0f0f0f0full) <<  4), \
	(x) = (((x) & 0xccccccccccccccccull) >>  2) | (((x) & 0x3333333333333333ull) <<  2)  \
)

#define __k63_revc_num(y, x, l, mask)					       \
(									       \
	(x) = (y), __k63_lshift(x, 128 - ((l) << 1)),			       \
	__reverse_bit_order64((x).bin[0]), __reverse_bit_order64((x).bin[1]),  \
	(x).bin[0] ^= (x).bin[1],					       \
	(x).bin[1] ^= (x).bin[0],					       \
	(x).bin[0] ^= (x).bin[1],					       \
	(x).bin[0] ^= 0xffffffffffffffffull, (x).bin[0] &= (mask).bin[0],      \
	(x).bin[1] ^= 0xffffffffffffffffull, (x).bin[1] &= (mask).bin[1]       \
)

#define __k63_rev_num(y, x, l)	\
(									       \
	(x) = (y), __k63_lshift(x, 128 - ((l) << 1)),			       \
	(x).bin[0] ^= (x).bin[1],					       \
	(x).bin[1] ^= (x).bin[0],					       \
	(x).bin[0] ^= (x).bin[1],					       \
	__reverse_bit_order64((x).bin[0]), __reverse_bit_order64((x).bin[1])   \
)

static inline uint32_t *k63_init_binseq(k63key_t key, int l)
{
	k63key_t tmp;
	__k63_rev_num(key, tmp, l);
	int m = (l + 15) / 16;
	uint32_t *ret = calloc(m, sizeof(uint32_t));
	uint32_t mask = (uint32_t)-1;
	int i = 0, k = 0;
	while (i < m) {
		ret[i++] = tmp.bin[k] & mask;
		if (i < m)
			ret[i++] = (tmp.bin[k++] >> 32) & mask;
	}
	return ret;
}

static inline void deb_dump_bin_seq(const char *label, uint32_t *bin, gint_t len)
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

static void k63_internal_build(int ksize, struct k63_idhash_t *edict,
			struct k63_idhash_t *ndict, struct asm_graph_t *g);

static void k63_count_edge(struct opt_count_t *opt, struct asm_graph_t *graph,
		struct k63_idhash_t *edict, struct k63_idhash_t *ndict);

void build_asm_graph_from_k63(struct opt_count_t *opt, int ksize,
			struct k63hash_t *kmer_hash, struct asm_graph_t *ret_g)
{
	struct k63_idhash_t *edict, *ndict;
	edict = calloc(1, sizeof(struct k63_idhash_t));
	ndict = calloc(1, sizeof(struct k63_idhash_t));
	k63_convert_table(kmer_hash, edict);
	k63hash_destroy(kmer_hash);

	k63_idhash_init(ndict, opt->hash_size - 1, 0);
	__VERBOSE("|----Condensing node and edge\n");
	k63_internal_build(ksize, edict, ndict, ret_g);
	test_asm_graph(ret_g);

	__VERBOSE("|----Estimating edges coveraage\n");
	k63_count_edge(opt, ret_g, edict, ndict);
	__VERBOSE("\n");
	k63_idhash_clean(edict);
	k63_idhash_clean(ndict);
	free(edict);
	free(ndict);
}

static void k63_internal_build(int ksize, struct k63_idhash_t *edict,
			struct k63_idhash_t *ndict, struct asm_graph_t *g)
{
	uint8_t adj_fw, adj_rv;
	int deg_fw, deg_rv;
	int lmc = (ksize << 1) - 2;
	k63key_t kmask;
	kmask.bin[0] = (uint64_t)-1;
	kmask.bin[1] = (1ull << ((ksize << 1) - 64)) - 1;
	gint_t n_v, n_e;

	n_v = n_e = 0;
	kmint_t i;
	for (i =  0; i < edict->size; ++i) {
		if (edict->flag[i] == KMFLAG_EMPTY)
			continue;
		adj_fw = IDHASH_ADJ(edict, i) & 0xf;
		adj_rv = IDHASH_ADJ(edict, i) >> 4;
		deg_fw = __bin_degree4(adj_fw);
		deg_rv = __bin_degree4(adj_rv);
		if ((deg_fw == 1 && deg_rv == 1) || deg_fw + deg_rv == 0)
			continue;
		kmint_t k = k63_idhash_put(ndict, IDHASH_KEY(edict, i));
		IDHASH_ID(ndict, k) = n_v++;
		n_e += (deg_fw + deg_rv);
	}

	assert(n_v == ndict->n_item);
	struct asm_node_t *nodes = calloc(n_v * 2, sizeof(struct asm_node_t));
	struct asm_edge_t *edges = calloc(n_e, sizeof(struct asm_edge_t));

	n_e = n_v = 0;
	for (i = 0; i < ndict->size; ++i) {
		if (ndict->flag[i] == KMFLAG_EMPTY)
			continue;
		k63key_t node_knum, node_krev, cur_knum, cur_krev;
		gint_t id_fw, id_rv;
		int c, cur_c;
		uint32_t *e_seq;
		gint_t *adj;
		gint_t e_len;

		node_knum = IDHASH_KEY(ndict, i);
		__k63_revc_num(node_knum, node_krev, ksize, kmask);
		gint_t node_id = IDHASH_ID(ndict, i);
		id_fw = node_id * 2;
		id_rv = node_id * 2 + 1;
		nodes[id_fw].rc_id = id_rv;
		nodes[id_rv].rc_id = id_fw;
		kmint_t k = k63_idhash_get(edict, node_knum);
		assert(k != IDHASH_END(edict));
		adj_fw = IDHASH_ADJ(edict, k) & 0xf;
		adj_rv = IDHASH_ADJ(edict, k) >> 4;
		deg_fw = __bin_degree4(adj_fw);
		deg_rv = __bin_degree4(adj_rv);
		nodes[id_fw].adj = malloc(deg_fw * sizeof(gint_t));
		nodes[id_rv].adj = malloc(deg_rv * sizeof(gint_t));
		nodes[id_fw].deg = deg_fw;
		nodes[id_rv].deg = deg_rv;

		deg_fw = deg_rv = 0;
		/* forward kmer's edges */
		for (c = 0; c < 4; ++c) {
			if (!((adj_fw >> c) & 1))
				continue;
			e_seq = k63_init_binseq(node_knum, ksize);
			e_len = ksize;
			cur_knum = node_knum;
			cur_krev = node_krev;
			cur_c = c;
			do {
				if (!(e_len & 15)) {
					e_seq = realloc(e_seq,
					((e_len >> 4) + 1) * sizeof(uint32_t));
					e_seq[e_len >> 4] = 0;
				}
				e_seq[e_len >> 4] |=
					(uint32_t)cur_c << ((e_len & 15) << 1);
				++e_len;
				__k63_lshift2(cur_knum);
				__k63_and(cur_knum, kmask);
				cur_knum.bin[0] |= cur_c;
				__k63_rshift2(cur_krev);
				cur_krev.bin[1] |= (uint64_t)(cur_c ^ 3) << (lmc - 64);
				// cur_knum = ((cur_knum << 2) & kmask) | cur_c;
				// cur_krev = (cur_krev >> 2) | ((k31key_t)(cur_c ^ 3) << lmc);
				kmint_t j;
				if (__k63_lt(cur_knum, cur_krev)) {
					k = k63_idhash_get(ndict, cur_knum);
					if (k == IDHASH_END(ndict)) {
						j = k63_idhash_get(edict, cur_knum);
						assert(j != IDHASH_END(edict));
						IDHASH_ID(edict, j) = n_e;
						assert(__bin_degree4(
							IDHASH_ADJ(edict, j) & 0xf) == 1);
						cur_c = __bin_only4(IDHASH_ADJ(edict, j) & 0xf);
					}
				} else {
					k = k63_idhash_get(ndict, cur_krev);
					if (k == IDHASH_END(ndict)) {
						j = k63_idhash_get(edict, cur_krev);
						assert(j != IDHASH_END(edict));
						assert(__bin_degree4(
							IDHASH_ADJ(edict, j) >> 4) == 1);
						cur_c = __bin_only4(IDHASH_ADJ(edict, j) >> 4);
					}
				}
			} while (k == IDHASH_END(ndict));
			edges[n_e].seq = e_seq;
			edges[n_e].seq_len = e_len;
			edges[n_e].source = id_fw;
			if (__k63_lt(cur_knum, cur_krev))
				edges[n_e].target = IDHASH_ID(ndict, k) * 2;
			else
				edges[n_e].target = IDHASH_ID(ndict, k) * 2 + 1;
			nodes[id_fw].adj[deg_fw++] = n_e;
			++n_e;
		}

		/* reverse kmer's edges */
		for (c = 0; c < 4; ++c) {
			if (!((adj_rv >> c) & 1))
				continue;
			e_seq = k63_init_binseq(node_krev, ksize);
			e_len = ksize;
			cur_knum = node_krev;
			cur_krev = node_knum;
			cur_c = c;
			do {
				if (!(e_len & 15)) {
					e_seq = realloc(e_seq,
					((e_len >> 4) + 1) * sizeof(uint32_t));
					e_seq[e_len >> 4] = 0;
				}
				e_seq[e_len >> 4] |=
					(uint32_t)cur_c << ((e_len & 15) << 1);
				++e_len;
				__k63_lshift2(cur_knum);
				__k63_and(cur_knum, kmask);
				cur_knum.bin[0] |= cur_c;
				__k63_rshift2(cur_krev);
				cur_krev.bin[1] |= (uint64_t)(cur_c ^ 3) << (lmc - 64);
				// cur_knum = ((cur_knum << 2) & kmask) | cur_c;
				// cur_krev = (cur_krev >> 2) | ((k31key_t)(cur_c ^ 3) << lmc);
				kmint_t j;
				if (__k63_lt(cur_knum, cur_krev)) {
					k = k63_idhash_get(ndict, cur_knum);
					if (k == IDHASH_END(ndict)) {
						j = k63_idhash_get(edict, cur_knum);
						assert(j != IDHASH_END(edict));
						IDHASH_ID(edict, j) = n_e;
						assert(__bin_degree4(
							IDHASH_ADJ(edict, j) & 0xf) == 1);
						cur_c = __bin_only4(IDHASH_ADJ(edict, j) & 0xf);
					}
				} else {
					k = k63_idhash_get(ndict, cur_krev);
					if (k == IDHASH_END(ndict)) {
						j = k63_idhash_get(edict, cur_krev);
						assert(j != IDHASH_END(edict));
						assert(__bin_degree4(
							IDHASH_ADJ(edict, j) >> 4) == 1);
						cur_c = __bin_only4(IDHASH_ADJ(edict, j) >> 4);
					}
				}
			} while (k == IDHASH_END(ndict));
			edges[n_e].seq = e_seq;
			edges[n_e].seq_len = e_len;
			edges[n_e].source = id_rv;
			if (__k63_lt(cur_knum, cur_krev))
				edges[n_e].target = IDHASH_ID(ndict, k) * 2;
			else
				edges[n_e].target = IDHASH_ID(ndict, k) * 2 + 1;
			nodes[id_rv].adj[deg_rv++] = n_e;
			++n_e;
		}
		++n_v;
	}

	gint_t e, e_rc, v, v_rc;
	for (e = 0; e < n_e; ++e) {
		edges[e].rc_id = -1;
		v = edges[e].target;
		v_rc = nodes[v].rc_id;
		int k;
		for (k = 0; k < nodes[v_rc].deg; ++k) {
			e_rc = nodes[v_rc].adj[k];
			if (edges[e_rc].target == nodes[edges[e].source].rc_id
				&& asm_is_edge_rc(edges[e].seq, edges[e].seq_len,
					edges[e_rc].seq, edges[e_rc].seq_len)) {
				edges[e].rc_id = e_rc;
				edges[e_rc].rc_id = e;
				break;
			}
		}
		if (edges[e].rc_id == -1) {
			deb_dump_bin_seq("edge seq: ", edges[e].seq, edges[e].seq_len);
			fprintf(stderr, "source = %lld; source reverse complement = %lld\n",
				edges[e].source, nodes[edges[e].source].rc_id);
			fprintf(stderr, "target = %lld; target reverse complement = %lld\n",
				v, v_rc);
			fprintf(stderr, "target rc deg = %lld\n", nodes[v_rc].deg);
			for (k = 0; k < nodes[v_rc].deg; ++k) {
				e_rc = nodes[v_rc].adj[k];
				fprintf(stderr, "edges: (id = %lld) (%lld) -> (%lld)\n",
					e_rc, edges[e_rc].source, edges[e_rc].target);
				deb_dump_bin_seq("edge seq: ", edges[e_rc].seq, edges[e_rc].seq_len);
			}
		}
		assert(edges[e].rc_id != -1);
	}

	g->n_e = n_e;
	g->n_v = n_v * 2;
	g->nodes = nodes;
	g->edges = edges;
	g->ksize = ksize;
}

static inline gint_t asm_find_edge(struct asm_edge_t *edges, gint_t *adj, int l,
					int deg, int c)
{
	int i, cur_c;
	gint_t e_id, ret;
	ret = -1;
	for (i = 0; i < deg; ++i) {
		e_id = adj[i];
		cur_c = __bin_seq_get_char(edges[e_id].seq, l);
		if (cur_c == c) {
			if (ret != -1)
				return -2;
			ret = e_id;
		}
	}
	return ret;
}

static void inc_count_edge(k63key_t knum, k63key_t krev, int c, struct asm_graph_t *g,
		struct k63_idhash_t *edict, struct k63_idhash_t *ndict)
{
	kmint_t e_pos, n_pos;
	gint_t e_id, e_id_rc, n_id;
	if (__k63_lt(knum, krev)) {
		n_pos = k63_idhash_get(ndict, knum);
		if (n_pos != IDHASH_END(ndict)) {
			n_id = IDHASH_ID(ndict, n_pos) * 2;
			e_id = asm_find_edge(g->edges, g->nodes[n_id].adj, g->ksize,
						g->nodes[n_id].deg, c);
			if (e_id >= 0) {
				atomic_add_and_fetch64(&(g->edges[e_id].count), 1);
				e_id_rc = g->edges[e_id].rc_id;
				atomic_add_and_fetch64(&(g->edges[e_id_rc].count), 1);
			}
		} else {
			e_pos = k63_idhash_get(edict, knum);
			if (e_pos == IDHASH_END(edict))
				return;
			e_id = IDHASH_ID(edict, e_pos);
			atomic_add_and_fetch64(&(g->edges[e_id].count), 1);
			e_id_rc = g->edges[e_id].rc_id;
			atomic_add_and_fetch64(&(g->edges[e_id_rc].count), 1);
		}
	} else {
		n_pos = k63_idhash_get(ndict, krev);
		if (n_pos != IDHASH_END(ndict)) {
			n_id = IDHASH_ID(ndict, n_pos) * 2 + 1;
			e_id = asm_find_edge(g->edges, g->nodes[n_id].adj, g->ksize,
						g->nodes[n_id].deg, c);
			if (e_id != -1) {
				atomic_add_and_fetch64(&(g->edges[e_id].count), 1);
				e_id_rc = g->edges[e_id].rc_id;
				atomic_add_and_fetch64(&(g->edges[e_id_rc].count), 1);
			}
		} else {
			e_pos = k63_idhash_get(edict, krev);
			if (e_pos == IDHASH_END(edict))
				return;
			e_id = IDHASH_ID(edict, e_pos);
			atomic_add_and_fetch64(&(g->edges[e_id].count), 1);
			e_id_rc = g->edges[e_id].rc_id;
			atomic_add_and_fetch64(&(g->edges[e_id_rc].count), 1);
		}
	}
}

static void count_edge_read(struct read_t *r, int ksize, struct asm_graph_t *graph,
			struct k63_idhash_t *edict, struct k63_idhash_t *ndict)
{
	int i, last, ci, ck, len, lmc, kedge;
	char *seq;
	len = r->len;
	seq = r->seq;

	k63key_t knum, krev, pknum, pkrev, kmask;
	kmask.bin[0] = (uint64_t)-1;
	kmask.bin[1] = (1ull << ((ksize << 1) - 64)) - 1;
	knum = krev = (k63key_t){{0ull, 0ull}};
	last = 0;
	lmc = (ksize - 1) << 1;
	kedge = ksize + 1;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		__k63_lshift2(knum); __k63_and(knum, kmask);
		__k63_rshift2(krev);
		// knum = (knum << 2) & kmask;
		// krev = krev >> 2;
		if (ci < 4) {
			knum.bin[0] |= ci;
			krev.bin[1] |= (uint64_t)(ci ^ 3) << (lmc - 64);
			// knum |= ci;
			// krev |= (k31key_t)(ci ^ 3) << lmc;
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

static void *count_edge_worker(void *data)
{
	struct edgecount_bundle_t *bundle = (struct edgecount_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct asm_graph_t *graph = bundle->graph;
	struct k63_idhash_t *edict, *ndict;
	edict = bundle->edict;
	ndict = bundle->ndict;

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

static void k63_count_edge(struct opt_count_t *opt, struct asm_graph_t *graph,
		struct k63_idhash_t *edict, struct k63_idhash_t *ndict)
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
		worker_bundles[i].ksize = graph->ksize;
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
