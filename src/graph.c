
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "attribute.h"
#include "dqueue.h"
#include "get_buffer.h"
#include "graph.h"
#include "io_utils.h"
#include "kmer_count.h"
#include "kmhash.h"
#include "kseq.h"
#include "utils.h"
#include "verbose.h"

static void dump_seq(uint64_t num, char *seq, int len)
{
	seq[len] = '\0';
	while (len) {
		seq[--len] = nt4_char[num & 3];
		num >>= 2;
	}
}

#define __get_revc_num(y, x, l, mask)					       \
(	(x) = (y) << (64 - ((l) << 1)),					       \
	(x) = (((x) & 0xffffffff00000000ull) >> 32) | (((x) & 0x00000000ffffffffull) << 32), \
	(x) = (((x) & 0xffff0000ffff0000ull) >> 16) | (((x) & 0x0000ffff0000ffffull) << 16), \
	(x) = (((x) & 0xff00ff00ff00ff00ull) >>  8) | (((x) & 0x00ff00ff00ff00ffull) <<  8), \
	(x) = (((x) & 0xf0f0f0f0f0f0f0f0ull) >>  4) | (((x) & 0x0f0f0f0f0f0f0f0full) <<  4), \
	(x) = (((x) & 0xccccccccccccccccull) >>  2) | (((x) & 0x3333333333333333ull) <<  2), \
	(x) ^= 0xffffffffffffffffull, (x) &= (mask))

#define __off_bit(a, i) ((a)[(i) >> 5] &= (~(1 << ((i) & 31))))
#define __on_bit(a, i) ((a)[(i) >> 5] |= (1 << ((i) & 31)))
#define __get_bit(a, i) (1 & ((a)[(i) >> 5] >> ((i) & 31)))

#define __degree(e) ((!!(e)[0]) + (!!(e)[1]) + (!!(e)[2]) + (!!(e)[3]))

/* We assume that there is only one out going edge */
#define __only_edge(e) ((!!(e)[1]) + (!!(e)[2]) * 2 + (!!(e)[3]) * 3)

void reduce_graph(struct opt_count_t *opt, khash_t(kvert) *h, int16_t *e)
{
	struct graph_t *ret_g, g;
	ret_g = calloc(1, sizeof(struct graph_t));

	int m_v, m_e, mchain, lchain, old_lchain, k, uk, c;
	int ksize, node_id, u_node, v_node;
	uint64_t kmask, idx, ridx, u_idx, v_idx, u_ridx, v_ridx, tmp;
	ksize = opt->kmer_size;
	kmask = (1ull << (ksize << 1)) - 1;

	uint64_t *chain;
	int16_t *adj, *v_radj;
	g.n_k = kh_size(h);
	g.n_v = g.n_e = 0;
	g.kmer_chain_id = malloc(g.n_k * sizeof(int));
	g.chain_kmer = malloc(g.n_k * sizeof(uint64_t));
	g.n_k = 0;
	m_v = 0x1000;
	m_e = 0x1000;

	mchain = 128;
	chain = malloc(mchain * sizeof(uint64_t));

	g.kmer_count = malloc(m_v * sizeof(int));
	g.chain_head = malloc(m_v * sizeof(int));

	uint32_t *visited;
	visited = calloc((g.n_k + 31) / 32, sizeof(uint32_t));

	khint_t i, ik;
	for (i = kh_begin(h); i != kh_end(h); ++i) {
		if (!kh_exist(h, i))
			continue;
		idx = kh_key(h, i);
		node_id = kh_value(h, i).idx;
		if (__get_bit(visited, node_id))
			continue;
		__get_revc_num(idx, ridx, ksize, kmask);
		assert(idx <= ridx);

		if (g.n_v + 1 > m_v) {
			m_v <<= 1;
			g.kmer_count = realloc(g.kmer_count, (m_v + 1) * sizeof(int));
			g.chain_head = realloc(g.chain_head, (m_v + 1) * sizeof(int));
		}

		chain = g.chain_kmer + g.chain_head[g.n_v];

		chain[0] = idx;
		lchain = 1;
		g.kmer_chain_id[node_id] = g.n_v;
		g.kmer_count[g.n_v] = kh_value(h, i).cnt;
		__on_bit(visited, node_id);

		// test the forward going
		u_node = node_id;
		u_idx = idx;
		u_ridx = ridx;

		while (1) {
			adj = e + (u_node * 8 + 4 * (u_idx > u_ridx));
			if (__degree(adj) == 1) {
				c = __only_edge(adj);
				v_idx = ((u_idx << 2) & kmask) | c;
				v_ridx = (u_ridx >> 2) | ((uint64_t)(c ^ 3) << ((ksize << 1) - 2));
				if (v_idx < v_ridx) {
					ik = kh_get(kvert, h, v_idx);
					assert(ik != kh_end(h));
					v_node = kh_value(h, ik).idx;
					v_radj = e + (v_node * 8 + 4);
				} else {
					ik = kh_get(kvert, h, v_ridx);
					assert(ik != kh_end(h));
					v_node = kh_value(h, ik).idx;
					v_radj = e + (v_node * 8);
				}
				// Check if node v is belonged to another chain
				if (__get_bit(visited, v_node))
					break;
				// edge is 1-1
				if (__degree(v_radj) != 1)
					break;
				assert(u_ridx == (((v_ridx << 2) & kmask) | __only_edge(v_radj)));

				chain[lchain++] = v_idx;
				g.kmer_chain_id[v_node] = g.n_v;
				g.kmer_count[g.n_v] += kh_value(h, ik).cnt;
				__on_bit(visited, v_node);
				u_node = v_node;
				u_idx = v_idx;
				u_ridx = v_ridx;
			} else {
				break;
			}
		}

		old_lchain = lchain;
		u_node = node_id;
		u_idx = ridx;
		u_ridx = idx;

		// test the reverse edge
		while (1) {
			adj = e + (u_node * 8 + 4 * (u_idx > u_ridx));
			if (__degree(adj) == 1) {
				c = __only_edge(adj);
				v_idx = ((u_idx << 2) & kmask) | c;
				v_ridx = (u_ridx >> 2) | ((uint64_t)(c ^ 3) << ((ksize << 1) - 2));
				if (v_idx < v_ridx) {
					ik = kh_get(kvert, h, v_idx);
					assert(ik != kh_end(h));
					v_node = kh_value(h, ik).idx;
					v_radj = e + (v_node * 8 + 4);
				} else {
					ik = kh_get(kvert, h, v_ridx);
					assert(ik != kh_end(h));
					v_node = kh_value(h, ik).idx;
					v_radj = e + (v_node * 8);
				}
				// Check if node v is belonged to another chain
				if (__get_bit(visited, v_node))
					break;
				// edge is 1-1
				if (__degree(v_radj) != 1)
					break;
				assert(u_ridx == (((v_ridx << 2) & kmask) | __only_edge(v_radj)));

				chain[lchain++] = v_idx;
				g.kmer_chain_id[v_node] = g.n_v;
				g.kmer_count[g.n_v] += kh_value(h, ik).cnt;
				__on_bit(visited, v_node);
				u_node = v_node;
				u_idx = v_idx;
				u_ridx = v_ridx;
			} else {
				break;
			}
		}

		if (lchain > old_lchain) {
			// reverse the forward part
			for (k = 0; k < (old_lchain >> 1); ++k) {
				__get_revc_num(chain[k], tmp, ksize, kmask);
				__get_revc_num(chain[old_lchain - k - 1], chain[k], ksize, kmask);
				chain[old_lchain - k - 1] = tmp;
			}
			if (old_lchain & 1) {
				__get_revc_num(chain[(old_lchain + 1) >> 1], tmp, ksize, kmask);
				chain[(old_lchain + 1) >> 1] = tmp;
			}
		}

		g.chain_head[g.n_v + 1] = g.chain_head[g.n_v] + lchain;
		++g.n_v;
	}

	g.fhead = calloc(g.n_v + 1, sizeof(int));
	g.rhead = calloc(g.n_v + 1, sizeof(int));

	// Count edge
	for (k = 0; k < g.n_v; ++k) {
		// forward edge
		u_idx = g.chain_kmer[g.chain_head[k + 1]];
		__get_revc_num(u_idx, u_ridx, ksize, kmask);
		if (u_idx < u_ridx) {
			ik = kh_get(kvert, h, u_idx);
			assert(ik != kh_end(h));
			u_node = kh_value(h, ik).idx;
			adj = e + (u_node * 8);
		} else {
			ik = kh_get(kvert, h, u_ridx);
			assert(ik != kh_end(h));
			u_node = kh_value(h, ik).idx;
			adj = e + (u_node * 8 + 4);
		}
		g.fhead[k] = __degree(adj);

		// reverse edge
		u_ridx = g.chain_kmer[g.chain_head[k]];
		__get_revc_num(u_ridx, u_idx, ksize, kmask);
		if (u_idx < u_ridx) {
			ik = kh_get(kvert, h, u_idx);
			assert(ik != kh_end(h));
			u_node = kh_value(h, ik).idx;
			adj = e + (u_node * 8);
		} else {
			ik = kh_get(kvert, h, u_ridx);
			assert(ik != kh_end(h));
			u_node = kh_value(h, ik).idx;
			adj = e + (u_node * 8 + 4);
		}
		g.rhead[k] = __degree(adj);
	}

	for (k = 1; k <= g.n_v; ++k) {
		g.fhead[k] += g.fhead[k - 1];
		g.rhead[k] += g.rhead[k - 1];
	}

	g.fadj = malloc(g.fhead[g.n_v] * sizeof(int));
	g.radj = malloc(g.rhead[g.n_v] * sizeof(int));

	for (k = 0; k < g.n_v; ++k) {
		// forward edge
		u_idx = g.chain_kmer[g.chain_head[k + 1] - 1];
		__get_revc_num(u_idx, u_ridx, ksize, kmask);
		if (u_idx < u_ridx) {
			ik = kh_get(kvert, h, u_idx);
			assert(ik != kh_end(h));
			u_node = kh_value(h, ik).idx;
			adj = e + (u_node * 8);
		} else {
			ik = kh_get(kvert, h, u_ridx);
			assert(ik != kh_end(h));
			u_node = kh_value(h, ik).idx;
			adj = e + (u_node * 8 + 4);
		}
		for (c = 3; c >= 0; --c) {
			if (adj[c]) {
				v_idx = ((u_idx << 2) & kmask) | c;
				v_ridx = (u_ridx >> 2) | ((uint64_t)(c ^ 3) << ((ksize << 1) - 2));
				if (v_idx < v_ridx) {
					ik = kh_get(kvert, h, v_idx);
					assert(ik != kh_end(h));
					v_node = kh_value(h, ik).idx;
				} else {
					ik = kh_get(kvert, h, v_ridx);
					assert(ik != kh_end(h));
					v_node = kh_value(h, ik).idx;
				}
				uk = g.kmer_chain_id[v_node];
				if (v_idx == g.chain_kmer[g.chain_head[uk]])
					g.fadj[g.fhead[k]--] = uk + 1;
				else if (v_ridx == g.chain_kmer[g.chain_head[uk + 1] - 1])
					g.fadj[g.fhead[k]--] = -(uk + 1);
				else
					assert(0);
			}
		}

		// reverse edge
		u_ridx = g.chain_kmer[g.chain_head[k]];
		__get_revc_num(u_ridx, u_idx, ksize, kmask);
		if (u_idx < u_ridx) {
			ik = kh_get(kvert, h, u_idx);
			assert(ik != kh_end(h));
			u_node = kh_value(h, ik).idx;
			adj = e + (u_node * 8);
		} else {
			ik = kh_get(kvert, h, u_ridx);
			assert(ik != kh_end(h));
			u_node = kh_value(h, ik).idx;
			adj = e + (u_node * 8 + 4);
		}
		for (c = 3; c >= 0; --c) {
			if (adj[c]) {
				v_idx = ((u_idx << 2) & kmask) | c;
				v_ridx = (u_ridx >> 2) | ((uint64_t)(c ^ 3) << ((ksize << 1) - 2));
				if (v_idx < v_ridx) {
					ik = kh_get(kvert, h, v_idx);
					assert(ik != kh_end(h));
					v_node = kh_value(h, ik).idx;
				} else {
					ik = kh_get(kvert, h, v_ridx);
					assert(ik != kh_end(h));
					v_node = kh_value(h, ik).idx;
				}
				uk = g.kmer_chain_id[v_node];
				if (v_idx == g.chain_kmer[g.chain_head[uk]])
					g.radj[g.rhead[k]--] = uk + 1;
				else if (v_ridx == g.chain_kmer[g.chain_head[uk + 1] - 1])
					g.radj[g.rhead[k]--] = -(uk + 1);
				else
					assert(0);
			}
		}
	}

	int seq_len, m_len;
	char *seq;
	m_len = 128;
	seq = malloc(m_len);

	// dumping graph
	char dump_path[1024];
	strcpy(dump_path, opt->out_dir);
	strcat(dump_path, "/graph_reduced.gfa");
	FILE *fp = xfopen(dump_path, "w");

	for (k = 0; k < g.n_v; ++k) {
		seq_len = (g.chain_head[k + 1] - g.chain_head[k] - 1) + ksize;
		if (seq_len + 1 > m_len) {
			m_len = seq_len + 1;
			seq = realloc(seq, m_len + 1);
		}
		dump_seq(g.chain_kmer[g.chain_head[k]], seq, ksize);
		seq_len = ksize;
		for (c = g.chain_head[k] + 1; c < g.chain_head[k + 1]; ++c)
			seq[seq_len++] = nt4_char[g.chain_kmer[c] & 3];
		seq[seq_len] = '\0';
		fprintf(fp, "S\t%d\t%s\tKC:i:%d\n", k + 1, seq, g.kmer_count[k]);
	}

	for (k = 0; k < g.n_v; ++k) {
		for (c = g.fhead[k]; c < g.fhead[k + 1]; ++c) {
			fprintf(fp, "L\t%d\t+\t%d\t%c\t%dM\n",
				k + 1, g.fadj[c] > 0 ? g.fadj[c] : -g.fadj[c],
				g.fadj[c] > 0 ? '+' : '-', ksize - 1);
		}
		for (c = g.rhead[k]; c < g.rhead[k + 1]; ++c) {
			fprintf(fp, "L\t%d\t+\t%d\t%c\t%dM\n",
				k + 1, g.radj[c] > 0 ? g.radj[c] : -g.radj[c],
				g.radj[c] > 0 ? '+' : '-', ksize - 1);
		}
	}
	fclose(fp);
}


struct e_bundle_t {
	struct dqueue_t *q;
	khash_t(kvert) *h;
	int16_t *e;
	int ksize;
};

static struct pair_buffer_t *init_pair_buffer()
{
	struct pair_buffer_t *ret = malloc(sizeof(struct pair_buffer_t));
	ret->buf1 = malloc(BUF_SIZE + 1);
	ret->buf2 = malloc(BUF_SIZE + 1);
	return ret;
}

static void free_pair_buffer(struct pair_buffer_t *p)
{
	if (!p) return;
	free(p->buf1);
	free(p->buf2);
	free(p);
}

static struct dqueue_t *init_dqueue_PE(int cap)
{
	struct dqueue_t *ret = init_dqueue(cap);
	struct pair_buffer_t *p;
	int i;
	for (i = 0; i < cap; ++i) {
		p = init_pair_buffer();
		d_enqueue_out(ret, p);
	}
	return ret;
}

static void *producer_worker(void *data)
{
	struct producer_bundle_t *bundle = (struct producer_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct gb_pair_data *input_stream = bundle->stream;
	struct pair_buffer_t *own_buf = init_pair_buffer();
	struct pair_buffer_t *external_buf;
	int64_t offset;
	int64_t n_frag = 0;

	while ((offset = gb_get_pair(input_stream, &own_buf->buf1, &own_buf->buf2)) != -1) {
		own_buf->input_format = input_stream->type;
		external_buf = d_dequeue_out(q);
		d_enqueue_in(q, own_buf);
		own_buf = external_buf;
		++n_frag;
	}
	free_pair_buffer(own_buf);

	int cur;
	pthread_barrier_wait(bundle->barrier);
	while (1) {
		pthread_mutex_lock(bundle->lock);
		cur = *(bundle->n_consumer);
		if (*(bundle->n_consumer) > 0)
			--*(bundle->n_consumer);
		pthread_mutex_unlock(bundle->lock);
		if (cur == 0)
			break;
		external_buf = d_dequeue_out(q);
		free_pair_buffer(external_buf);
		d_enqueue_in(q, NULL);
	}

	pthread_exit(NULL);
}

void add_edge(struct read_t *r, struct e_bundle_t *bundle)
{
	khash_t(kvert) *h = bundle->h;
	int16_t *e = bundle->e;

	int i, k, last, ci, ck, len, lmc;
	char *seq;
	len = r->len;
	seq = r->seq;
	kmkey_t knum, krev, pknum, pkrev, kmask;
	khint_t ki, kk;


	k = bundle->ksize;
	kmask = ((kmkey_t)1 << (k << 1)) - 1;
	knum = krev = pknum = pkrev = 0;
	last = 0;
	lmc = (k - 1) << 1;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		knum = (knum << 2) & kmask;
		krev = (krev >> 2) & kmask;
		if (ci < 4) {
			knum |= ci;
			krev |= (kmkey_t)(ci ^ 3) << lmc;
			++last;
		} else {
			last = 0;
		}
		if (last >= k + 1) { // k + 1 for an edge
			ck = nt4_table[(int)seq[i - k]] ^ 3;
			// insert forward edge
			if (pknum < pkrev) {
				ki = kh_get(kvert, h, pknum);
			} else {
				ki = kh_get(kvert, h, pkrev);
				ci += 4;
			}
			if (knum < krev) {
				kk = kh_get(kvert, h, knum);
				ck += 4;
			} else {
				kk = kh_get(kvert, h, krev);
			}
			if (ki != kh_end(h) && kk != kh_end(h)) {
				__sync_fetch_and_add(e + (kh_value(h, ki).idx * 8 + ci), 1);
				__sync_fetch_and_add(e + (kh_value(h, kk).idx * 8 + ck), 1);
			}
		}
		pknum = knum;
		pkrev = krev;
	}

}

void *edge_worker(void *data)
{
	struct e_bundle_t *bundle = (struct e_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	// khash_t(kvert) *h = bundle->h;
	// int16_t *e = bundle->e;

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

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
		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, buf1, &pos1) :
				get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, buf2, &pos2) :
				get_read_from_fa(&read2, buf2, &pos2);


			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file\n");

			add_edge(&read1, bundle);
			add_edge(&read2, bundle);
			// count_kmer_read(&read1, bundle);
			// count_kmer_read(&read2, bundle);

			if (rc1 == READ_END)
				break;
		}
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);

}

int16_t *get_edges(struct opt_count_t *opt, khash_t(kvert) *h)
{
	int nvert = kh_size(h);
	__VERBOSE("Number of vertices: %d\n", nvert);
	int16_t *edges = calloc(nvert * 8, sizeof(int16_t));

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	struct dqueue_t *q;
	q = init_dqueue_PE(opt->n_threads * 2);
	int n_consumer;
	n_consumer = opt->n_threads;

	struct producer_bundle_t *producer_bundles;
	pthread_t *producer_threads;

	producer_bundles = malloc(opt->n_files * sizeof(struct producer_bundle_t));
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));

	pthread_mutex_t producer_lock;
	pthread_barrier_t producer_barrier;

	pthread_mutex_init(&producer_lock, NULL);
	pthread_barrier_init(&producer_barrier, NULL, opt->n_files);

	int i;
	for (i = 0; i < opt->n_files; ++i) {
		struct gb_pair_data *data = calloc(1, sizeof(struct gb_pair_data));
		gb_pair_init(data, opt->files_1[i], opt->files_2[i]);

		producer_bundles[i].n_consumer = &n_consumer;
		producer_bundles[i].stream = (void *)data;
		producer_bundles[i].q = q;
		producer_bundles[i].barrier = &producer_barrier;
		producer_bundles[i].lock = &producer_lock;
		pthread_create(producer_threads + i, &attr, producer_worker, producer_bundles + i);
	}

	struct e_bundle_t *worker_bundles;
	pthread_t *worker_threads;

	worker_bundles = malloc(opt->n_threads * sizeof(struct e_bundle_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = q;
		worker_bundles[i].h = h;
		worker_bundles[i].e = edges;
		// worker_bundles[i].n_threads = opt->n_threads;
		// worker_bundles[i].thread_no = i;
		// worker_bundles[i].barrier_hash = &barrier_hash;
		// worker_bundles[i].lock_count = &lock_count;
		// worker_bundles[i].lock_hash = h->locks + i;
		// worker_bundles[i].global_stat = &result;

		worker_bundles[i].ksize = opt->kmer_size;
		pthread_create(worker_threads + i, &attr, edge_worker, worker_bundles + i);
	}

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	return edges;
}
