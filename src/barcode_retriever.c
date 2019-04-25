#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "fastq_producer.h"
#include "io_utils.h"
#include "k31hash.h"
#include "k63hash.h"
#include "k31_count.h"
#include "k63_count.h"
#include "khash.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

struct edge_coor_t {
	gint_t idx;
	gint_t pos;
	struct edge_coor_t *next;
};

#define __not_null(x) ((x).idx != -1)
#define MIN_HEAD_LEN 5

KHASH_INIT(k63_dict, k63key_t, struct edge_coor_t *, 1, __hash_k63, __k63_equal)

KHASH_INIT(k31_dict, k31key_t, struct edge_coor_t *, 1, __hash_k31, __k31_equal)

struct bctrie_bundle_t {
	struct dqueue_t *q;
	struct asm_graph_t *graph;
	void *dict;
	int64_t *n_reads;
};

void init_barcode_map(struct asm_graph_t *g, int buck_len);
void k31_build_naive_index(struct asm_graph_t *g, struct opt_build_t *opt,
						khash_t(k31_dict) *dict);
void k63_build_naive_index(struct asm_graph_t *g, struct opt_build_t *opt,
						khash_t(k63_dict) *dict);
static void k31_retrieve_barcode(struct asm_graph_t *g, struct opt_build_t *opt,
					khash_t(k31_dict) *dict);
static void k63_retrieve_barcode(struct asm_graph_t *g, struct opt_build_t *opt,
					khash_t(k63_dict) *dict);

gint_t count_bc(struct barcode_hash_t *t)
{
	gint_t ret = 0;
	kmint_t i;
	for (i = 0; i < t->size; ++i) {
		if (t->keys[i] == (uint64_t)-1)
			continue;
		ret += t->cnts[i] > 30;
	}
	return ret;
}

gint_t count_shared_bc(struct barcode_hash_t *t1, struct barcode_hash_t *t2)
{
	gint_t ret = 0;
	gint_t i, k;
	for (i = 0; i < t1->size; ++i) {
		if (t1->keys[i] == (uint64_t)-1)
			continue;
		k = barcode_hash_get(t2, t1->keys[i]);
		if (k != BARCODE_HASH_END(t2))
			ret += __min(t1->cnts[i], t2->cnts[k]) > 30;
	}
	return ret;
}

static inline void desc_sort(int *b, int *e)
{
	int *i, *j, t;
	for (i = b + 1; i < e; ++i) {
		if (*i > *(i - 1)) {
			t = *i;
			for (j = i; j > b && t > *(j - 1); --j)
				*j = *(j - 1);
			*j = t;
		}
	}
}

static inline int get_n90(int *a, int n)
{
	int s, sum, i;
	s = sum = 0;
	for (i = 0; i < n; ++i)
		sum += a[i];
	for (i = 0; i < n && s * 100 < sum * 90; ++i)
		s += a[i];
	return i;
}

void barcode_bin_profiling(struct asm_graph_t *g, gint_t *bx_bin)
{
	for (int i = 0; i < g->n_e; ++i){
		gint_t len1 = get_edge_len(g->edges + i);
		bx_bin[i] = (len1 + g->bin_size - 1) / g->bin_size;
	}
}

double get_barcode_ratio_small(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	gint_t len1, len2, n1, n2, i, k;
	len1 = get_edge_len(g->edges + e1);
	len2 = get_edge_len(g->edges + e2);
	if (len1 < MIN_FRAG_LEN || len2 < MIN_FRAG_LEN)
		return -1.0;
	n1 = (len1 + g->bin_size - 1) / g->bin_size;
	n2 = (len2 + g->bin_size - 1) / g->bin_size;
	n1 = __min(n1, 10);
	n2 = __min(n2, 10);
	gint_t *s1, *s2;
	s1 = alloca(n1 * sizeof(gint_t));
	s2 = alloca(n2 * sizeof(gint_t));
	for (i = 0; i < n1; ++i)
		s1[i] = count_bc(g->edges[e1].bucks + i);
	for (k = 0; k < n2; ++k)
		s2[k] = count_bc(g->edges[e2].bucks + k);
	double s = 0;
	gint_t cnt = 0;
	for (i = 0; i < n1; ++i) {
		if (s1[i] < MIN_UNIQUE_BARCODE)
			continue;
		for (k = 0; k < n2; ++k) {
			if (s2[k] < MIN_UNIQUE_BARCODE)
				continue;
			gint_t ret = count_shared_bc(g->edges[e1].bucks + i,
							g->edges[e2].bucks + k);
			s += ret * 1.0 / (s1[i] + s2[k]);
			++cnt;
		}
	}
	if (cnt == 0)
		return -1.0;
	return (s / cnt);
}

double get_barcode_ratio(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	gint_t len1, len2, n1, n2, i, k;
	len1 = get_edge_len(g->edges + e1);
	len2 = get_edge_len(g->edges + e2);
	if (len1 < MIN_CONTIG_BARCODE || len2 < MIN_CONTIG_BARCODE)
		return -1.0;
	n1 = (len1 + g->bin_size / 10) / g->bin_size;
	n2 = (len2 + g->bin_size / 10) / g->bin_size;
	n1 = __min(n1, 4);
	n2 = __min(n2, 4);
	// if (n1 * n2 < 5)
	// 	return -1.0;
	gint_t *s1, *s2;
	s1 = alloca(n1 * sizeof(gint_t));
	s2 = alloca(n2 * sizeof(gint_t));
	for (i = 0; i < n1; ++i)
		s1[i] = count_bc(g->edges[e1].bucks + i);
	for (k = 0; k < n2; ++k)
		s2[k] = count_bc(g->edges[e2].bucks + k);
	double s = 0;
	gint_t cnt = 0;
	for (i = 1; i < n1; ++i) {
		if (s1[i] < MIN_UNIQUE_BARCODE)
			continue;
		for (k = 1; k < n2; ++k) {
			if (s2[k] < MIN_UNIQUE_BARCODE)
				continue;
			gint_t ret = count_shared_bc(g->edges[e1].bucks + i,
						g->edges[e2].bucks + k);
			s += ret * 1.0 / (s1[i] + s2[k]);
			++cnt;
		}
	}
	if (cnt < 4)
		return -1.0;
	return (s  / cnt);
}

int test_edge_barcode(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	gint_t h, n1, n2, len1, len2;
	len1 = get_edge_len(g->edges + e1);
	len2 = get_edge_len(g->edges + e2);
	if (len1 < 1500 || len2 < 1500)
		return -1;
	n1 = (len1 + g->bin_size - 1) / g->bin_size;
	n2 = (len2 + g->bin_size - 1) / g->bin_size;
	h = __min(n1, 20) * __min(n2, 20);
	/* length is too short for barcoding */
	if (h < 5)
		return -1;
	int *s = calloc(n1 * n2, sizeof(int));
	gint_t i, k;
	for (i = 0; i < __min(n1, 20); ++i) {
		for (k = 0; k < __min(n2, 20); ++k) {
			s[i * __min(n2, 20) + k] = count_shared_bc(g->edges[e1].bucks + i,
					g->edges[e2].bucks + k);
		}
	}
	desc_sort(s, s + h);
	// for (i = 0; i < h; ++i)
	// 	fprintf(stdout, "%ld\n", s[i]);
	k = get_n90(s, h);
	// fprintf(stdout, "h = %ld; k = %ld; s[k  -1] = %d\n", h, k, s[k - 1]);
	//return k * 2 > h && s[k - 1] >= 5 ? 1 : 0;
	if ( k * 2 > h && s[k - 1] >= 5){
		__VERBOSE("N90 Pos: %ld\n" , k);
		__VERBOSE("Value at N90: %d\n" , s[k - 1]);
	}
	return k * 2 > h && s[k - 1] >= 5 ? 1 : 0;
}

int test_edge_barcode2(struct asm_graph_t *g, gint_t e1, gint_t e2, 
			gint_t *bx_bin, uint32_t *score)
{
	gint_t h, len1, len2, n1, n2;
	len1 = get_edge_len(g->edges + e1);
	len2 = get_edge_len(g->edges + e2);
	n1 = bx_bin[e1];
	n2 = bx_bin[e2];
	if (len1 < 1500 || len2 < 1500)
		return -1;
	h = __min(n1, MIN_HEAD_LEN) * __min(n2, MIN_HEAD_LEN);
	/* length is too short for barcoding */
	if (h < 5)
		return -1;
	int *s = calloc(n1 * n2, sizeof(int));
	uint32_t sum = 0;
	gint_t i, k;
	for (i = 0; i < __min(n1, MIN_HEAD_LEN); ++i) {
		for (k = 0; k < __min(n2, MIN_HEAD_LEN); ++k) {
			s[i * __min(n2, MIN_HEAD_LEN) + k] = count_shared_bc(g->edges[e1].bucks + i,
					g->edges[e2].bucks + k);
			sum += s[i * __min(n2, MIN_HEAD_LEN) + k];
		}
	}
	desc_sort(s, s + h);
	k = get_n90(s, h);
	if (s[k - 1] >= 7){
		*score = sum;
	}
	return s[k - 1] >= 7 ? 1 : 0;
}

void print_test_barcode_edge(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	assert(e1 < g->n_e && e2 < g->n_e);
	fprintf(stdout, "Print table (%ld <-> %ld)\n", e1, e2);
	gint_t n1, n2, len1, len2, i, k;
	len1 = get_edge_len(g->edges + e1);
	len2 = get_edge_len(g->edges + e2);
	if (len1 < 2000 || len2 < 2000)
		return;
	n1 = (len1 + g->bin_size / 10) / g->bin_size;
	n2 = (len2 + g->bin_size / 10) / g->bin_size;
	n1 = __min(n1, 4);
	n2 = __min(n2, 4);
	fprintf(stdout, "bin e1: ");
	for (i = 0; i < n1; ++i) {
		fprintf(stdout, i + 1 == n1 ? "%ld\n" : "%ld ",
			count_bc(g->edges[e1].bucks + i));
	}
	fprintf(stdout, "bin e2: ");
	for (i = 0; i < n2; ++i) {
		fprintf(stdout, i + 1 == n2 ? "%ld\n" : "%ld ",
			count_bc(g->edges[e2].bucks + i));
	}
	for (i = 0; i < n1; ++i) {
		for (k = 0; k < n2; ++k) {
			gint_t s = count_shared_bc(g->edges[e1].bucks + i,
					g->edges[e2].bucks + k);
			fprintf(stdout, k + 1 == n2 ? "%ld\n" : "%ld,", s);
		}
	}
	fprintf(stdout, "\n");
}

void print_test_barcode_edge2(struct asm_graph_t *g, gint_t e1, gint_t e2, 
				gint_t *bx_bin)
{
	assert(e1 < g->n_e && e2 < g->n_e);
	fprintf(stdout, "Print table (%ld <-> %ld)\n", e1, e2);
	gint_t n1, n2, len1, len2;
	len1 = get_edge_len(g->edges + e1);
	len2 = get_edge_len(g->edges + e2);
	n1 = bx_bin[e1];
	n2 = bx_bin[e2];
	gint_t i, k;
	for (i = 0; i < n1; ++i) {
		for (k = 0; k < n2; ++k) {
			gint_t s = count_shared_bc(g->edges[e1].bucks + i,
					g->edges[e2].bucks + k);
			fprintf(stdout, k + 1 == n2 ? "%ld\n" : "%ld,", s);
		}
	}
}
void construct_barcode_map(struct asm_graph_t *g, struct opt_build_t *opt)
{
	init_barcode_map(g, opt->split_len);
	if (g->ksize < 32) {
		khash_t(k31_dict) *edict = kh_init(k31_dict);
		k31_build_naive_index(g, opt, edict);
		k31_retrieve_barcode(g, opt, edict);
		kh_destroy(k31_dict, edict);
	} else if (g->ksize >= 32 && g->ksize < 64) {
		khash_t(k63_dict) *edict = kh_init(k63_dict);
		k63_build_naive_index(g, opt, edict);
		k63_retrieve_barcode(g, opt, edict);
		kh_destroy(k63_dict, edict);
	}
}

void init_barcode_map(struct asm_graph_t *g, int buck_len)
{
	gint_t i, e;
	g->bin_size = buck_len;
	for (e = 0; e < g->n_e; ++e) {
		// gint_t e_rc = g->edges[e].rc_id;
		// if (e > e_rc) {
		// 	/* at first, we imaginarily ceil the length of edges to
		// 	 *                        nearest multiple of split_len
		// 	 * then we can share the barcode hash between the two
		// 	 *                        reverse complemented sequences
		// 	 */
		// 	g->edges[e].bucks = g->edges[e_rc].bucks;
		// 	continue;
		// }
		pthread_mutex_init(&(g->edges[e].lock), NULL);
		uint32_t len = get_edge_len(g->edges + e);
		gint_t n_bucks = (len + buck_len - 1) / buck_len;
		g->edges[e].bucks = calloc(n_bucks, sizeof(struct barcode_hash_t));
		for (i = 0; i < n_bucks; ++i) {
			barcode_hash_init(g->edges[e].bucks + i, 4);
			// g->edges[e].bucks[i].lock = &(g->edges[e].lock);
		}
	}
}

void test_bucket(khash_t(k31_dict) *dict)
{
	khiter_t it;
	uint32_t sum = 0;
	for (it = kh_begin(dict); it != kh_end(dict); ++it) {
		if (!kh_exist(dict, it))
			continue;
		struct edge_coor_t *p;
		p = kh_value(dict, it);
		while (p != NULL) {
			++sum;
			p = p->next;
		}
	}
	fprintf(stderr, "Number of kmer: %u\n", kh_size(dict));
	fprintf(stderr, "Number of position: %u\n", sum);
	fprintf(stderr, "Mean positon/kmer: %f\n", sum * 1.0 / kh_size(dict));
}

void k31_build_naive_index(struct asm_graph_t *g, struct opt_build_t *opt,
						khash_t(k31_dict) *dict)
{
	gint_t e, i, k, last;
	k31key_t knum, krev, kmask;
	knum = krev = 0;
	kmask = ((k31key_t)1 << (g->ksize << 1)) - 1;
	int lmc = (g->ksize << 1) - 2;
	for (e = 0; e < g->n_e; ++e) {
		gint_t e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		uint32_t p = 0;
		gint_t slen = get_edge_len(g->edges + e);
		for (i = k = last = 0; i + 1 < g->edges[e].seq_len; ++i, ++k) {
			uint32_t c = __binseq_get(g->edges[e].seq, i);
			knum = ((knum << 2) & kmask) | c;
			krev = (krev >> 2) | ((k31key_t)(c ^ 3) << lmc);
			++last;
			if (last >= g->ksize && i >= g->ksize) {
				khiter_t it;
				int ret;
				if (knum <= krev)
					it = kh_put(k31_dict, dict, knum, &ret);
				else
					it = kh_put(k31_dict, dict, krev, &ret);
				if (ret == 1)
					kh_value(dict, it) = NULL;
				struct edge_coor_t *coor = calloc(1, sizeof(struct edge_coor_t));
				coor->idx = e;
				coor->pos = k;
				coor->next = kh_value(dict, it);
				kh_value(dict, it) = coor;
				coor = calloc(1, sizeof(struct edge_coor_t));
				coor->idx = e_rc;
				coor->pos = slen - k - 1;
				coor->next = kh_value(dict, it);
				kh_value(dict, it) = coor;
				// if (ret == 1) {
				// 	kh_value(dict, it) = (struct edge_coor_t){e, k};
				// } else
				// 	kh_value(dict, it) = (struct edge_coor_t){(gint_t)-1, (gint_t)-1};
			}
			if (p < g->edges[e].n_holes &&
				i == g->edges[e].p_holes[p]) {
				last = 0;
				k += g->edges[e].l_holes[p];
				++p;
			}
		}
	}
	test_bucket(dict);
}

void k63_build_naive_index(struct asm_graph_t *g, struct opt_build_t *opt,
							khash_t(k63_dict) *dict)
{
	gint_t e, i, k, last;
	k63key_t knum, krev, kmask;
	kmask.bin[0] = (uint64_t)-1;
	kmask.bin[1] = (1ull << ((g->ksize << 1) - 64)) - 1;
	knum = krev = (k63key_t){{0ull, 0ull}};
	int lmc = (g->ksize - 1) << 1;
	for (e = 0; e < g->n_e; ++e) {
		gint_t e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		uint32_t p = 0;
		gint_t slen = get_edge_len(g->edges + e);
		for (i = k = last = 0; i < g->edges[e].seq_len; ++i, ++k) {
			uint32_t c = __binseq_get(g->edges[e].seq, i);
			__k63_lshift2(knum); __k63_and(knum, kmask);
			knum.bin[0] |= c;
			__k63_rshift2(krev);
			krev.bin[1] |= (uint64_t)(c ^ 3) << (lmc - 64);
			++last;
			if (last >= g->ksize) {
				khiter_t it;
				int ret;
				if (__k63_lt(knum, krev))
					it = kh_put(k63_dict, dict, knum, &ret);
				else
					it = kh_put(k63_dict, dict, krev, &ret);
				if (ret == 1)
					kh_value(dict, it) = NULL;
				struct edge_coor_t *coor = calloc(1, sizeof(struct edge_coor_t));
				coor->idx = e;
				coor->pos = k;
				coor->next = kh_value(dict, it);
				kh_value(dict, it) = coor;
				coor = calloc(1, sizeof(struct edge_coor_t));
				coor->idx = e_rc;
				coor->pos = slen - k - 1;
				coor->next = kh_value(dict, it);
				kh_value(dict, it) = coor;
				// if (ret == 1)
				// 	kh_value(dict, it) = (struct edge_coor_t){e, k};
				// else
				// 	kh_value(dict, it) = (struct edge_coor_t){(gint_t)-1, (gint_t)-1};
			}
			if (p < g->edges[e].n_holes &&
				i == g->edges[e].p_holes[p]) {
				last = 0;
				k += g->edges[e].l_holes[p];
				++p;
			}
		}
	}
}

static inline uint64_t convert_barcode_info(struct read_t *r)
{
	char *s = r->info;
	int i, k, len = 0;
	uint64_t ret = 0;
	for (i = 0; s[i]; ++i) {
		if (strncmp(s + i, "BX:Z:", 5) == 0) {
			for (k = i + 5; s[k] && !__is_sep(s[k]); ++k) {
				ret = ret * 5 + nt4_table[(int)s[k]];
				++len;
			}
		}
	}
	assert(len == 18);
	return ret;
}

static void k31_barcode_retrieve_read(struct read_t *r, struct asm_graph_t *g,
					khash_t(k31_dict) *dict)
{
	uint64_t barcode;
	barcode = convert_barcode_info(r);
	int i, last, ci, len, lmc, ksize;
	char *seq;
	len = r->len;
	seq = r->seq;
	ksize = g->ksize;

	k31key_t knum, krev, kmask;
	kmask = ((k31key_t)1 << (ksize << 1)) - 1;
	knum = krev = 0;
	last = 0;
	lmc = (ksize - 1) << 1;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		knum = (knum << 2) & kmask;
		krev = krev >> 2;
		if (ci < 4) {
			knum |= ci;
			krev |= (k31key_t)(ci ^ 3) << lmc;
			++last;
		} else {
			last = 0;
		}
		if (last >= ksize) {
			khiter_t k;
			if (knum < krev) {
				k = kh_get(k31_dict, dict, knum);
			} else {
				k = kh_get(k31_dict, dict, krev);
			}
			if (k != kh_end(dict)) {
				struct edge_coor_t *coor;
				gint_t e, pos, bin_pos;
				coor = kh_value(dict, k);
				while (coor) {
					e = coor->idx;
					pos = coor->pos;
					bin_pos = pos / g->bin_size;
					pthread_mutex_lock(&(g->edges[e].lock));
					barcode_hash_inc_count(g->edges[e].bucks + bin_pos, barcode);
					pthread_mutex_unlock(&(g->edges[e].lock));
					coor = coor->next;
				}
			}
		}
	}
}

static void *k31_barcode_retriever(void *data)
{
	struct bctrie_bundle_t *bundle = (struct bctrie_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct asm_graph_t *graph = bundle->graph;
	khash_t(k31_dict) *dict = (khash_t(k31_dict) *)(bundle->dict);

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	int64_t n_reads;
	int64_t *gcnt_reads;
	gcnt_reads = bundle->n_reads;

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
			k31_barcode_retrieve_read(&read1, graph, dict);
			k31_barcode_retrieve_read(&read2, graph, dict);

			if (rc1 == READ_END)
				break;
		}
		n_reads = atomic_add_and_fetch64(gcnt_reads, n_reads);
		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

static void k31_retrieve_barcode(struct asm_graph_t *g, struct opt_build_t *opt,
					khash_t(k31_dict) *dict)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	int i;

	struct producer_bundle_t *producer_bundles;
	producer_bundles = init_fastq_PE(opt->n_threads, opt->n_files,
						opt->files_1, opt->files_2);

	struct bctrie_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct bctrie_bundle_t));

	int64_t n_reads;
	n_reads = 0;

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].graph = g;
		worker_bundles[i].dict = dict;
		worker_bundles[i].n_reads = &n_reads;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, k31_barcode_retriever,
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

static void k63_barcode_retrieve_read(struct read_t *r, struct asm_graph_t *g,
					khash_t(k63_dict) *dict)
{
	uint64_t barcode;
	barcode = convert_barcode_info(r);
	int i, last, ci, len, lmc, ksize;
	char *seq;
	len = r->len;
	seq = r->seq;
	ksize = g->ksize;

	k63key_t knum, krev, kmask;
	kmask.bin[0] = (uint64_t)-1;
	kmask.bin[1] = (1ull << ((g->ksize << 1) - 64)) - 1;
	knum = krev = (k63key_t){{0ull, 0ull}};
	lmc = (g->ksize - 1) << 1;
	last = 0;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		__k63_lshift2(knum); __k63_and(knum, kmask);
		__k63_rshift2(krev);
		if (ci < 4) {
			knum.bin[0] |= ci;
			krev.bin[1] |= (uint64_t)(ci ^ 3) << (lmc - 64);
			++last;
		} else {
			last = 0;
		}
		if (last >= ksize) {
			khiter_t k;
			if (__k63_lt(knum, krev)) {
				k = kh_get(k63_dict, dict, knum);
			} else {
				k = kh_get(k63_dict, dict, krev);
			}
			if (k != kh_end(dict)) {
				struct edge_coor_t *coor;
				gint_t e, pos, bin_pos;
				coor = kh_value(dict, k);
				while (coor) {
					e = coor->idx;
					pos = coor->pos;
					bin_pos = pos / g->bin_size;
					pthread_mutex_lock(&(g->edges[e].lock));
					barcode_hash_inc_count(g->edges[e].bucks + bin_pos, barcode);
					pthread_mutex_unlock(&(g->edges[e].lock));
					coor = coor->next;
				}
			}
		}
	}
}

static void *k63_barcode_retriever(void *data)
{
	struct bctrie_bundle_t *bundle = (struct bctrie_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct asm_graph_t *graph = bundle->graph;
	khash_t(k63_dict) *dict = (khash_t(k63_dict) *)(bundle->dict);

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	int64_t n_reads;
	int64_t *gcnt_reads;
	gcnt_reads = bundle->n_reads;

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
			k63_barcode_retrieve_read(&read1, graph, dict);
			k63_barcode_retrieve_read(&read2, graph, dict);

			if (rc1 == READ_END)
				break;
		}
		n_reads = atomic_add_and_fetch64(gcnt_reads, n_reads);
		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

static void k63_retrieve_barcode(struct asm_graph_t *g, struct opt_build_t *opt,
					khash_t(k63_dict) *dict)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	int i;

	struct producer_bundle_t *producer_bundles;
	// producer_bundles = init_fastq_PE(opt);
	producer_bundles = init_fastq_PE(opt->n_threads, opt->n_files,
						opt->files_1, opt->files_2);

	struct bctrie_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct bctrie_bundle_t));

	int64_t n_reads;
	n_reads = 0;

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].graph = g;
		worker_bundles[i].dict = dict;
		worker_bundles[i].n_reads = &n_reads;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, k63_barcode_retriever,
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
