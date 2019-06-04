#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "fastq_producer.h"
#include "io_utils.h"
#include "kmer.h"
#include "kmhash.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

#define __not_null(x) ((x).idx != -1)
#define MIN_HEAD_LEN 5
#define MAX_N_BUCK		6

// KHASH_INIT(k63_dict, k63key_t, struct edge_data_t, 1, __hash_k63, __k63_equal)

// KHASH_INIT(k31_dict, k31key_t, struct edge_data_t, 1, __hash_k31, __k31_equal)

struct bccount_bundle_t {
	struct asm_graph_t *g;
	struct dqueue_t *q;
	int64_t *n_reads;
	void *dict;
	void (*read_process_func)(struct read_t *, uint64_t, struct bccount_bundle_t *);
	uint64_t (*barcode_calculator)(struct read_t *, struct read_t *);
	int need_count;
};

uint64_t ust_get_barcode(struct read_t *r1, struct read_t *r2)
{
	char *s = r1->info;
	if (s == NULL)
		return (uint64_t)-1;
	int i, k, len = 0;
	uint64_t ret = 0;
	for (i = 0; s[i]; ++i) {
		if (strncmp(s + i, "BX:Z:", 5) == 0) {
			for (k = i + 5; s[k] && !__is_sep(s[k]); ++k) {
				ret = ret * 5 + nt4_table[(int)s[k]];
				++len;
			}
			break;
		}
	}
	if (len != 18)
		return (uint64_t)-1;
	return ret;
}

uint64_t x10_get_barcode(struct read_t *r1, struct read_t *r2)
{
	char *s = r1->seq;
	assert(r1->len >= 23); /* 16 bp barcode and 7 bp UMI */
	uint64_t ret = 0;
	int i;
	for (i = 0; i < 16; ++i)
		ret = ret * 5 + nt4_table[(int)s[i]];
	r1->seq += 23;
	return ret;
}

uint64_t (*barcode_calculators[])(struct read_t *, struct read_t *) = {ust_get_barcode, x10_get_barcode};

void init_barcode_map(struct asm_graph_t *g, int buck_len, int is_small);
void kmer_build_index(struct kmhash_t *h, struct asm_graph_t *g, int is_small);
// void k31_build_index(struct asm_graph_t *g, struct opt_proc_t *opt,
// 					khash_t(k31_dict) *dict, int is_small);
// void k63_build_index(struct asm_graph_t *g, struct opt_proc_t *opt,
// 					khash_t(k63_dict) *dict, int is_small);
// void k31_read_iterator(struct read_t *r, uint64_t barcode,
// 					struct bccount_bundle_t *bundle);
// void k63_read_iterator(struct read_t *r, uint64_t barcode,
// 					struct bccount_bundle_t *bundle);
void barcode_start_count(struct opt_proc_t *opt, struct bccount_bundle_t *ske);

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
			s += ret * 1.0 / (s1[i] + s2[k] - ret);
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
			s += ret * 1.0 / (s1[i] + s2[k] - ret);
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

// void debug_build_index(khash_t(k31_dict) *dict)
// {
// 	gint_t s = 0;
// 	khiter_t it;
// 	for (it = kh_begin(dict); it != kh_end(dict); ++it) {
// 		if (kh_exist(dict, it))
// 			s += kh_value(dict, it).n;
// 	}
// 	__VERBOSE("Number of indexed position: %ld\n", s);
// }

void construct_barcode_map_ust(struct opt_proc_t *opt, struct asm_graph_t *g,
					int is_small, int need_count)
{
	init_barcode_map(g, opt->split_len, is_small);
	struct bccount_bundle_t ske;
	ske.g = g;
	extern uint64_t (*barcode_calculators[])(struct read_t *, struct read_t *);
	ske.barcode_calculator = barcode_calculators[opt->lib_type];
	// ske.barcode_calculator = ust_get_barcode;
	ske.need_count = need_count;
	struct kmhash_t edict;
	kmhash_init(&edict, opt->hash_size, g->ksize + 1, KM_AUX_POS);
	kmer_build_index(&edict, g, is_small);
	ske.dict = &edict;
	barcode_start_count(opt, &ske);
	// if (g->ksize < 32) {
	// 	khash_t(k31_dict) *edict = kh_init(k31_dict);
	// 	k31_build_index(g, opt, edict, is_small);
	// 	ske.dict = edict;
	// 	ske.read_process_func = k31_read_iterator;
	// } else if (g->ksize >= 32 && g->ksize < 64) {
	// 	khash_t(k63_dict) *edict = kh_init(k63_dict);
	// 	k63_build_index(g, opt, edict, is_small);
	// 	ske.dict = edict;
	// 	ske.read_process_func = k63_read_iterator;
	// }
	// barcode_start_count(opt, &ske);
	// if (g->ksize < 32) {
	// 	khash_t(k31_dict) *edict = ske.dict;
	// 	khiter_t it;
	// 	for (it = kh_begin(edict); it != kh_end(edict); ++it) {
	// 		if (kh_exist(edict, it))
	// 			free(kh_value(edict, it).e);
	// 	}
	// 	kh_destroy(k31_dict, edict);
	// } else {
	// 	khash_t(k63_dict) *edict = ske.dict;
	// 	khiter_t it;
	// 	for (it = kh_begin(edict); it != kh_end(edict); ++it) {
	// 		if (kh_exist(edict, it))
	// 			free(kh_value(edict, it).e);
	// 	}
	// 	kh_destroy(k63_dict, edict);
	// }
}

void init_barcode_map(struct asm_graph_t *g, int buck_len, int is_small)
{
	gint_t i, e;
	g->bin_size = buck_len;
	for (e = 0; e < g->n_e; ++e) {
		pthread_mutex_init(&(g->edges[e].lock), NULL);
		uint32_t len = get_edge_len(g->edges + e);
		gint_t n_bucks = (len + buck_len - 1) / buck_len;
		if (is_small)
			n_bucks = __min(n_bucks, MAX_N_BUCK);
		g->edges[e].bucks = calloc(n_bucks, sizeof(struct barcode_hash_t));
		for (i = 0; i < n_bucks; ++i) {
			barcode_hash_init(g->edges[e].bucks + i, 4);
		}
	}
}

// void test_bucket(khash_t(k31_dict) *dict)
// {
// 	khiter_t it;
// 	uint32_t sum = 0;
// 	for (it = kh_begin(dict); it != kh_end(dict); ++it) {
// 		if (!kh_exist(dict, it))
// 			continue;
// 		struct edge_coor_t *p;
// 		p = kh_value(dict, it);
// 		while (p != NULL) {
// 			++sum;
// 			p = p->next;
// 		}
// 	}
// 	fprintf(stderr, "Number of kmer: %u\n", kh_size(dict));
// 	fprintf(stderr, "Number of position: %u\n", sum);
// 	fprintf(stderr, "Mean positon/kmer: %f\n", sum * 1.0 / kh_size(dict));
// }

void kmer_build_index(struct kmhash_t *h, struct asm_graph_t *g, int is_small)
{
	int ksize, word_size;
	ksize = g->ksize + 1;
	word_size = (ksize + 3) >> 2;
	uint8_t *knum, *krev;
	knum = alloca(word_size);
	krev = alloca(word_size);
	gint_t e, e_rc, i, k, last;
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		uint32_t p = 0;
		memset(knum, 0, word_size);
		memset(krev, 0, word_size);
		gint_t slen, seq_len, edge_len;
		edge_len = get_edge_len(g->edges + e);
		if (is_small)
			slen = __min(edge_len, g->bin_size * MAX_N_BUCK);
		else
			slen = edge_len;
		seq_len = g->edges[e].seq_len;
		for (i = k = last = 0; i < seq_len && k < slen; ++i, ++k) {
			uint32_t c = __binseq_get(g->edges[e].seq, i);
			km_shift_append(knum, ksize, word_size, c);
			km_shift_append_rv(krev, ksize, word_size, c ^ 3);
			++last;
			if (last >= ksize) {
				kmint_t it;
				if (km_cmp(knum, krev, word_size) <= 0)
					it = kmhash_put(h, knum);
				else
					it = kmhash_put(h, knum);
				struct edge_data_t *b = &KMHASH_POS(h, it);
				if ((b->n & (b->n + 1)) == 0) {
					if (b->n == 0)
						b->e = calloc(2, sizeof(struct edge_idx_t));
					else
						b->e = realloc(b->e, (b->n << 1) * sizeof(struct edge_idx_t));
				}
				b->e[b->n].idx = e;
				b->e[b->n].pos = k;
				++b->n;
				b->e[b->n].idx = e_rc;
				b->e[b->n].pos = edge_len - k - ksize;
				++b->n;
			}
			if (p < g->edges[e].n_holes && i == g->edges[e].p_holes[p]) {
				last = 0;
				k += g->edges[e].l_holes[p++];
			}
		}
	}
}

// void k31_build_index(struct asm_graph_t *g, struct opt_proc_t *opt,
// 					khash_t(k31_dict) *dict, int is_small)
// {
// 	gint_t e, i, k, last;
// 	k31key_t knum, krev, kmask;
// 	knum = krev = 0;
// 	kmask = ((k31key_t)1 << (g->ksize << 1)) - 1;
// 	int lmc = (g->ksize << 1) - 2;
// 	for (e = 0; e < g->n_e; ++e) {
// 		uint32_t p = 0;
// 		gint_t slen, seq_len;
// 		slen = get_edge_len(g->edges + e);
// 		if (is_small)
// 			slen = __min(slen, g->bin_size * MAX_N_BUCK);
// 		seq_len = g->edges[e].seq_len;
// 		for (i = k = last = 0; i < seq_len && k < slen; ++i, ++k) {
// 			uint32_t c = __binseq_get(g->edges[e].seq, i);
// 			knum = ((knum << 2) & kmask) | c;
// 			krev = (krev >> 2) | ((k31key_t)(c ^ 3) << lmc);
// 			++last;
// 			if (last >= g->ksize) {
// 				khiter_t it;
// 				int ret;
// 				if (knum <= krev)
// 					it = kh_put(k31_dict, dict, knum, &ret);
// 				else
// 					it = kh_put(k31_dict, dict, krev, &ret);
// 				struct edge_data_t *b = &kh_value(dict, it);
// 				if (ret == 1) {
// 					b->e = NULL;
// 					b->n = 0;
// 				}
// 				b->e = realloc(b->e, (b->n + 1) * sizeof(struct edge_idx_t));
// 				b->e[b->n].idx = e;
// 				b->e[b->n].pos = k;
// 				++b->n;
// 			}
// 			if (p < g->edges[e].n_holes &&
// 				i == g->edges[e].p_holes[p]) {
// 				last = 0;
// 				k += g->edges[e].l_holes[p];
// 				++p;
// 			}
// 		}
// 	}
// }

// void k63_build_index(struct asm_graph_t *g, struct opt_proc_t *opt,
// 					khash_t(k63_dict) *dict, int is_small)
// {
// 	gint_t e, i, k, last;
// 	k63key_t knum, krev, kmask;
// 	kmask.bin[0] = (uint64_t)-1;
// 	kmask.bin[1] = (1ull << ((g->ksize << 1) - 64)) - 1;
// 	knum = krev = (k63key_t){{0ull, 0ull}};
// 	int lmc = (g->ksize << 1) - 2;
// 	for (e = 0; e < g->n_e; ++e) {
// 		uint32_t p = 0;
// 		gint_t slen, seq_len;
// 		slen = get_edge_len(g->edges + e);
// 		if (is_small)
// 			slen = __min(slen, g->bin_size * MAX_N_BUCK);
// 		seq_len = g->edges[e].seq_len;
// 		for (i = k = last = 0; i < seq_len && k < slen; ++i, ++k) {
// 			uint32_t c = __binseq_get(g->edges[e].seq, i);
// 			__k63_lshift2(knum); __k63_and(knum, kmask);
// 			knum.bin[0] |= c;
// 			__k63_rshift2(krev);
// 			krev.bin[1] |= (uint64_t)(c ^ 3) << (lmc - 64);
// 			++last;
// 			if (last >= g->ksize) {
// 				khiter_t it;
// 				int ret;
// 				if (__k63_lt(knum, krev))
// 					it = kh_put(k63_dict, dict, knum, &ret);
// 				else
// 					it = kh_put(k63_dict, dict, krev, &ret);
// 				struct edge_data_t *b = &kh_value(dict, it);
// 				if (ret == 1) {
// 					b->e = NULL;
// 					b->n = 0;
// 				}
// 				b->e = realloc(b->e, (b->n + 1) * sizeof(struct edge_idx_t));
// 				b->e[b->n].idx = e;
// 				b->e[b->n].pos = k;
// 				++b->n;
// 			}
// 			if (p < g->edges[e].n_holes &&
// 				i == g->edges[e].p_holes[p]) {
// 				last = 0;
// 				k += g->edges[e].l_holes[p];
// 				++p;
// 			}
// 		}
// 	}
// }

void bcread_iterator(struct read_t *r, uint64_t barcode, struct bccount_bundle_t *bundle)
{
	struct asm_graph_t *g = bundle->g;
	struct kmhash_t *dict = bundle->dict;

	int ksize, word_size, i, k, last, len;
	ksize = g->ksize + 1;
	word_size = (ksize + 3) >> 2;
	uint8_t *knum, *krev;
	uint8_t ci;
	knum = alloca(word_size);
	krev = alloca(word_size);
	char *seq = r->seq;
	len = r->len;

	last = 0;
	for (i = 0; i < len; ++i) {
		ci = nt4_table[(int)seq[i]];
		if (ci < 4) {
			km_shift_append(knum, ksize, word_size, ci);
			km_shift_append_rv(krev, ksize, word_size, ci ^ 3);
			++last;
		} else {
			last = 0;
		}
		if (last >= ksize) {
			kmint_t k;
			if (km_cmp(knum, krev, word_size) <= 0)
				k = kmhash_get(dict, knum);
			else
				k = kmhash_get(dict, krev);
			if (k == KMHASH_END(dict))
				continue;
			struct edge_idx_t *p = KMHASH_POS(dict, k).e;
			gint_t n, e, bin, j;
			n = KMHASH_POS(dict, k).n;
			for (j = 0; j < n; ++j) {
				e = p[j].idx;
				bin = p[j].pos / g->bin_size;
				pthread_mutex_lock(&(g->edges[e].lock));
				barcode_hash_inc_count(g->edges[e].bucks + bin, barcode);
				pthread_mutex_unlock(&(g->edges[e].lock));
			}
			if (bundle->need_count) {
				gint_t prev_e = -1;
				for (j = 0; j < n; ++j) {
					e = p[j].idx;
					if (e != prev_e)
						atomic_add_and_fetch64(&(g->edges[e].count), 1);
					prev_e = e;
				}
			}
		}
	}
}

// void k31_read_iterator(struct read_t *r, uint64_t barcode,
// 						struct bccount_bundle_t *bundle)
// {
// 	struct asm_graph_t *g = bundle->g;
// 	khash_t(k31_dict) *dict = bundle->dict;
// 	int i, last, ci, len, lmc, ksize;
// 	char *seq;
// 	len = r->len;
// 	seq = r->seq;
// 	ksize = g->ksize;

// 	k31key_t knum, krev, kmask;
// 	kmask = ((k31key_t)1 << (ksize << 1)) - 1;
// 	knum = krev = 0;
// 	last = 0;
// 	lmc = (ksize - 1) << 1;
// 	for (i = 0; i < len; ++i) {
// 		ci = nt4_table[(int)seq[i]];
// 		knum = (knum << 2) & kmask;
// 		krev = krev >> 2;
// 		if (ci < 4) {
// 			knum |= ci;
// 			krev |= (k31key_t)(ci ^ 3) << lmc;
// 			++last;
// 		} else {
// 			last = 0;
// 		}
// 		if (last >= ksize) {
// 			khiter_t k;
// 			if (knum < krev) {
// 				k = kh_get(k31_dict, dict, knum);
// 			} else {
// 				k = kh_get(k31_dict, dict, krev);
// 			}
// 			if (k != kh_end(dict)) {
// 				struct edge_idx_t *p = kh_value(dict, k).e;
// 				gint_t n, e, bin, j;
// 				n = kh_value(dict, k).n;
// 				for (j = 0; j < n; ++j) {
// 					e = p[j].idx;
// 					bin = p[j].pos / g->bin_size;
// 					pthread_mutex_lock(&(g->edges[e].lock));
// 					barcode_hash_inc_count(g->edges[e].bucks + bin, barcode);
// 					pthread_mutex_unlock(&(g->edges[e].lock));
// 				}
// 				if (bundle->need_count) {
// 					gint_t prev_e = -1;
// 					for (j = 0; j < n; ++j) {
// 						e = p[j].idx;
// 						if (e != prev_e)
// 							atomic_add_and_fetch64(&(g->edges[e].count), 1);
// 						prev_e = e;
// 					}
// 				}
// 			}
// 		}
// 	}
// }

// void k63_read_iterator(struct read_t *r, uint64_t barcode,
// 						struct bccount_bundle_t *bundle)
// {
// 	struct asm_graph_t *g = bundle->g;
// 	khash_t(k63_dict) *dict = bundle->dict;
// 	int i, last, ci, len, lmc, ksize;
// 	char *seq;
// 	len = r->len;
// 	seq = r->seq;
// 	ksize = g->ksize;

// 	k63key_t knum, krev, kmask;
// 	kmask.bin[0] = (uint64_t)-1;
// 	kmask.bin[1] = (1ull << ((g->ksize << 1) - 64)) - 1;
// 	knum = krev = (k63key_t){{0ull, 0ull}};
// 	lmc = (g->ksize - 1) << 1;
// 	last = 0;
// 	for (i = 0; i < len; ++i) {
// 		ci = nt4_table[(int)seq[i]];
// 		__k63_lshift2(knum); __k63_and(knum, kmask);
// 		__k63_rshift2(krev);
// 		if (ci < 4) {
// 			knum.bin[0] |= ci;
// 			krev.bin[1] |= (uint64_t)(ci ^ 3) << (lmc - 64);
// 			++last;
// 		} else {
// 			last = 0;
// 		}
// 		if (last >= ksize) {
// 			khiter_t k;
// 			if (__k63_lt(knum, krev)) {
// 				k = kh_get(k63_dict, dict, knum);
// 			} else {
// 				k = kh_get(k63_dict, dict, krev);
// 			}
// 			if (k != kh_end(dict)) {
// 				struct edge_idx_t *p = kh_value(dict, k).e;
// 				gint_t n, e, bin, j;
// 				n = kh_value(dict, k).n;
// 				for (j = 0; j < n; ++j) {
// 					e = p[j].idx;
// 					bin = p[j].pos / g->bin_size;
// 					pthread_mutex_lock(&(g->edges[e].lock));
// 					barcode_hash_inc_count(g->edges[e].bucks + bin, barcode);
// 					pthread_mutex_unlock(&(g->edges[e].lock));
// 				}
// 				if (bundle->need_count) {
// 					gint_t prev_e = -1;
// 					for (j = 0; j < n; ++j) {
// 						e = p[j].idx;
// 						if (e != prev_e)
// 							atomic_add_and_fetch64(&(g->edges[e].count), 1);
// 						prev_e = e;
// 					}
// 				}
// 			}
// 		}
// 	}
// }

void *barcode_retriever(void *data)
{
	struct bccount_bundle_t *bundle = (struct bccount_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format;

	int64_t n_reads;
	int64_t *gcnt_reads;
	uint64_t barcode;
	gcnt_reads = bundle->n_reads;
	// void (*read_process)(struct read_t *, uint64_t, struct bccount_bundle_t *) = bundle->read_process_func;
	uint64_t (*barcode_calculator)(struct read_t *, struct read_t *) = bundle->barcode_calculator;

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
			barcode = barcode_calculator(&read1, &read2);
			if (barcode != (uint64_t)-1) {
				bcread_iterator(&read1, barcode, bundle);
				bcread_iterator(&read2, barcode, bundle);
			}

			if (rc1 == READ_END)
				break;
		}
		n_reads = atomic_add_and_fetch64(gcnt_reads, n_reads);
		__VERBOSE("\rNumber of process read:    %lld", (long long)n_reads);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

void barcode_start_count(struct opt_proc_t *opt, struct bccount_bundle_t *ske)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int i;
	struct producer_bundle_t *producer_bundles;
	producer_bundles = init_fastq_PE(opt->n_threads, opt->n_files,
						opt->files_1, opt->files_2);

	struct bccount_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct bccount_bundle_t));
	int64_t n_reads = 0;
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].n_reads = &n_reads;
		worker_bundles[i].g = ske->g;
		worker_bundles[i].dict = ske->dict;
		// worker_bundles[i].read_process_func = ske->read_process_func;
		worker_bundles[i].barcode_calculator = ske->barcode_calculator;
		worker_bundles[i].need_count = ske->need_count;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_PE_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, barcode_retriever,
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

