#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "atomic.h"
#include "KMC_reader.h"
#include "kmer.h"
#include "kmer_build.h"
#include "kmhash.h"
#include "utils.h"
#include "verbose.h"

#define __bin_degree4(e) (((e) & 1) + (((e) >> 1) & 1) + (((e) >> 2) & 1) + (((e) >> 3) & 1))

#define __bin_only4(e) ((((e) >> 1) & 1) * 1 + (((e) >> 2) & 1) * 2 + (((e) >> 3) & 1) * 3)

#define __kmerseq_get(seq, k) (((seq)[(k) >> 2] >> (((k) & 3) << 1)) & (uint32_t)3)

struct kmbuild_bundle_t {
	struct kmhash_t *h;
	int ksize;
	uint8_t *k1;
	uint8_t *k2;
	uint8_t *k1_rc;
	uint8_t *k2_rc;
};

void build_asm_graph_from_kmhash(int n_threads, int ksize,
				struct kmhash_t *h, struct asm_graph_t *g);
void assign_edge_kmer_count(struct opt_proc_t *opt, struct kmhash_t *h,
							struct asm_graph_t *g);
void build_edge_kmer_index(struct kmhash_t *h, struct asm_graph_t *g);
void build_edge_kmer_index_multi(struct opt_proc_t *opt, struct kmhash_t *h, struct asm_graph_t *g);

/* ------------------------------>
 *  T A C C A C T G G G A T T C A
 *  4 3 2 1 0 9 8 7 6 5 4 3 2 1 0
 * 987654321098765432109876543210
 *      3       2       1       0
 */

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

static inline int is_seq_rc_kmer(uint8_t *seq1, uint8_t *seq2, int l)
{
	uint8_t c1, c2;
	int i, k;
	for (i = 0 ; i < l; ++i) {
		k = l - i - 1;
		c1 = __kmerseq_get(seq1, i);
		c2 = __kmerseq_get(seq2, k);
		if (c1 != (c2 ^ 3))
			return 0;
	}
	return 1;
}

static inline void dump_kmer(uint8_t *kmer, int len, char *s)
{
	int i;
	for (i = 0; i < len; ++i)
		s[len - i - 1] = nt4_char[(uint32_t)((kmer[i >> 2] >> ((i & 3) << 1)) & 0x3)];
	s[len] = '\0';
}

static inline int check_kmer(uint8_t *kmer, int len, int l)
{
	int k = ((len - 1) & 3) - 1;
	return ((kmer[l - 1] & (~(((uint8_t)1 << (k << 1)) - 1))) == 0);
}

static inline int check_edge(struct kmhash_t *h, uint8_t *knum, int ksize, int word_size)
{
	kmint_t i, k;
	uint8_t *krev, *cur_knum, *cur_krev;
	krev = alloca(word_size);
	cur_knum = alloca(word_size);
	cur_krev = alloca(word_size);
	i = kmhash_get(h, knum);
	assert(i != KMHASH_END(h));
	km_get_rc(krev, knum, ksize, word_size);
	assert(check_kmer(krev, ksize, word_size));
	int c;
	for (c = 0; c < 4; ++c) {
		if ((KMHASH_ADJ(h, i) >> c) & 1) {
			memcpy(cur_knum, knum, word_size);
			memcpy(cur_krev, krev, word_size);
			km_shift_append(cur_knum, ksize, word_size, c);
			km_shift_append_rv(cur_krev, ksize, word_size, c ^ 3);
			assert(check_kmer(cur_knum, ksize, word_size));
			assert(check_kmer(cur_krev, ksize, word_size));
			if (km_cmp(cur_knum, cur_krev, word_size) <= 0)
				k = kmhash_get(h, cur_knum);
			else
				k = kmhash_get(h, cur_krev);
			if (k == KMHASH_END(h)) {
				__VERBOSE("sad 1");
				char *s = alloca(ksize + 1);
				__VERBOSE("c = %d\n", c);
				dump_kmer(knum, ksize, s);
				__VERBOSE("knum = %s\n", s);
				dump_kmer(krev, ksize, s);
				__VERBOSE("krev = %s\n", s);
				dump_kmer(cur_knum, ksize, s);
				__VERBOSE("nnum = %s\n", s);
				dump_kmer(cur_krev, ksize, s);
				__VERBOSE("nrev = %s\n", s);
				return 0;
			}
		}

		if ((KMHASH_ADJ(h, i) >> (c + 4)) & 1) {
			memcpy(cur_knum, krev, word_size);
			memcpy(cur_krev, knum, word_size);
			km_shift_append(cur_knum, ksize, word_size, c);
			km_shift_append_rv(cur_krev, ksize, word_size, c ^ 3);
			assert(check_kmer(cur_knum, ksize, word_size));
			assert(check_kmer(cur_krev, ksize, word_size));
			if (km_cmp(cur_knum, cur_krev, word_size) <= 0)
				k = kmhash_get(h, cur_knum);
			else
				k = kmhash_get(h, cur_krev);
			if (k == KMHASH_END(h)) {
				__VERBOSE("sad 2");
				char *s = alloca(ksize + 1);
				__VERBOSE("c = %d\n", c);
				dump_kmer(knum, ksize, s);
				__VERBOSE("knum = %s\n", s);
				dump_kmer(krev, ksize, s);
				__VERBOSE("krev = %s\n", s);
				dump_kmer(cur_knum, ksize, s);
				__VERBOSE("nnum = %s\n", s);
				dump_kmer(cur_krev, ksize, s);
				__VERBOSE("nrev = %s\n", s);
				return 0;
			}
		}
	}
	return 1;
}

void test_kmhash(struct kmhash_t *h, int ksize)
{
	kmint_t i, k;
	int c, word_size = (ksize + 3) >> 2;
	uint8_t *knum, *krev, *cur_knum, *cur_krev;
	krev = alloca(word_size);
	cur_knum = alloca(word_size);
	cur_krev = alloca(word_size);
	for (i = 0; i < h->size; ++i) {
		if (KM_IS_EMPTY(h, i))
			continue;
		knum = KMHASH_KEY(h, i);
		assert(check_kmer(knum, ksize, word_size));
		k = kmhash_get(h, knum);
		if (k != i) {
			__ERROR("Test fail k != i");
		}
		km_get_rc(krev, knum, ksize, word_size);
		assert(check_kmer(krev, ksize, word_size));
		for (c = 0; c < 4; ++c) {
			if ((KMHASH_ADJ(h, i) >> c) & 1) {
				memcpy(cur_knum, knum, word_size);
				memcpy(cur_krev, krev, word_size);
				km_shift_append(cur_knum, ksize, word_size, c);
				km_shift_append_rv(cur_krev, ksize, word_size, c ^ 3);
				assert(check_kmer(cur_knum, ksize, word_size));
				assert(check_kmer(cur_krev, ksize, word_size));
				if (km_cmp(cur_knum, cur_krev, word_size) <= 0)
					k = kmhash_get(h, cur_knum);
				else
					k = kmhash_get(h, cur_krev);
			}

			if ((KMHASH_ADJ(h, i) >> (c + 4)) & 1) {
				memcpy(cur_knum, krev, word_size);
				memcpy(cur_krev, knum, word_size);
				km_shift_append(cur_knum, ksize, word_size, c);
				km_shift_append_rv(cur_krev, ksize, word_size, c ^ 3);
				assert(check_kmer(cur_knum, ksize, word_size));
				assert(check_kmer(cur_krev, ksize, word_size));
				if (km_cmp(cur_knum, cur_krev, word_size) <= 0)
					k = kmhash_get(h, cur_knum);
				else
					k = kmhash_get(h, cur_krev);
				if (k == KMHASH_END(h)) {
					char *s = alloca(ksize + 1);
					__VERBOSE("c = %d\n", c);
					dump_kmer(knum, ksize, s);
					__VERBOSE("knum = %s\n", s);
					dump_kmer(krev, ksize, s);
					__VERBOSE("krev = %s\n", s);
					dump_kmer(cur_knum, ksize, s);
					__VERBOSE("nnum = %s\n", s);
					dump_kmer(cur_krev, ksize, s);
					__VERBOSE("nrev = %s\n", s);
				}
			}
		}
	}
}

void split_kmer_from_kedge_multi(int thread_no, uint8_t *kedge, uint32_t count, void *data)
{
	struct kmbuild_bundle_t *bundle = (struct kmbuild_bundle_t *)data;
	struct kmhash_t *h = bundle->h;
	int ksize = bundle->ksize;
	int word_size = (ksize + 3) >> 2;
	// char *s = alloca(ksize + 2);
	// dump_kmer(kedge, ksize + 1, s);
	// __VERBOSE("kedge = %s\n", s);
	// if (s[0] != 'A')
	// 	__VERBOSE("ok kmer = %s\n", s);
	uint8_t *k1, *k2, *k1_rc, *k2_rc;
	kmint_t ik1, ik2;
	int c1, c2;
	// k1 = bundle->k1;
	// k2 = bundle->k2;
	// k1_rc = bundle->k1_rc;
	// k2_rc = bundle->k2_rc;
	k1 = alloca(word_size);
	k2 = alloca(word_size);
	k1_rc = alloca(word_size);
	k2_rc = alloca(word_size);
	kedge_get_left(k1, kedge, ksize, word_size);
	kedge_get_right(k2, kedge, ksize, word_size);
	km_get_rc(k1_rc, k1, ksize, word_size);
	km_get_rc(k2_rc, k2, ksize, word_size);

	c1 = kedge[0] & 0x3;
	c2 = (kedge[ksize >> 2] >> ((ksize & 0x3) << 1)) ^ 0x3;
	// assert(c2 == (nt4_table[(int)s[0]] ^ 3));
	// if (c2 != 3)
	// 	__VERBOSE("c1 = %d; c2 = %d\n", c1, c2);
	if (km_cmp(k1, k1_rc, word_size) <= 0) {
		kmhash_set_adj_multi(h, k1, c1, h->locks + thread_no);
		// ik1 = kmhash_put(h, k1);
		// kmhash_set_adj(h, ik1, c1);
	} else {
		kmhash_set_adj_multi(h, k1_rc, c1 + 4, h->locks + thread_no);
		// ik1 = kmhash_put(h, k1_rc);
		// kmhash_set_adj(h, ik1, c1 + 4);
	}

	if (km_cmp(k2, k2_rc, word_size) <= 0) {
		kmhash_set_adj_multi(h, k2, c2 + 4, h->locks + thread_no);
		// ik2 = kmhash_put(h, k2);
		// kmhash_set_adj(h, ik2, c2 + 4);
	} else {
		kmhash_set_adj_multi(h, k2_rc, c2, h->locks + thread_no);
		// ik2 = kmhash_put(h, k2_rc);
		// kmhash_set_adj(h, ik2, c2);
	}
}

void split_kmer_from_kedge(uint8_t *kedge, uint32_t count, void *data)
{
	struct kmbuild_bundle_t *bundle = (struct kmbuild_bundle_t *)data;
	struct kmhash_t *h = bundle->h;
	int ksize = bundle->ksize;
	int word_size = (ksize + 3) >> 2;
	// char *s = alloca(ksize + 2);
	// dump_kmer(kedge, ksize + 1, s);
	// if (s[0] != 'A')
	// 	__VERBOSE("ok kmer = %s\n", s);
	uint8_t *k1, *k2, *k1_rc, *k2_rc;
	kmint_t ik1, ik2;
	int c1, c2;
	k1 = bundle->k1;
	k2 = bundle->k2;
	k1_rc = bundle->k1_rc;
	k2_rc = bundle->k2_rc;
	kedge_get_left(k1, kedge, ksize, word_size);
	kedge_get_right(k2, kedge, ksize, word_size);
	km_get_rc(k1_rc, k1, ksize, word_size);
	km_get_rc(k2_rc, k2, ksize, word_size);

	c1 = kedge[0] & 0x3;
	c2 = (kedge[ksize >> 2] >> ((ksize & 0x3) << 1)) ^ 0x3;
	// assert(c2 == (nt4_table[(int)s[0]] ^ 3));
	// if (c2 != 3)
	// 	__VERBOSE("c1 = %d; c2 = %d\n", c1, c2);
	if (km_cmp(k1, k1_rc, word_size) <= 0) {
		ik1 = kmhash_put(h, k1);
		kmhash_set_adj(h, ik1, c1);
	} else {
		ik1 = kmhash_put(h, k1_rc);
		kmhash_set_adj(h, ik1, c1 + 4);
	}

	if (km_cmp(k2, k2_rc, word_size) <= 0) {
		ik2 = kmhash_put(h, k2);
		kmhash_set_adj(h, ik2, c2 + 4);
	} else {
		ik2 = kmhash_put(h, k2_rc);
		kmhash_set_adj(h, ik2, c2);
	}
}

void kmbuild_bundle_init(struct kmbuild_bundle_t *b, struct kmhash_t *h,
								int ksize)
{
	b->h = h;
	b->ksize = ksize;
	int l = (ksize + 3) >> 2;
	b->k1 = malloc(l);
	b->k2 = malloc(l);
	b->k1_rc = malloc(l);
	b->k2_rc = malloc(l);
}

void kmbuild_bundle_destroy(struct kmbuild_bundle_t *b)
{
	free(b->k1);
	free(b->k2);
	free(b->k1_rc);
	free(b->k2_rc);
	b->k1 = b->k2 = b->k1_rc = b->k2_rc = NULL;
}

void retrieve_kmer_from_kedge(struct opt_proc_t *opt, int ksize,
						struct kmhash_t *kmer_table)
{
	char *database_path = alloca(strlen(opt->out_dir) + 20);
	sprintf(database_path, "%s/KMC_%d_count", opt->out_dir, ksize + 1);
	char *kmc_pre = alloca(strlen(database_path) + 10);
	char *kmc_suf = alloca(strlen(database_path) + 10);
	strcpy(kmc_pre, database_path); strcat(kmc_pre, ".kmc_pre");
	strcpy(kmc_suf, database_path); strcat(kmc_suf, ".kmc_suf");
	struct kmc_info_t kmc_inf;
	KMC_read_prefix(kmc_pre, &kmc_inf);
	struct kmbuild_bundle_t bundle;
	kmbuild_bundle_init(&bundle, kmer_table, ksize);
	KMC_retrieve_kmer_multi(kmc_suf, opt->n_threads, &kmc_inf, (void *)(&bundle), split_kmer_from_kedge_multi);
	// KMC_retrive_kmer(kmc_suf, &kmc_inf, (void *)(&bundle), split_kmer_from_kedge);
	kmbuild_bundle_destroy(&bundle);
	destroy_kmc_info(&kmc_inf);
}

uint8_t fake_knum[8] = {178, 235, 171, 128, 178, 235, 170, 0};
uint8_t fake_krev[8] = {127, 21, 69, 92, 127, 5, 69, 28};

void find_bugs(struct kmhash_t *h, int ksize)
{
	kmint_t i, k;
	int word_size = (ksize + 3) >> 2;
	if (km_cmp(fake_knum, fake_krev, word_size) <= 0) {
		i = kmhash_get(h, fake_knum);
	} else {
		i = kmhash_get(h, fake_krev);
	}
	assert(i != KMHASH_END(h));
	__VERBOSE("[deb] adj = %hhu\n", KMHASH_ADJ(h, i));
}

void build_asm_graph_KMC(struct opt_proc_t *opt, int ksize, struct asm_graph_t *g)
{
	struct kmhash_t kmer_table;
	__VERBOSE("|-------- Retrieve kmer nodes\n");
	/* FIXME: parallel */
	kmhash_init(&kmer_table, opt->hash_size, (ksize + 3) >> 2, KM_AUX_ADJ, opt->n_threads);
	retrieve_kmer_from_kedge(opt, ksize, &kmer_table);
	// find_bugs(&kmer_table, ksize);
	__VERBOSE("Number of kmer = %lu\n", kmer_table.n_item);
	// test_kmhash(&kmer_table, ksize);
	__VERBOSE("|-------- Build kmer graph from kmer table\n");
	build_asm_graph_from_kmhash(opt->n_threads, ksize, &kmer_table, g);
	__VERBOSE("Number of nodes: %ld; Number of edges: %ld\n", g->n_v, g->n_e);

	kmhash_destroy(&kmer_table);
	__VERBOSE("|-------- Build edge index\n");
	kmhash_init(&kmer_table, opt->hash_size, (ksize + 4) >> 2, KM_AUX_IDX, opt->n_threads);
	/* FIXME: parallel */
	build_edge_kmer_index_multi(opt, &kmer_table, g);
	// build_edge_kmer_index(&kmer_table, g);
	__VERBOSE("|-------- Assign edge kmer count\n");
	/* FIXME: parallel */
	assign_edge_kmer_count(opt, &kmer_table, g);
	kmhash_destroy(&kmer_table);
}

void graph_build_KMC(struct opt_proc_t *opt, int ksize, struct asm_graph_t *g)
{
	set_time_now();
	__VERBOSE("|----Estimating kmer\n");
	char **tmp_files = alloca(opt->n_files * 2 * sizeof(char *));
	memcpy(tmp_files, opt->files_1, opt->n_files * sizeof(char *));
	memcpy(tmp_files + opt->n_files, opt->files_2, opt->n_files * sizeof(char *));
	KMC_build_kmer_database(ksize + 1, opt->out_dir, opt->n_threads, opt->mmem,
						opt->n_files * 2, tmp_files);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "Estimating kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("----Building assembly graph\n");
	build_asm_graph_KMC(opt, ksize, g);
	__VERBOSE_LOG("TIMER", "Building graph time: %.3f\n", sec_from_prev_time());
}

void build_local_graph(int ksize, int n_threads, int mmem, int n_files,
			char **files_1, char **files_2, char *work_dir,
						struct asm_graph_t *g)
{
	set_time_now();
	__VERBOSE("Building local graph\n");
	__VERBOSE("|---- Counting kmer\n");
	char **tmp_files = alloca(n_files * 2 * sizeof(char *));
	memcpy(tmp_files, files_1, n_files * sizeof(char *));
	memcpy(tmp_files + n_files, files_2, n_files * sizeof(char *));
	KMC_build_kmer_database(ksize + 1, work_dir, n_threads, mmem,
							n_files * 2, tmp_files);
	__VERBOSE("\n");
	__VERBOSE_LOG("TIMER", "KMC count kmer time: %.3f\n", sec_from_prev_time());
	set_time_now();

	__VERBOSE("|---- Retrieving kmer from KMC database\n");
	struct kmhash_t kmer_table;
	char *kmc_pre = alloca(strlen(work_dir) + 50);
	char *kmc_suf = alloca(strlen(work_dir) + 50);
	sprintf(kmc_pre, "%s/KMC_%d_count.kmc_pre", work_dir, ksize + 1);
	sprintf(kmc_suf, "%s/KMC_%d_count.kmc_suf", work_dir, ksize + 1);
	struct kmc_info_t kmc_inf;
	KMC_read_prefix(kmc_pre, &kmc_inf);
	struct kmbuild_bundle_t bundle;
	kmbuild_bundle_init(&bundle, &kmer_table, ksize);
	KMC_retrieve_kmer_multi(kmc_suf, n_threads, &kmc_inf, (void *)(&bundle),
						split_kmer_from_kedge_multi);
	kmbuild_bundle_destroy(&bundle);
	destroy_kmc_info(&kmc_inf);
	uint64_t size = kmer_table.size;
	__VERBOSE_LOG("BUILD", "Number of kmer: %lu\n", kmer_table.n_item);

	__VERBOSE("|-------- Build kmer graph from kmer table\n");
	build_asm_graph_from_kmhash(n_threads, ksize, &kmer_table, g);
	__VERBOSE("Number of nodes: %ld; Number of edges: %ld\n", g->n_v, g->n_e);

	kmhash_destroy(&kmer_table);
	__VERBOSE("|-------- Build edge index\n");
	kmhash_init(&kmer_table, opt->hash_size, (ksize + 4) >> 2, KM_AUX_IDX, opt->n_threads);
	/* FIXME: parallel */
	build_edge_kmer_index_multi(opt, &kmer_table, g);
	// build_edge_kmer_index(&kmer_table, g);
	__VERBOSE("|-------- Assign edge kmer count\n");
	/* FIXME: parallel */
	assign_edge_kmer_count(opt, &kmer_table, g);
	kmhash_destroy(&kmer_table);

	__VERBOSE("|---- Building graph connection\n");
	build_asm_graph_from_kmhash(n_threads, ksize, &kmer_table, g);
	__VERBOSE_LOG("BUILD", "Number of nodes: %ld; Number of edges: %ld\n", g->n_v, g->n_e);
	kmhash_destroy(&kmer_table);

	__VERBOSE("|---- Assigning edge count\n");
	kmhash_init(&kmer_table, size, (ksize + 4) >> 2, KM_AUX_IDX, n_threads);
	build_edge_kmer_index_multi(opt, &kmer_table, g);

	assign_edge_kmer_count(opt, &kmer_table, g);
	kmhash_destroy(&kmer_table);
}

struct kmedge_bundle_t {
	struct kmhash_t *h;
	struct asm_graph_t *g;
};

void assign_count_kedge(uint8_t *kmer, uint32_t count, void *data)
{
	struct kmedge_bundle_t *bundle = (struct kmedge_bundle_t *)data;
	struct kmhash_t *h = bundle->h;
	struct asm_graph_t *g = bundle->g;
	int ksize = g->ksize + 1;
	int word_size = (ksize + 3) >> 2;
	kmint_t k = kmhash_get(h, kmer);
	if (k == KMHASH_END(h)) {
		return;
	}
	gint_t e = KMHASH_IDX(h, k);
	g->edges[e].count += count;
	g->edges[g->edges[e].rc_id].count += count;
}

void assign_count_kedge_multi(int thread_no, uint8_t *kmer, uint32_t count, void *data)
{
	struct kmedge_bundle_t *bundle = (struct kmedge_bundle_t *)data;
	struct kmhash_t *h = bundle->h;
	struct asm_graph_t *g = bundle->g;
	int ksize = g->ksize + 1;
	int word_size = (ksize + 3) >> 2;
	kmint_t k = kmhash_get(h, kmer);
	if (k == KMHASH_END(h)) {
		return;
	}
	gint_t e = KMHASH_IDX(h, k);
	atomic_add_and_fetch64(&g->edges[e].count, count);
	atomic_add_and_fetch64(&g->edges[g->edges[e].rc_id].count, count);
	// g->edges[e].count += count;
	// g->edges[g->edges[e].rc_id].count += count;
}

void assign_edge_kmer_count(struct opt_proc_t *opt, struct kmhash_t *h,
							struct asm_graph_t *g)
{
	char *database_path = alloca(strlen(opt->out_dir) + 20);
	sprintf(database_path, "%s/KMC_%d_count", opt->out_dir, g->ksize + 1);
	char *kmc_pre = alloca(strlen(database_path) + 10);
	char *kmc_suf = alloca(strlen(database_path) + 10);
	strcpy(kmc_pre, database_path); strcat(kmc_pre, ".kmc_pre");
	strcpy(kmc_suf, database_path); strcat(kmc_suf, ".kmc_suf");
	struct kmc_info_t kmc_inf;
	KMC_read_prefix(kmc_pre, &kmc_inf);
	struct kmedge_bundle_t bundle;
	bundle.h = h;
	bundle.g = g;
	KMC_retrieve_kmer_multi(kmc_suf, opt->n_threads, &kmc_inf, (void *)(&bundle), assign_count_kedge_multi);
	// KMC_retrive_kmer(kmc_suf, &kmc_inf, (void *)(&bundle), assign_count_kedge);
	destroy_kmc_info(&kmc_inf);
}

static inline void deb_dump_bin_seq(const char *label, uint32_t *bin, uint32_t len)
{
	char *seq;
	uint32_t i;
	seq = malloc(len + 1);
	for (i = 0; i < len; ++i) {
		seq[i] = nt4_char[(bin[i >> 4] >> ((i & 15) << 1)) & (uint32_t)0x3];
	}
	seq[len] = '\0';
	fprintf(stderr, "%s%s\n", label, seq);
	free(seq);
}

struct iedge_bundle_t {
	struct asm_graph_t *g;
	struct kmhash_t *h;
	gint_t lo_e;
	gint_t hi_e;
	pthread_mutex_t *lock;
};

void *build_edge_index_worker(void *data)
{
	struct iedge_bundle_t *bundle = (struct iedge_bundle_t *)data;
	struct kmhash_t *h = bundle->h;
	struct asm_graph_t *g = bundle->g;
	gint_t lo_e, hi_e;
	lo_e = bundle->lo_e;
	hi_e = bundle->hi_e;
	// __VERBOSE("hi_e = %ld; lo_e = %ld\n", hi_e, lo_e);
	pthread_mutex_t *lock = bundle->lock;
	int ksize, word_size;
	ksize = g->ksize + 1;
	word_size = (ksize + 3) >> 2;
	uint8_t *knum, *krev;
	knum = alloca(word_size);
	krev = alloca(word_size);
	gint_t e, e_rc;
	for (e = lo_e; e < hi_e; ++e) {
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		uint32_t i;
		memset(knum, 0, word_size);
		memset(krev, 0, word_size);
		for (i = 0; i < g->edges[e].seq_len; ++i) {
			uint32_t c = __binseq_get(g->edges[e].seq, i);
			km_shift_append(knum, ksize, word_size, c);
			km_shift_append_rv(krev, ksize, word_size, c ^ 3);
			if (i + 1 < (uint32_t)ksize)
				continue;
			kmint_t k;
			if (km_cmp(knum, krev, word_size) <= 0) {
				kmhash_set_idx_multi(h, knum, e, lock);
				// k = kmhash_put(h, knum);
				// KMHASH_IDX(h, k) = e;
			} else {
				kmhash_set_idx_multi(h, krev, e, lock);
				// k = kmhash_put(h, krev);
				// KMHASH_IDX(h, k) = e;
			}
		}
	}
	return NULL;
}

void build_edge_kmer_index_multi(struct opt_proc_t *opt, struct kmhash_t *h, struct asm_graph_t *g)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	struct iedge_bundle_t *bundle;
	bundle = calloc(opt->n_threads, sizeof(struct iedge_bundle_t));

	gint_t e, e_rc, cur_e, cap, sum = 0;
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		sum += g->edges[e].seq_len - g->ksize;
	}
	cap = sum / opt->n_threads + 1;
	int k = 0;
	sum = 0;
	bundle[0].lo_e = 0;
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		sum += g->edges[e].seq_len - g->ksize;
		if (sum >= cap * (k + 1)) {
			bundle[k].hi_e = e + 1;
			bundle[k + 1].lo_e = e + 1;
			++k;
		}
	}
	bundle[opt->n_threads - 1].hi_e = g->n_e;

	for (k = 0; k < opt->n_threads; ++k) {
		bundle[k].h = h;
		bundle[k].g = g;
		bundle[k].lock = h->locks + k;
	}

	pthread_t *threads = calloc(opt->n_threads, sizeof(pthread_t));
	for (k = 0; k < opt->n_threads; ++k)
		pthread_create(threads + k, &attr, build_edge_index_worker, bundle + k);
	for (k = 0; k < opt->n_threads; ++k)
		pthread_join(threads[k], NULL);
	__VERBOSE("Number of (k+1)-mer: %lu\n", h->n_item);
	free(threads);
	free(bundle);
}

void build_edge_kmer_index(struct kmhash_t *h, struct asm_graph_t *g)
{
	int ksize, word_size;
	ksize = g->ksize + 1;
	word_size = (ksize + 3) >> 2;
	uint8_t *knum, *krev;
	knum = alloca(word_size);
	krev = alloca(word_size);
	gint_t e, e_rc;
	gint_t n_pos = 0;
	for (e = 0; e < g->n_e; ++e) {
		e_rc = g->edges[e].rc_id;
		if (e > e_rc)
			continue;
		uint32_t i;
		memset(knum, 0, word_size);
		memset(krev, 0, word_size);
		for (i = 0; i < g->edges[e].seq_len; ++i) {
			uint32_t c = __binseq_get(g->edges[e].seq, i);
			km_shift_append(knum, ksize, word_size, c);
			km_shift_append_rv(krev, ksize, word_size, c ^ 3);
			if (i + 1 < (uint32_t)ksize)
				continue;
			++n_pos;
			kmint_t k;
			if (km_cmp(knum, krev, word_size) <= 0) {
				k = kmhash_put(h, knum);
				KMHASH_IDX(h, k) = e;
			} else {
				k = kmhash_put(h, krev);
				KMHASH_IDX(h, k) = e;
			}
		}
	}
}

struct kmgraph_bundle_t {
	struct asm_graph_t *g;
	struct kmhash_t *h;
	int ksize;
	gint_t lo_e;
	kmint_t l;
	kmint_t r;
};

static void asm_init_edge(struct asm_edge_t *e, const uint8_t *kmer, int ksize)
{
	int m = (ksize + 15) / 16, i, k;
	uint32_t c;
	e->seq = calloc(m, sizeof(uint32_t));
	for (i = 0; i < ksize; ++i) {
		k = ksize - i - 1;
		c = __kmerseq_get(kmer, k);
		__binseq_set(e->seq, i, c);
	}
	e->seq_len = ksize;
}

static void asm_append_edge_char(struct asm_edge_t *e, uint32_t c)
{
	if (!(e->seq_len & 15)) {
		e->seq = realloc(e->seq, ((e->seq_len >> 4) + 1) * sizeof(uint32_t));
		e->seq[e->seq_len >> 4] = 0;
	}
	e->seq[e->seq_len >> 4] |= c << ((e->seq_len & 15) << 1);
	++e->seq_len;
}

static void *build_graph_worker(void *data)
{
	struct kmgraph_bundle_t *bundle = (struct kmgraph_bundle_t *)data;
	struct asm_graph_t *g = bundle->g;
	struct kmhash_t *h = bundle->h;
	int ksize = bundle->ksize;
	gint_t lo_e = bundle->lo_e;
	kmint_t it_l = bundle->l;
	kmint_t it_r = bundle->r;

	struct asm_node_t *nodes = g->nodes;
	struct asm_edge_t *edges = g->edges;
	int word_size = (ksize + 3) >> 2;
	uint8_t *knum, *krev, *cur_knum, *cur_krev;
	krev = alloca(word_size);
	cur_knum = alloca(word_size);
	cur_krev = alloca(word_size);
	// __VERBOSE("it_l = %ld; it_r = %ld\n", it_l, it_r);

	kmint_t i, k;
	for (i = it_l; i < it_r; ++i) {
		if (KM_IS_EMPTY(h, i))
			continue;
		uint8_t adj_fw, adj_rv, cur_adj_fw, cur_adj_rv;
		int deg_fw, deg_rv, cur_deg_fw, cur_deg_rv;
		gint_t knum_idx, krev_idx;
		int c, cur_c;

		adj_fw = KMHASH_ADJ(h, i) & 0xf;
		adj_rv = KMHASH_ADJ(h, i) >> 4;
		deg_fw = __bin_degree4(adj_fw);
		deg_rv = __bin_degree4(adj_rv);
		if (deg_fw == 1 && deg_rv == 1) /* non-branching nodes */
			continue;
		knum = KMHASH_KEY(h, i);
		km_get_rc(krev, knum, ksize, word_size);
		knum_idx = KMHASH_IDX(h, i) * 2;
		krev_idx = KMHASH_IDX(h, i) * 2 + 1;

		deg_fw = deg_rv = 0;
		for (c = 0; c < 4; ++c) {
			if (!((adj_fw >> c) & 1))
				continue;
			struct asm_edge_t *e = edges + lo_e;
			asm_init_edge(e, knum, ksize);
			memcpy(cur_knum, knum, word_size);
			memcpy(cur_krev, krev, word_size);
			cur_c = c;
			do {
				asm_append_edge_char(e, cur_c);
				km_shift_append(cur_knum, ksize, word_size, cur_c);
				km_shift_append_rv(cur_krev, ksize, word_size, cur_c ^ 3);
				if (km_cmp(cur_knum, cur_krev, word_size) <= 0) {
					k = kmhash_get(h, cur_knum);
					assert(k != KMHASH_END(h));
					cur_adj_fw = KMHASH_ADJ(h, k) & 0xf;
					cur_adj_rv = KMHASH_ADJ(h, k) >> 4;
					cur_deg_fw = __bin_degree4(cur_adj_fw);
					cur_deg_rv = __bin_degree4(cur_adj_rv);
					if (cur_deg_fw == 1 && cur_deg_rv == 1)
						cur_c = __bin_only4(cur_adj_fw);
				} else {
					k = kmhash_get(h, cur_krev);
					assert(k != KMHASH_END(h));
					cur_adj_fw = KMHASH_ADJ(h, k) & 0xf;
					cur_adj_rv = KMHASH_ADJ(h, k) >> 4;
					cur_deg_fw = __bin_degree4(cur_adj_fw);
					cur_deg_rv = __bin_degree4(cur_adj_rv);
					if (cur_deg_fw == 1 && cur_deg_rv == 1)
						cur_c = __bin_only4(cur_adj_rv);
				}
			} while (cur_deg_fw == 1 && cur_deg_rv == 1);
			e->source = knum_idx;
			if (km_cmp(cur_knum, cur_krev, word_size) <= 0)
				e->target = KMHASH_IDX(h, k) * 2;
			else
				e->target = KMHASH_IDX(h, k) * 2 + 1;
			nodes[knum_idx].adj[deg_fw++] = lo_e++;
		}

		for (c = 0; c < 4; ++c) {
			if (!((adj_rv >> c) & 1))
				continue;
			struct asm_edge_t *e = edges + lo_e;
			asm_init_edge(e, krev, ksize);
			memcpy(cur_knum, krev, word_size);
			memcpy(cur_krev, knum, word_size);
			cur_c = c;
			do {
				asm_append_edge_char(e, cur_c);
				km_shift_append(cur_knum, ksize, word_size, cur_c);
				km_shift_append_rv(cur_krev, ksize, word_size, cur_c ^ 3);
				if (km_cmp(cur_knum, cur_krev, word_size) <= 0) {
					k = kmhash_get(h, cur_knum);
					assert(k != KMHASH_END(h));
					cur_adj_fw = KMHASH_ADJ(h, k) & 0xf;
					cur_adj_rv = KMHASH_ADJ(h, k) >> 4;
					cur_deg_fw = __bin_degree4(cur_adj_fw);
					cur_deg_rv = __bin_degree4(cur_adj_rv);
					if (cur_deg_fw == 1 && cur_deg_rv == 1)
						cur_c = __bin_only4(cur_adj_fw);
				} else {
					k = kmhash_get(h, cur_krev);
					assert(k != KMHASH_END(h));
					cur_adj_fw = KMHASH_ADJ(h, k) & 0xf;
					cur_adj_rv = KMHASH_ADJ(h, k) >> 4;
					cur_deg_fw = __bin_degree4(cur_adj_fw);
					cur_deg_rv = __bin_degree4(cur_adj_rv);
					if (cur_deg_fw == 1 && cur_deg_rv == 1)
						cur_c = __bin_only4(cur_adj_rv);
				}
			} while (cur_deg_fw == 1 && cur_deg_rv == 1);
			e->source = krev_idx;
			if (km_cmp(cur_knum, cur_krev, word_size) <= 0)
				e->target = KMHASH_IDX(h, k) * 2;
			else
				e->target = KMHASH_IDX(h, k) * 2 + 1;
			nodes[krev_idx].adj[deg_rv++] = lo_e++;
		}
	}
	return NULL;
}

void build_asm_graph_from_kmhash(int n_threads, int ksize,
				struct kmhash_t *h, struct asm_graph_t *g)
{
	kmhash_alloc_aux(h, KM_AUX_IDX);
	gint_t n_v, n_e;
	kmint_t i;
	n_v = n_e = 0;
	uint8_t adj_fw, adj_rv;
	int deg_fw, deg_rv;
	for (i = 0; i < h->size; ++i) {
		if (KM_IS_EMPTY(h, i))
			continue;
		adj_fw = KMHASH_ADJ(h, i) & 0xf;
		adj_rv = KMHASH_ADJ(h, i) >> 4;
		deg_fw = __bin_degree4(adj_fw);
		deg_rv = __bin_degree4(adj_rv);
		// __VERBOSE("deg_fw = %d; deg_rv = %d\n", deg_fw, deg_rv);
		if (deg_fw == 1 && deg_rv == 1)
			continue;
		n_e += (deg_fw + deg_rv);
		KMHASH_IDX(h, i) = n_v++;
	}

	struct asm_node_t *nodes = calloc(n_v * 2, sizeof(struct asm_node_t));
	struct asm_edge_t *edges = calloc(n_e, sizeof(struct asm_edge_t));
	g->nodes = nodes;
	g->edges = edges;
	g->ksize = ksize;
	g->aux_flag = 0;
	g->n_e = n_e;
	g->n_v = n_v * 2;
	g->bin_size = 0;

	struct kmgraph_bundle_t *bundle;
	bundle = calloc(n_threads, sizeof(struct kmgraph_bundle_t));

	kmint_t cap = h->size / n_threads + 1;
	n_e = 0;
	int k;
	for (k = 0; k < n_threads; ++k) {
		bundle[k].h = h;
		bundle[k].g = g;
		bundle[k].ksize = ksize;
		bundle[k].l = k * cap;
		bundle[k].r = __min((k + 1) * cap, h->size);
		bundle[k].lo_e = n_e;
		for (i = bundle[k].l; i < bundle[k].r; ++i) {
			if (KM_IS_EMPTY(h, i))
				continue;
			adj_fw = KMHASH_ADJ(h, i) & 0xf;
			adj_rv = KMHASH_ADJ(h, i) >> 4;
			deg_fw = __bin_degree4(adj_fw);
			deg_rv = __bin_degree4(adj_rv);
			if (deg_fw == 1 && deg_rv == 1)
				continue;
			n_e += (deg_fw + deg_rv);
			gint_t knum_idx, krev_idx;
			knum_idx = KMHASH_IDX(h, i) * 2;
			krev_idx = KMHASH_IDX(h, i) * 2 + 1;
			nodes[knum_idx].rc_id = krev_idx;
			nodes[krev_idx].rc_id = knum_idx;
			nodes[knum_idx].adj = malloc(deg_fw * sizeof(gint_t));
			nodes[krev_idx].adj = malloc(deg_rv * sizeof(gint_t));
			nodes[knum_idx].deg = deg_fw;
			nodes[krev_idx].deg = deg_rv;
		}
	}
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_t *threads = calloc(n_threads, sizeof(pthread_t));
	for (k = 0; k < n_threads; ++k)
		pthread_create(threads + k, &attr, build_graph_worker, bundle + k);
	for (k =0; k < n_threads; ++k)
		pthread_join(threads[k], NULL);
	free(threads);
	free(bundle);

	/* link reverse complement edges */
	gint_t e, e_rc, v, v_rc;
	for (e = 0; e < n_e; ++e) {
		edges[e].rc_id = -1;
		v = edges[e].target;
		v_rc = nodes[v].rc_id;
		for (k = 0; k < nodes[v_rc].deg; ++k) {
			e_rc = nodes[v_rc].adj[k];
			if (edges[e_rc].target == nodes[edges[e].source].rc_id
				&& is_seq_rc(edges[e].seq, edges[e].seq_len,
					edges[e_rc].seq, edges[e_rc].seq_len)) {
				edges[e].rc_id = e_rc;
				edges[e_rc].rc_id = e;
				break;
			}
		}
		assert(edges[e].rc_id != -1);
	}

	gint_t u, u_rc;
	for (u = 0; u < g->n_v; ++u) {
		u_rc = g->nodes[u].rc_id;
		if (g->nodes[u].deg + g->nodes[u_rc].deg == 0)
			assert(0);
	}
}

