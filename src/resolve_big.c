//
// Created by che on 31/10/2019.
//

#include "assembly_graph.h"
#include "khash.h"
#include "resolve_big.h"
#include<stdio.h>
#include <assert.h>
#include "kmhash.h"
#include "radix_sort.h"
#include "sort_read.h"
#include "buffer_file_wrapper.h"
#include "io_utils.h"
#include "utils.h"
#include "KMC_reader.h"
#include "kmer_build.h"
#include "log.h"
#include "barcode_resolve2.h"
#include "atomic.h"

#if defined(_MSC_VER)

#define FORCE_INLINE	__forceinline

#include <stdlib.h>

#define ROTL32(x,y)	_rotl(x,y)
#define ROTL64(x,y)	_rotl64(x,y)

#define BIG_CONSTANT(x) (x)

// Other compilers

#else	// defined(_MSC_VER)

#define	FORCE_INLINE inline __attribute__((always_inline))

inline uint32_t rotl32(uint32_t x, int8_t r)
{
	return (x << r) | (x >> (32 - r));
}

static inline uint64_t rotl64(uint64_t x, int8_t r)
{
	return (x << r) | (x >> (64 - r));
}

#define	ROTL32(x,y)	rotl32(x,y)
#define ROTL64(x,y)	rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)

#endif // !defined(_MSC_VER)

static inline uint64_t getblock64(const uint64_t * p, int i)
{
	return p[i];
}

static inline uint64_t fmix64(uint64_t k)
{
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xff51afd7ed558ccd);
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
	k ^= k >> 33;

	return k;
}

inline int get_nu(const uint32_t *seq, int pos)
{
	return (seq[pos >> 4] >> ((pos & 15) << 1)) & 3;
}


void test(struct asm_edge_t *e, khash_t(big_kmer_count) *table)
{
	int new_len = MIN(81, e->seq_len);
	uint8_t *new_seq = calloc((new_len + 3) >> 2, 1);
	for (int i = 0; i < new_len; i++) {
		new_seq[i >> 2] |= (e->seq[i >> 4] >> ((i & 15) << 1) & 3) << ((i & 3) << 1);
	}
	gint_t hash_key = MurmurHash3_x64_64(new_seq, (new_len+3) >> 2);
	khint_t k = kh_get(big_kmer_count, table, hash_key);
	int value;
	if (k == kh_end(table))
		value = 0;
	else
		value = kh_value(table, k);
	log_warn("this must not 0: %d", value);
}

void print_compress_seq(uint8_t *new_seq, int len)
{
	int to_char[4] = {65, 67,71,84};
	for (int i= 0; i < len; i++) {
		printf("%c", to_char[new_seq[i>>2] >> ((i&3)<<1) &3]);
	}
	printf("\n");
}

/**
 *						a	    e         b
 *					 o----1-ksize-o--------o-ksize-1-------o
 * @param a
 * @param b
 * @param e
 * @return
 */
int get_new_seq_count(struct asm_edge_t *a, struct asm_edge_t *b, struct asm_edge_t *e, int ksize,
					  khash_t(big_kmer_count) *table)
{
	test(e, table); //todo 0 remove this
	int new_len = e->seq_len + 2;
	uint8_t *new_seq = calloc((new_len + 3) >> 2, 1);
	int x_a = get_nu(a->seq, a->seq_len - ksize - 1);
	new_seq[0] |= x_a;
	for (int i = 0; i < e->seq_len; i++) {
		new_seq[(i + 1) >> 2] |= (e->seq[i >> 4] >> ((i & 15) << 1) & 3) << (((i + 1) & 3) << 1);
	}
	int x_b = get_nu(b->seq, ksize);
	new_seq[((new_len + 3) >> 2) -1] |= x_b << ((new_len & 3) << 1);
//	print_compress_seq(new_seq, new_len); //todo 0 remove this
	int value = 0;
	gint_t hash_key = MurmurHash3_x64_64(new_seq, (new_len+3)>>2);
	khint_t k = kh_get(big_kmer_count, table, hash_key);
	if (k == kh_end(table))
		value = 0;
	else
		value = kh_value(table, k);
	return value;
}

int append_barcode_contig(struct barcode_hash_t *bc, khash_t(union_barcode) *uni_bar)
{
	int count = 0;
	for (uint32_t i = 0; i < bc->size; i++) {
		if (bc->keys[i] == (uint64_t) -1) {
			continue;
		}
		count++;
		int ret = 0;
		kh_put(union_barcode, uni_bar, bc->keys[i], &ret);
		assert(ret != -1);
	}
	assert(count == bc->n_item);
	return count;
}

int is_case_2_1_2(struct asm_graph_t *g, int i_e)
{
	int source = g->edges[i_e].source;
	if (source == -1)
		return 0;
	if (g->edges[i_e].seq_len > 150) {
		return 0;
	}
	int target = g->edges[i_e].target;
	int source_rc = g->nodes[source].rc_id;
	if (g->nodes[target].deg != 2 || g->nodes[source_rc].deg != 2) {
		return 0;
	}
	int i_a0 = g->nodes[g->nodes[source].rc_id].adj[0];
	int i_a1 = g->nodes[g->nodes[source].rc_id].adj[1];
	int i_o0 = g->nodes[target].adj[0];
	int i_o1 = g->nodes[target].adj[1];
	struct asm_edge_t *a0 = &g->edges[i_a0];
	struct asm_edge_t *a1 = &g->edges[i_a1];
	struct asm_edge_t *o0 = &g->edges[i_o0];
	struct asm_edge_t *o1 = &g->edges[i_o1];
	struct asm_edge_t *e = &g->edges[i_e];

	log_warn("len of: a0 %d a1 %d e %d o0 %d o1%d\n", a0->seq_len, a1->seq_len, e->seq_len, o0->seq_len, o1->seq_len);
	if (a0->rc_id == i_o0 || a0->rc_id == i_o1 || a1->rc_id == i_o0 || a1->rc_id == i_o1)
		return 0;
	return 1;
}

int dfs_partition(struct asm_graph_t *g, int *partition, int i_e, int index_par, khash_t(union_barcode) *uni_bar,
				  int *total_barcodes, int *total_length_component, int *count_resolvable)
{
//	log_warn("dfs from %d", i_e);
	if (partition[i_e] != -1)
		return 0;
	if (g->edges[i_e].seq_len >= CONTIG_PARTITION_LEN) {
		log_warn("bc0&1&2 %d %d %d", g->edges[i_e].barcodes[0].n_item , g->edges[i_e].barcodes[1].n_item, g->edges[i_e].barcodes[2].n_item);

//		assert(g->edges[i_e].barcodes[0].n_item <= g->edges[i_e].barcodes[2].n_item);
//		assert(g->edges[i_e].barcodes[0].n_item <= g->edges[i_e].barcodes[1].n_item);
		*total_barcodes += append_barcode_contig(&g->edges[i_e].barcodes[2], uni_bar);//todo 1 all barcode
		return 0;
	}
	*count_resolvable += is_case_2_1_2(g, i_e);
	*total_length_component += g->edges[i_e].seq_len;
//	if ((*uni_bar)->size > 1500000 ){
//		return 0;
//	}
	int res = 1;
	int x_rc = g->edges[i_e].rc_id;
	partition[i_e] = index_par;
	partition[x_rc] = index_par;
	//todo 1 build all barcode

	int source = g->edges[i_e].source;
	int source_rc = g->nodes[source].rc_id;
	int target = g->edges[i_e].target;

	for (int i = 0; i < g->nodes[source_rc].deg; i++) {
		int i_e_rc = g->nodes[source_rc].adj[i];
		res += dfs_partition(g, partition, i_e_rc, index_par, uni_bar, total_barcodes, total_length_component, count_resolvable);
	}
	for (int i = 0; i < g->nodes[target].deg; i++) {
		int i_e = g->nodes[target].adj[i];
		res += dfs_partition(g, partition, i_e, index_par, uni_bar, total_barcodes, total_length_component, count_resolvable);
	}
	return res;
}

void get_pos(khash_t(bcpos) *dict, khash_t(union_barcode) *bc, struct read_index_t **read_pos, int *n_barcodes)
{
	log_warn("start get pos");
	int count = 0;
	struct read_index_t *pos = NULL;
	for (khint_t it = kh_begin(bc); it != kh_end(bc); it++){
		if (!(kh_exist(bc, it))) {
			continue;
		}

		gint_t key = kh_key(bc, it);
		khiter_t k = kh_get(bcpos, dict, key);
		pos = realloc(pos, (count +1) * sizeof(struct read_index_t));
		pos[count] = kh_value(dict, k);
		count++;
		if (count > MAX_BARCODE_REGION) {
			break;
		}
	}
	log_warn("start sort read index");
	rs_sort(read_index, pos, pos + count);
	*n_barcodes = count;
	*read_pos = pos;
}

void km_count_bundle_init(struct km_count_bundle_t *km_count_bundle, khash_t(big_kmer_count) *kmer_count_table,
						  int ksize)
{
	km_count_bundle->kmer_count_table = kmer_count_table;
	km_count_bundle->ksize = ksize;
	pthread_mutex_init( &km_count_bundle->lock, NULL);
}

static inline uint64_t MurmurHash3_x64_64(const uint8_t *data, const int len)
{
	int n_blocks = len >> 4;
	uint64_t h1 = BIG_CONSTANT(0x13097);
	uint64_t h2 = BIG_CONSTANT(0x13097);

	const uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
	const uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

	const uint64_t *blocks = (const uint64_t *)(data);

	int i;
	for (i = 0; i < n_blocks; ++i) {
		uint64_t k1 = getblock64(blocks, i << 1);
		uint64_t k2 = getblock64(blocks, (i << 1) + 1);
		k1 *= c1; k1 = ROTL64(k1, 31); k1 *= c2; h1 ^= k1;
		h1 = ROTL64(h1, 27); h1 += h2; h1 = h1 * 5 + 0x52dce729;
		k2 *= c2; k2 = ROTL64(k2, 33); k2 *= c1; h2 ^= k2;
		h2 = ROTL64(h2, 31); h2 += h1; h2 = h2 * 5 + 0x38495ab5;
	}

	const uint8_t *tail = (const uint8_t *)(data + (n_blocks << 4));
	uint64_t k1 = 0;
	uint64_t k2 = 0;

	switch (len & 15) {
		case 15: k2 ^= ((uint64_t)tail[14]) << 48;
		case 14: k2 ^= ((uint64_t)tail[13]) << 40;
		case 13: k2 ^= ((uint64_t)tail[12]) << 32;
		case 12: k2 ^= ((uint64_t)tail[11]) << 24;
		case 11: k2 ^= ((uint64_t)tail[10]) << 16;
		case 10: k2 ^= ((uint64_t)tail[ 9]) << 8;
		case  9: k2 ^= ((uint64_t)tail[ 8]) << 0;
			k2 *= c2; k2  = ROTL64(k2,33); k2 *= c1; h2 ^= k2;

		case  8: k1 ^= ((uint64_t)tail[ 7]) << 56;
		case  7: k1 ^= ((uint64_t)tail[ 6]) << 48;
		case  6: k1 ^= ((uint64_t)tail[ 5]) << 40;
		case  5: k1 ^= ((uint64_t)tail[ 4]) << 32;
		case  4: k1 ^= ((uint64_t)tail[ 3]) << 24;
		case  3: k1 ^= ((uint64_t)tail[ 2]) << 16;
		case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;
		case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;
			k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
	}

	h1 ^= len; h2 ^= len;

	h1 += h2;
	h2 += h1;

	h1 = fmix64(h1);
	h2 = fmix64(h2);

	h1 += h2;
	h2 += h1;

	return h2;
}

void count_to_big_table(int thread_no, uint8_t *kmer, uint32_t count, void *data)
{
	struct km_count_bundle_t *bundle = (struct km_count_bundle_t *) data;
	khash_t(big_kmer_count) *table = bundle->kmer_count_table;

	int ksize = bundle->ksize;
	uint64_t k = MurmurHash3_x64_64(kmer, (ksize + 3) / 4);
	khint_t h = kh_get(big_kmer_count, table, k);
	if (h == kh_end(table)) {
		int tmp;
		pthread_mutex_lock(&bundle->lock);
		h = kh_put(big_kmer_count, table, k, &tmp);
		kh_value(table, h) = 0;
		pthread_mutex_unlock(&bundle->lock);
	}
	int t = kh_value(table, h);
	atomic_add_and_fetch32( &kh_value(table, h), count);
//	log_warn("value before add %d %d after %d",t, count, kh_value(table, h));
}

void count_kmer_from_read_file(struct read_path_t *read_file, khash_t(big_kmer_count) *kmer_count_table, int ksize,
							   int n_threads, char *work_dir, int mmem)
{
	char **tmp_files = alloca(2 * sizeof(char *));
	tmp_files[0] = read_file->R1_path;
	tmp_files[1] = read_file->R2_path;
	KMC_build_kmer_database(ksize + 1, work_dir, n_threads, mmem, 2, tmp_files);
	log_debug("|---- Retrieving kmer from KMC database");
	struct kmc_info_t kmc_inf;
	char *kmc_pre = alloca(strlen(work_dir) + 50);
	char *kmc_suf = alloca(strlen(work_dir) + 50);
	sprintf(kmc_pre, "%s/KMC_%d_count.kmc_pre", work_dir, ksize + 1);
	sprintf(kmc_suf, "%s/KMC_%d_count.kmc_suf", work_dir, ksize + 1);
	KMC_read_prefix(kmc_pre, &kmc_inf);

	struct km_count_bundle_t km_count_bundle;
	km_count_bundle_init(&km_count_bundle, kmer_count_table, ksize);
	KMC_retrieve_kmer_multi(kmc_suf, n_threads, &kmc_inf,
							(void *) (&km_count_bundle), count_to_big_table);
	km_count_bundle_destroy(&km_count_bundle);
}

void filter_read_build_kmer(struct read_path_t *ori_read, khash_t(bcpos) *dict,
							khash_t(union_barcode) *bc, khash_t(big_kmer_count) *kmer_count_table,
							int min_ksize, int max_ksize, int n_threads, int mmem, int partition_index)
{
	struct read_index_t *pos = NULL;
	int n_barcodes =0;
	get_pos(dict, bc, &pos, &n_barcodes);

	struct read_path_t *ans = calloc(1, sizeof(struct read_path_t));
	ans->R1_path = calloc(100, 1);
	ans->R1_path = str_concate(ori_read->R1_path, "aaa.fq");
	ans->R2_path = calloc(100, 1);
	ans->R2_path = str_concate(ori_read->R2_path, "aaa.fq");
	log_warn("%s %s", ans->R1_path, ans->R2_path);

	struct buffered_file_t fo1, fo2;
	bf_open(&fo1, ans->R1_path, "wb", SIZE_16MB);
	bf_open(&fo2, ans->R2_path, "wb", SIZE_16MB);
	FILE *fi1 = xfopen(ori_read->R1_path, "rb");
	FILE *fi2 = xfopen(ori_read->R2_path, "rb");
	int64_t m_buf, len;
	m_buf = 0x100;
	char *buf = malloc(m_buf);
	log_warn("Start read data from %d barcodes", n_barcodes);
	n_barcodes = MIN(n_barcodes, 5000);
	for (int i = 0; i < n_barcodes; ++i) {
//		log_warn("r1len %d r2len %d r1off %d r2off %d", pos[i].r1_len, pos[i].r2_len, pos[i].r1_offset, pos[i].r2_offset);
		len = __max(pos[i].r1_len, pos[i].r2_len);
		if (len > m_buf) {
			m_buf = len;
			buf = realloc(buf, m_buf);
		}
		fseek(fi1, pos[i].r1_offset, SEEK_SET);
		xfread(buf, 1, pos[i].r1_len, fi1);
		bf_write(&fo1, buf, pos[i].r1_len);
		fseek(fi2, pos[i].r2_offset, SEEK_SET);
		xfread(buf, 1, pos[i].r2_len, fi2);
		bf_write(&fo2, buf, pos[i].r2_len);
	}
	fclose(fi1);
	fclose(fi2);
	bf_close(&fo1);
	bf_close(&fo2);
	free(buf);
	free(pos);

	char *work_dir = calloc(100, 1);
	int last = 0;
	for (int i = 0; i < strlen(ans->R1_path); i++) {
		if (ans->R1_path[i] == '/') {
			last = i;
		}
	}
	for (int i = 0; i < last; i++) {
		work_dir[i] = ans->R1_path[i];
	}
	char *tmp = calloc(100, 1);
	sprintf(tmp, "/%d", partition_index);
	work_dir = str_concate(work_dir, tmp);
	log_warn("work dir %s ", work_dir);
	mkdir(work_dir, 0755);

	for (int ksize = min_ksize; ksize < max_ksize; ksize++) {
		char *this_work_dir = calloc(100, 1);
		sprintf(this_work_dir, "%s/%d", work_dir, ksize);
		log_warn("this work dir %s ", this_work_dir);
		mkdir(this_work_dir, 0755);
		count_kmer_from_read_file(ans, kmer_count_table, ksize, n_threads, this_work_dir, mmem);
	}
	log_warn("kmer table size %d", kmer_count_table->size);
}

void dfs_resolve(struct asm_graph_t *g, int x, int partition_index, khash_t(big_kmer_count) *table, int *partition)
{
	if (partition[x] != partition_index ){
		return;
	}
	resolve_using_big_kmer(g, x, table);
//	resolve_using_pair_kmer(g, x, table);
	partition[x] = 0;

	int source = g->edges[x].source;
	int source_rc = g->nodes[source].rc_id;
	int target = g->edges[x].target;

	for (int i = 0; i < g->nodes[source_rc].deg; i++) {
		int i_e_rc = g->nodes[source_rc].adj[i];
		dfs_resolve(g, i_e_rc, partition_index, table, partition);
	}
	for (int i = 0; i < g->nodes[target].deg; i++) {
		int i_e = g->nodes[target].adj[i];
		dfs_resolve(g, i_e, partition_index, table, partition);
	}
}

void partition_graph(struct read_path_t *ori_read, struct asm_graph_t *g, int *partition, int n_threads,
					 int mmem,  int *n_partition)
{
	khash_t(bcpos) *dict = kh_init(bcpos);
	construct_read_index(ori_read, dict);

	for (int i = 0; i < (int) g->n_e; i++)
		partition[i] = -1;
	int count = 0;
	for (int i = 0; i < (int) g->n_e; i++)
		if (g->edges[i].seq_len < CONTIG_PARTITION_LEN) {
			log_warn("Dfs from %d", i);
			khash_t(union_barcode) *bc = kh_init(union_barcode);
			khash_t(big_kmer_count) *kmer_count_table = kh_init(big_kmer_count);
			int total_barcodes = 0, total_length_component = 0, count_resolvable=0;
			int n_edges = dfs_partition(g, partition, i, count, bc, &total_barcodes, &total_length_component, &count_resolvable);
			count++;

			log_warn("n edges %d n barcodes %d total length %d", n_edges, total_barcodes, total_length_component);
			if (total_barcodes > 10000) {
				//todo 1 work for big case
				continue;
			}
			if (total_barcodes > 10 && total_length_component/2 > MIN_COMPONENT && count_resolvable > 2) {
				filter_read_build_kmer(ori_read, dict, bc, kmer_count_table, 55, 70, n_threads,
									   mmem, count); //todo 0 choose range better
				dfs_resolve(g, i, count - 1, kmer_count_table, partition);
				//todo 0 not break
				break;
			}
		}
	*n_partition = count;
}

/**
 *
 * @return
 */
int resolve_using_big_kmer(struct asm_graph_t *g, int i_e,
						   khash_t(big_kmer_count) *table)
{
	if (!is_case_2_1_2(g, i_e))  {
		return -1;
	}
	int source = g->edges[i_e].source;
	int target = g->edges[i_e].target;
	int i_e_rc = g->edges[i_e].rc_id;
	int i_a0 = g->nodes[g->nodes[source].rc_id].adj[0];
	int i_a1 = g->nodes[g->nodes[source].rc_id].adj[1];
	int i_o0 = g->nodes[target].adj[0];
	int i_o1 = g->nodes[target].adj[1];
	struct asm_edge_t *a0 = &g->edges[i_a0];
	struct asm_edge_t *a1 = &g->edges[i_a1];
	struct asm_edge_t *o0 = &g->edges[i_o0];
	struct asm_edge_t *o1 = &g->edges[i_o1];
	struct asm_edge_t *e = &g->edges[i_e];

	assert(table != NULL);

	log_warn("table size %d", table->size);
	int c00 = get_new_seq_count(a0, o0, e, g->ksize, table);
	int c01 = get_new_seq_count(a0, o1, e, g->ksize, table);
	int c10 = get_new_seq_count(a1, o0, e, g->ksize, table);
	int c11 = get_new_seq_count(a1, o1, e, g->ksize, table);
	log_warn("count00 %d count01 %d count10 %d count11 %d", c00, c01,c10, c11);

	if (c00 > 0 && c11 > 0 && c00 + c11 > c10 + c01) {
		//todo 1 get new count proportion to path
		asm_join_edge3_wrapper(g, i_a0, i_e, i_o0, g->edges[i_e].count / 2);
		asm_join_edge3_wrapper(g, i_a1, i_e, i_o1, g->edges[i_e].count / 2);
		asm_remove_edge(g, i_e);
		asm_remove_edge(g, i_e_rc);
	} else if (c10 > 0 && c01 > 0 && c10 + c01 > c00 + c11) {
		asm_join_edge3_wrapper(g, i_a0, i_e, i_o1, g->edges[i_e].count / 2);
		asm_join_edge3_wrapper(g, i_a1, i_e, i_o0, g->edges[i_e].count / 2);
		asm_remove_edge(g, i_e);
		asm_remove_edge(g, i_e_rc);
	} else {
		return -1;
	}
	return 0;
}
