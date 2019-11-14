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
#include "build_hash_table.h"


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
	gint_t hash_key = MurmurHash3_x64_64(new_seq, (new_len + 3) >> 2);
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
	int to_char[4] = {65, 67, 71, 84};
	for (int i = 0; i < len; i++) {
		printf("%c", to_char[new_seq[i >> 2] >> ((i & 3) << 1) & 3]);
	}
	printf("\n");
}

int join_3_and_count(struct asm_edge_t *left, struct asm_edge_t *right, struct asm_edge_t *mid,
                     int big_ksize, int graph_ksize, khash_t(pair_kmer_count) *table)
{
	//todo 1 get more pair instead of only one pair right now
	int span_len = big_ksize;
	assert(mid->seq_len <= span_len - 2);
	int mid_len = mid->seq_len;
	int left_len = MIN(left->seq_len - graph_ksize, span_len - mid_len - 1);
	int right_len = MIN(right->seq_len - graph_ksize, span_len - mid_len - 1);
	int new_len = left_len + mid_len + right_len;
	if (new_len < span_len)
		return -1;

	uint8_t *new_seq = calloc((new_len + 3) >> 2, 1);
	copy_seq32_seq8(left->seq, left->seq_len - graph_ksize - left_len, new_seq, 0, left_len);
	copy_seq32_seq8(mid->seq, 0, new_seq, left_len, mid->seq_len);
	copy_seq32_seq8(right->seq, graph_ksize, new_seq, left_len + mid_len, right_len);
	print_u8_seq(new_seq, new_len);

	//todo 1: does not need to copy all middle seq

	int value = 0;
	for (int i = 0; i < new_len - span_len + 1; i++) {
		uint8_t *h1 = calloc((big_ksize + 3) >> 2, sizeof(uint8_t));
		get_seq(new_seq, i, big_ksize, h1);
//		print_u8_seq(h1, big_ksize);
		int64_t key;
		key = MurmurHash3_x64_64(h1, (big_ksize + 3) >> 2);

		khint_t k = kh_get(pair_kmer_count, table, key);
		if (k == kh_end(table))
			continue;
		int value1 = kh_value(table, k);
//		log_warn("value of get pair is %d", value1);
		value += value1;
	}
	return value;
}

int join_2_and_count(struct asm_edge_t *a, struct asm_edge_t *b,
                     int big_ksize, int graph_ksize, khash_t(pair_kmer_count) *table)
{
	int span_len = big_ksize;

	int a_len = MIN(a->seq_len - graph_ksize, span_len - 1);
	int b_len = MIN(b->seq_len, span_len - 1);
	int new_len = a_len + b_len;
	log_warn("len a %d len b %d len new %d", a->seq_len, b->seq_len, new_len);
	if (new_len < span_len)
		return -1;

	uint8_t *new_seq = calloc((new_len + 3) >> 2, 1);
	copy_seq32_seq8(a->seq, a->seq_len - graph_ksize - a_len, new_seq, 0, a_len);
	copy_seq32_seq8(b->seq, 0, new_seq, a_len, b_len);
	print_u8_seq(a->seq, a->seq_len);
	print_u8_seq(b->seq, b->seq_len);
	print_u8_seq(new_seq, new_len);
	printf("\n\n");
	int value = 0;
	for (int i = 0; i < new_len - span_len + 1; i++) {
		uint8_t *h1 = calloc((big_ksize + 3) >> 2, sizeof(uint8_t));
		get_seq(new_seq, i, big_ksize, h1);
		print_u8_seq(h1, big_ksize);
		int64_t key;
		key = MurmurHash3_x64_64(h1, (big_ksize + 3) >> 2);

		khint_t k = kh_get(pair_kmer_count, table, key);
		if (k == kh_end(table))
			continue;
		int value1 = kh_value(table, k);
//		log_warn("value of get pair is %d", value1);
		value += value1;
	}
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

//	log_warn("len of: a0 %d a1 %d e %d o0 %d o1%d\n", a0->seq_len, a1->seq_len, e->seq_len, o0->seq_len, o1->seq_len);
	if (a0->rc_id == i_o0 || a0->rc_id == i_o1 || a1->rc_id == i_o0 || a1->rc_id == i_o1)
		return 0;
	if (a0->rc_id == i_a1 || a1->rc_id == i_a0)
		return 0;
	return 1;
}

int dfs_partition(struct asm_graph_t *g, int *partition, int i_e, int index_par,
                  struct partition_information *part_infor)
{
//	log_warn("dfs from %d %d", i_e, partition[i_e]);
	if (partition[i_e]) {
		if (partition[i_e] != index_par) {
			log_error("dfs partition[i_e] %d, index :%d", partition[i_e], index_par);
		}
		return 0;
	}
	if (g->edges[i_e].seq_len >= CONTIG_PARTITION_LEN) {
//		assert(g->edges[i_e].barcodes[0].n_item <= g->edges[i_e].barcodes[2].n_item);
//		assert(g->edges[i_e].barcodes[0].n_item <= g->edges[i_e].barcodes[1].n_item);
//		*total_barcodes += append_barcode_contig(&g->edges[i_e].barcodes[2], uni_bar);//todo 1 all barcode
		return 0;
	}
	part_infor->total_len += g->edges[i_e].seq_len;
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
		res += dfs_partition(g, partition, i_e_rc, index_par, part_infor);
	}
	for (int i = 0; i < g->nodes[target].deg; i++) {
		int i_e = g->nodes[target].adj[i];
		res += dfs_partition(g, partition, i_e, index_par, part_infor);
	}
	return res;
}

void get_pos(khash_t(bcpos) *dict, khash_t(union_barcode) *bc, struct read_index_t **read_pos, int *n_barcodes)
{
	log_warn("start get pos");
	int count = 0;
	struct read_index_t *pos = NULL;
	for (khint_t it = kh_begin(bc); it != kh_end(bc); it++) {
		if (!(kh_exist(bc, it))) {
			continue;
		}

		gint_t key = kh_key(bc, it);
		khiter_t k = kh_get(bcpos, dict, key);
		assert(k != KMHASH_END(dict));
		pos = realloc(pos, (count + 1) * sizeof(struct read_index_t));
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
	pthread_mutex_init(&km_count_bundle->lock, NULL);
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
	atomic_add_and_fetch32(&kh_value(table, h), count);
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
	int n_barcodes = 0;
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

void dfs_resolve(struct asm_graph_t *g, int i_e, int partition_index, khash_t(pair_kmer_count) *table, int *partition,
                 int *resolve_stat)
{
//	log_warn("resol %d", i_e);
	if (g->edges[i_e].seq_len >= CONTIG_PARTITION_LEN) {
		return 0;
	}
	if (partition[i_e] != partition_index && partition[i_e] != -partition_index) {
		log_error("partition[i_e] %d, real :%d", partition[i_e], partition_index);
		assert(0);
	}
	if (partition[i_e] == -partition_index || g->edges[i_e].source == -1) {
		return;
	}
////	int t = 1;
//	if (t == NOT_212_CASE) {
//		t = resolve_202_using_big_kmer(g, i_e, table);
//		if (t == NOT_212_CASE) {
//			t = resolve_202_using_big_kmer(g, g->edges[i_e].rc_id, table);
//		}
//		resolve_stat[t]++;
//	}
	int x_rc = g->edges[i_e].rc_id;
	partition[i_e] = -partition[i_e];
	partition[x_rc] = -partition[x_rc];

	int source = g->edges[i_e].source;
	int source_rc = g->nodes[source].rc_id;
	int target = g->edges[i_e].target;

	for (int i = 0; i < g->nodes[source_rc].deg; i++) {
		int i_e_rc = g->nodes[source_rc].adj[i];
		dfs_resolve(g, i_e_rc, partition_index, table, partition, resolve_stat);
	}
	for (int i = 0; i < g->nodes[target].deg; i++) {
		int i_e = g->nodes[target].adj[i];
		dfs_resolve(g, i_e, partition_index, table, partition, resolve_stat);
	}
}

void partition_graph(struct read_path_t *ori_read, struct asm_graph_t *g, int *partition, int n_threads,
                     int mmem, int *n_partition, khash_t(pair_kmer_count) *kmer_pair_table)
{
//	khash_t(bcpos) *dict = kh_init(bcpos);
//	construct_read_index(ori_read, dict);

	for (int i = 0; i < (int) g->n_e; i++)
		partition[i] = 0;
	int count = 0;
	int *resolve_stat = calloc(4, 4);
	for (int i = 0; i < (int) g->n_e; i++)
		if (partition[i] == 0)
			if (g->edges[i].seq_len < CONTIG_PARTITION_LEN) {
				log_warn("Dfs from %d", i);
				count++;
//			khash_t(big_kmer_count) *kmer_count_table = kh_init(big_kmer_count);
				struct partition_information *part_infor
					= calloc(1, sizeof(struct partition_information));
				part_infor->union_barcodes = kh_init(union_barcode);
				int n_edges = dfs_partition(g, partition, i, count, part_infor);
//
//			log_warn("n edges %d n barcodes %d total length %d", n_edges, total_barcodes, total_length_component);
//			if (total_barcodes > 10000) {
//				//todo 1 work for big case
//				continue;
//			}
//			if (total_barcodes > 10 && total_length_component / 2 > MIN_COMPONENT && count_resolvable > 2) {
////				filter_read_build_kmer(ori_read, dict, bc, kmer_count_table, 55, 70, n_threads,
////									   mmem, count); //todo 0 choose range better
				dfs_resolve(g, i, count, kmer_pair_table, partition, resolve_stat);
				//todo 0 not break
//				break;
//			}
			}
//	for (int i = 0; i < g->n_e; i++) {
////		if (g->edges[i].seq_len < CONTIG_PARTITION_LEN)
//		assert(partition[i] <= 0);
//	}
	for (int i = 0 ; i < g->n_e; i++) {
		int t = resolve_212_using_big_kmer(g, i, kmer_pair_table);
		resolve_stat[t]++;
		if (t == NOT_212_CASE) {
			t = resolve_202_using_big_kmer(g, i, kmer_pair_table);
			resolve_stat[t]++;
		}
	}
	log_info("Resolve 2-2 case:");
	for (int i = 0; i < 4; i++) {
		if (i == NOT_HAVE_SPAN_KMER) {
			log_info("   Can't solve - not have span kmer: %d", resolve_stat[i]);
		}
		if (i == NOT_LONG_ENOUGH) {
			log_info("   Can't solve - kmer not long enough: %d", resolve_stat[i]);
		}
		if (i == 0) {
			log_info("   Solved: %d", resolve_stat[i]);
		}
	}
	*n_partition = count;
}

int resolve_212_using_big_kmer(struct asm_graph_t *g, int i_e,
                               khash_t(pair_kmer_count) *table)
{
	if (!is_case_2_1_2(g, i_e)) {
		return NOT_212_CASE;
	}
	int source = g->edges[i_e].source;
	int target = g->edges[i_e].target;
	int i_e_rc = g->edges[i_e].rc_id;
	int i_a0 = g->edges[g->nodes[g->nodes[source].rc_id].adj[0]].rc_id;
	int i_a1 = g->edges[g->nodes[g->nodes[source].rc_id].adj[1]].rc_id;
	int i_o0 = g->nodes[target].adj[0];
	int i_o1 = g->nodes[target].adj[1];
	struct asm_edge_t *a0 = &g->edges[i_a0];
	struct asm_edge_t *a1 = &g->edges[i_a1];
	struct asm_edge_t *o0 = &g->edges[i_o0];
	struct asm_edge_t *o1 = &g->edges[i_o1];
	struct asm_edge_t *e = &g->edges[i_e];

	assert(table != NULL);

	if (e->seq_len > DISTANCE_KMER + KMER_PAIR_SIZE - 2)
		return NOT_LONG_ENOUGH;
	int pair_ksize = KMER_PAIR_SIZE;
	int c00 = join_3_and_count(a0, o0, e, BIG_KSIZE, g->ksize, table);
	int c01 = join_3_and_count(a0, o1, e, BIG_KSIZE, g->ksize, table);
	int c10 = join_3_and_count(a1, o0, e, BIG_KSIZE, g->ksize, table);
	int c11 = join_3_and_count(a1, o1, e, BIG_KSIZE, g->ksize, table);
	log_warn("pair_ksize %d count00 %d count01 %d count10 %d count11 %d", pair_ksize, c00, c01, c10, c11);

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
		return NOT_HAVE_SPAN_KMER;
	}
	return 0;
}

/**
 *                   o
 *                   |
 *                   |  a0
 *                   |
 *                   v
 *                   o
 *                 /  \
 *                /    \
 *       e       /      \     a1
 *  o---------->o        o<-----------o
 *               \      /
 *                \    /
 *                 \  /
 *                  o
 *                  ^
 *                  |
 *                  | a2
 *                  |
 *                  |
 *                  o
 *
 */
int resolve_202_using_big_kmer(struct asm_graph_t *g, int i_e,
                               khash_t(pair_kmer_count) *table)
{
	int i_target = g->edges[i_e].target;
	if (i_target == -1)
		return NOT_202_CASE;
	int i_e_rc = g->edges[i_e].rc_id;
	struct asm_node_t *target = &g->nodes[g->edges[i_e].target];
	if (target->deg != 2)
		return NOT_202_CASE;
	int i_a0 = g->edges[target->adj[0]].rc_id;
	int i_a2 = g->edges[target->adj[1]].rc_id;

	struct asm_node_t *target_a0 = &g->nodes[g->edges[i_a0].target];
	struct asm_node_t *target_a2 = &g->nodes[g->edges[i_a2].target];
	if (target_a0->deg != 2 ||
	    target_a2->deg != 2)
		return NOT_202_CASE;

	int u = target_a0->adj[0];
	int v = target_a0->adj[1];
	if (u != i_e_rc && v != i_e_rc)
		return NOT_202_CASE;

	int i_a1_rc = u ^v ^i_e_rc;

	int x = target_a2->adj[0];
	int y = target_a2->adj[1];
	if (x != i_e_rc && y != i_e_rc)
		return NOT_202_CASE;

	if (i_e_rc ^ x ^ y ^ i_a1_rc)
		return NOT_202_CASE;

	gint_t i_a1 = g->edges[i_a1_rc].rc_id;
	if (i_e == i_a0 || i_e == i_a2 || i_e == i_a1_rc || i_a0 == i_a1 || i_a0 == i_a2 || i_a1 == i_a2)
		return NOT_202_CASE;

	gint_t i_a0_rc = g->edges[i_a0].rc_id;
	gint_t i_a2_rc = g->edges[i_a2].rc_id;
	struct asm_edge_t *a0_rc = &g->edges[i_a0_rc];
	struct asm_edge_t *a2_rc = &g->edges[i_a2_rc];
	struct asm_edge_t *e = &g->edges[i_e];
	struct asm_edge_t *a1 = &g->edges[i_a1];
	log_warn("edge len %d %d %d %d", e->seq_len, a0_rc->seq_len, a1->seq_len, a2_rc->seq_len);
	assert(a0_rc->source != -1);
	assert(a1->source != -1);
	assert(a2_rc->source != -1);
	assert(e->source != -1);
	int ce0 = join_2_and_count(e, a0_rc, BIG_KSIZE, g->ksize, table);
	int ce2 = join_2_and_count(e, a2_rc, BIG_KSIZE, g->ksize, table);
	int c10 = join_2_and_count(a1, a0_rc, BIG_KSIZE, g->ksize, table);
	int c12 = join_2_and_count(a1, a2_rc, BIG_KSIZE, g->ksize, table);


	log_warn("%d %d %d %d", ce0, ce2, c10, c12);
	if (ce0 > 0 && c12 > 0 && ce0 + c12 > ce2 + c10) {
		asm_join_edge_wrapper(g, i_e, i_a0_rc);
		asm_join_edge_wrapper(g, i_a1, i_a2_rc);
	} else if (ce2 > 0 && c10 > 0 && ce2 + c10 > ce0 + c12) {
		asm_join_edge_wrapper(g, i_e, i_a2_rc);
		asm_join_edge_wrapper(g, i_a1, i_a0_rc);
	} else {
		return NOT_HAVE_SPAN_KMER;
	}
	return 0;
}

void resolve_1_2(struct asm_graph_t *g, struct opt_proc_t *opt)
{
	int *partition = calloc(g->n_e, 4);
	khash_t(big_kmer_count) **partition_kmer_count = NULL;
	int n_partitions;
	struct read_path_t *ori_read = calloc(1, sizeof(struct read_path_t));
	assert(opt->n_files == 1); //todo 1 handle multiple files
	log_info("%s", opt->files_1[0]);
	ori_read->R1_path = opt->files_1[0];
	ori_read->R2_path = opt->files_2[0];
	ori_read->idx_path = opt->files_I[0];

	khash_t(pair_kmer_count) *table = kh_init(pair_kmer_count);
	build_pair_kmer_table(opt, table);

	partition_graph(ori_read, g, partition, opt->n_threads, opt->mmem, &n_partitions, table);
	free(partition);
}
