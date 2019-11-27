//
// Created by BioTuring on 2019-09-25.
//

#include <zlib.h>
#include <string.h>

#include "smart_load.h"
#include "sort_read.h"
#include "verbose.h"
#include "../radix_sort.h"
#include "../buffer_file_wrapper.h"
#include "../io_utils.h"
#include "kseq.h"
#include "../utils.h"
#include "../log.h"
#include "count_barcodes.h"
#include "minimizers.h"
#include "../assembly_graph.h"
#include "get_buffer.h"
static const char *bit_rep[16] = {
	[ 0] = "0000", [ 1] = "0001", [ 2] = "0010", [ 3] = "0011",
	[ 4] = "0100", [ 5] = "0101", [ 6] = "0110", [ 7] = "0111",
	[ 8] = "1000", [ 9] = "1001", [10] = "1010", [11] = "1011",
	[12] = "1100", [13] = "1101", [14] = "1110", [15] = "1111",
};

#define PRINT_BYTE(byte) printf("%s%s\n", bit_rep[(byte >> 4)& 0x0F], bit_rep[byte & 0x0F]);
#define PRINT_UINT64(byte) printf("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n", bit_rep[(byte >> 60) & 0x0F], bit_rep[(byte >> 56)& 0x0F], \
								bit_rep[(byte >> 52) & 0x0F], bit_rep[(byte >> 48) & 0x0F], \
								bit_rep[(byte >> 44) & 0x0F], bit_rep[(byte >> 40) & 0x0F], \
								bit_rep[(byte >> 36) & 0x0F], bit_rep[(byte >> 32) & 0x0F], \
								bit_rep[(byte >> 28) & 0x0F], bit_rep[(byte >> 24) & 0x0F], \
								bit_rep[(byte >> 20) & 0x0F], bit_rep[(byte >> 16) & 0x0F], \
								bit_rep[(byte >> 12) & 0x0F], bit_rep[(byte >> 8) & 0x0F],  \
								bit_rep[(byte >> 4) & 0x0F], bit_rep[(uint8_t)byte & 0x0F]);

/**
 * Construct the dictionary index of
 * barcode position in the barcode sorted file
 * @param rpath pointer to the read_path_t struct, which contains path to R1 and R2
 * @param h khash struct that whose barcodes indices. Will be used to look for position of any barcode
 */
void smart_construct_read_index(struct read_path_t *rpath, khash_t(bcpos) *h)
{
	FILE *fp = xfopen(rpath->idx_path, "rb");
	size_t byte_read;
	khint_t k;
	int ret;
	uint64_t barcode;
	char *buf = alloca(40);
	while ((byte_read = fread(buf, 1, 40, fp))) {
		if (byte_read != 40)
			__ERROR("Corrupted barcode in read index file");
		barcode = unpack_int64((uint8_t *)buf);
		k = kh_put(bcpos, h, barcode, &ret);
		if (ret != 1)
			__ERROR("Insert barcode failed");
		kh_value(h, k).r1_offset = unpack_int64((uint8_t *)buf + 8);
		kh_value(h, k).r2_offset = unpack_int64((uint8_t *)buf + 16);
		kh_value(h, k).r1_len = unpack_int64((uint8_t *)buf + 24);
		kh_value(h, k).r2_len = unpack_int64((uint8_t *)buf + 32);
	}
	fclose(fp);
}

/**
 * @brief Reads of one set of barcodes into string stream
 * @param ref read_path struct that contain path to original read1 and read2 (maybe barcode included)
 * @param dict khash_t(bcpost) hash table to lookup the offset of one barcode in the sorted data
 * @param shared array of barcodes pseudo hash values
 * @param n_shared number of barcodes in shared array
 * @param buf1 buffer stream for read 1
 * @param buf2 buffer stream for read 2
 * @param size1 size of the buffer stream for read 1
 * @param size2 size of the buffer stream for read 2
 */
void stream_filter_read(struct read_path_t *ref, khash_t(bcpos) *dict,
                 uint64_t *shared, int n_shared, char **buf1, char **buf2,
                 uint64_t *size1, uint64_t *size2)
{
	struct read_index_t *pos;
	pos = malloc(n_shared * sizeof(struct read_index_t));
	khiter_t k;
	int i;
	for (i = 0; i < n_shared; ++i) {
		k = kh_get(bcpos, dict, shared[i]);
		pos[i] = kh_value(dict, k);
	}
	rs_sort(read_index, pos, pos + n_shared); //TODO: get rid of this function
	FILE *fi1 = xfopen(ref->R1_path, "rb");
	FILE *fi2 = xfopen(ref->R2_path, "rb");
	uint64_t m_buf1 = 0, m_buf2 = 0;
	uint64_t prev1 = 0 , prev2 = 0;

	for (i = 0 ; i < n_shared; ++i) {
		m_buf1 += pos[i].r1_len;
		m_buf2 += pos[i].r2_len;
	}
	char *buf1_ = calloc(m_buf1 + 1, sizeof(char));
	char *buf2_ = calloc(m_buf2 + 1, sizeof(char));

	for (i = 0; i < n_shared; ++i) {
		fseek(fi1, pos[i].r1_offset, SEEK_SET);
		xfread(buf1_ + prev1, 1, pos[i].r1_len, fi1);
		fseek(fi2, pos[i].r2_offset, SEEK_SET);
		xfread(buf2_ + prev2, 1, pos[i].r2_len, fi2);
		prev1 += pos[i].r1_len ;
		prev2 += pos[i].r2_len ;
	}
	fclose(fi1);
	fclose(fi2);
	free(pos);
	*buf1 = buf1_;
	*buf2 = buf2_;
	*size1 = m_buf1;
	*size2 = m_buf2;
}
#define kh_set(kname, hash, key, val) ({int ret; k = kh_put(kname, hash,key,&ret); kh_value(hash,k) = val; ret;})
// shorthand way to get the key from hashtable or defVal if not found
#define kh_get_val(kname, hash, key, defVal) ({k=kh_get(kname, hash, key);(k!=kh_end(hash)?kh_val(hash,k):defVal);})
/**
 * @brief Seeking and loading the barcode sequences into the buffer stream
 * @param opt
 */
void smart_load_barcode(struct opt_proc_t *opt)
{
	struct read_path_t read_sorted_path;
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
		log_info("Reads are not sorted. Sort reads by barcode sequence...");
		sort_read(opt, &read_sorted_path);
	}
	uint64_t bx_encoded = barcode_hash_mini(opt->bx_str);
	log_info("Hashed barcode: %lu", bx_encoded);
	uint64_t bx[1] = {bx_encoded}; //43 15 mock barcode pseudo hash id here
	khash_t(bcpos) *bx_pos_dict = kh_init(bcpos);
	smart_construct_read_index(&read_sorted_path, bx_pos_dict); //load the barcode indices

	khint_t k_bx = kh_get(bcpos, bx_pos_dict, bx_encoded);          // query the hash table
	if (k_bx == kh_end(bx_pos_dict)) {
		log_error("Barcode does not exists!");
	} else {
		log_info("Barcode does exist. Getting reads");
	}

	char *buf1, *buf2;
	uint64_t m_buf1, m_buf2;
	stream_filter_read(&read_sorted_path, bx_pos_dict, bx, 1, &buf1, &buf2, &m_buf1, &m_buf2);
	printf("buf1: %s, buf2: %s", buf1, buf2);
	struct read_t r1, r2;
	int pos1 = 0, pos2 = 0;
	int n_reads = 0;
	struct mm_db_t *db1, *db2;
	struct mm_hits_t *hits1, *hits2, *hits;


	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	struct mm_db_edge_t *mm_edges_db = mm_index_edges(&g, MINIMIZERS_KMER, MINIMIZERS_WINDOW);
	kh_mm_pw_t *kh_pw = kh_init(mm_pw);
	khiter_t k;

	while (get_read_from_fq(&r1, buf1, &pos1) == READ_SUCCESS && get_read_from_fq(&r2, buf2, &pos2) == READ_SUCCESS ) {
		n_reads++;
		db1 = mm_index_char_str(r1.seq, MINIMIZERS_KMER, MINIMIZERS_WINDOW, r1.len);
		db2 = mm_index_char_str(r2.seq, MINIMIZERS_KMER, MINIMIZERS_WINDOW, r2.len);
		khiter_t k1, k2;

		hits1 = mm_hits_init();
		hits2 = mm_hits_init();

		mm_hits_cmp(db1, mm_edges_db, hits1, &g);
		mm_hits_cmp(db2, mm_edges_db, hits2, &g);

		log_info("Number of hits on R1: %d", hits1->n);
		log_info("Number of hits on R2: %d", hits2->n);
		for (k1 = kh_begin(hits1->edges); k1 != kh_end(hits1->edges); ++k1){
			if (!kh_exist(hits1->edges, k1))
				continue;
			for (k2 = kh_begin(hits2->edges); k2 != kh_end(hits2->edges); ++k2) {
				if (!kh_exist(hits2->edges, k2))
					continue;
				uint64_t key1 = kh_key(hits1->edges, k1);
				uint64_t key2 = kh_key(hits2->edges, k2);
				uint64_t hk = key1 * 100000000 + key2;
				if (kh_get_val(mm_pw, kh_pw, hk, -1) == -1)
					kh_set(mm_pw, kh_pw, hk, 1);
				else
					kh_set(mm_pw, kh_pw, hk, kh_get_val(mm_pw, kh_pw, hk, -1) + 1);

			}
		}
		mm_hits_destroy(hits1);
		mm_hits_destroy(hits2);
		mm_db_destroy(db1);
		mm_db_destroy(db2);
	}
	hits = mm_hits_init();
	log_info("Number of read-pairs in barcode %s: %d", opt->bx_str, n_reads);
	for (k = kh_begin(kh_pw); k != kh_end(kh_pw); ++k){
		if (kh_exist(kh_pw, k)){
			uint64_t key = kh_key(kh_pw, k);
			uint64_t u = key / 100000000;
			uint64_t v = key % 100000000;
			log_info("Read pair that mapped to u->v: %lu -> %lu: %d", u, v, kh_val(kh_pw, k));
			if (kh_val(kh_pw, k) > 1) {
				kh_set(mm_edges, hits->edges, u, 1);
				kh_set(mm_edges, hits->edges, v, 1);
			}
		}

	}

	free(buf1);
	free(buf2);
}

struct mm_hits_t *get_hits_from_barcode(char *bc, struct bc_hit_bundle_t *bc_hit_bundle)
{
	struct read_path_t *read_sorted_path = bc_hit_bundle->read_sorted_path;
	khash_t(bcpos) *bx_pos_dict = bc_hit_bundle->bc_pos_dict;
	struct asm_graph_t *g = bc_hit_bundle->g;
	struct mm_db_edge_t *mm_edges_db = bc_hit_bundle->mm_edges;

	uint64_t bx_encoded = barcode_hash_mini(bc);
	uint64_t bx[1] = {bx_encoded}; //43 15 mock barcode pseudo hash id here

	khint_t k_bx = kh_get(bcpos, bx_pos_dict, bx_encoded);          // query the hash table
	if (k_bx == kh_end(bx_pos_dict)) {
		log_error("Barcode does not exists!");
	}

	char *buf1, *buf2;
	uint64_t m_buf1, m_buf2;
	stream_filter_read(read_sorted_path, bx_pos_dict, bx, 1, &buf1, &buf2,
			&m_buf1, &m_buf2);

	struct read_t r1, r2;
	int pos1 = 0, pos2 = 0;
	int n_reads = 0;
	struct mm_db_t *db1, *db2;
	struct mm_hits_t *hits1, *hits2, *hits;

	kh_mm_pw_t *kh_pw = kh_init(mm_pw);
	khiter_t k;

	while (get_read_from_fq(&r1, buf1, &pos1) == READ_SUCCESS && get_read_from_fq(&r2, buf2, &pos2) == READ_SUCCESS ) {
		n_reads++;
		db1 = mm_index_char_str(r1.seq, MINIMIZERS_KMER, MINIMIZERS_WINDOW, r1.len);
		db2 = mm_index_char_str(r2.seq, MINIMIZERS_KMER, MINIMIZERS_WINDOW, r2.len);
		khiter_t k1, k2;

		hits1 = mm_hits_init();
		hits2 = mm_hits_init();

		mm_hits_cmp(db1, mm_edges_db, hits1, g);
		mm_hits_cmp(db2, mm_edges_db, hits2, g);

		log_info("Number of hits on R1: %d", hits1->n);
		log_info("Number of hits on R2: %d", hits2->n);
		for (k1 = kh_begin(hits1->edges); k1 != kh_end(hits1->edges); ++k1){
			if (!kh_exist(hits1->edges, k1))
				continue;
			for (k2 = kh_begin(hits2->edges); k2 != kh_end(hits2->edges); ++k2) {
				if (!kh_exist(hits2->edges, k2))
					continue;
				uint64_t key1 = kh_key(hits1->edges, k1);
				uint64_t key2 = kh_key(hits2->edges, k2);
				uint64_t hk = key1 * 100000000 + key2;
				if (kh_get_val(mm_pw, kh_pw, hk, -1) == -1)
					kh_set(mm_pw, kh_pw, hk, 1);
				else
					kh_set(mm_pw, kh_pw, hk, kh_get_val(mm_pw, kh_pw, hk, -1) + 1);

			}
		}
		mm_hits_destroy(hits1);
		mm_hits_destroy(hits2);
		mm_db_destroy(db1);
		mm_db_destroy(db2);
	}
	hits = mm_hits_init();
	for (k = kh_begin(kh_pw); k != kh_end(kh_pw); ++k){
		if (kh_exist(kh_pw, k)){
			uint64_t key = kh_key(kh_pw, k);
			uint64_t u = key / 100000000;
			uint64_t v = key % 100000000;
			log_info("Read pair that mapped to u->v: %lu -> %lu: %d", u, v, kh_val(kh_pw, k));
			if (kh_val(kh_pw, k) > 1) {
				kh_set(mm_edges, hits->edges, u, 1);
				kh_set(mm_edges, hits->edges, v, 1);
			}
		}
	}
	free(buf1);
	free(buf2);
	return hits;
}

void get_bc_hit_bundle(struct opt_proc_t *opt, struct bc_hit_bundle_t *bc_hit_bundle)
{
	struct read_path_t *read_sorted_path = calloc(1, sizeof(struct read_path_t));
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path->R1_path = opt->files_1[0];
		read_sorted_path->R2_path = opt->files_2[0];
		read_sorted_path->idx_path = opt->files_I[0];
	} else {
		log_info("Reads are not sorted. Sort reads by barcode sequence...");
		sort_read(opt, read_sorted_path);
	}
	khash_t(bcpos) *bc_pos_dict = kh_init(bcpos);
	smart_construct_read_index(read_sorted_path, bc_pos_dict); //load the barcode indices

	struct asm_graph_t *g = calloc(1, sizeof(struct asm_graph_t));
	assert(opt->in_file != NULL);
	load_asm_graph(g, opt->in_file);
	struct mm_db_edge_t *mm_edges = mm_index_edges(g, MINIMIZERS_KMER,
			MINIMIZERS_WINDOW);

	bc_hit_bundle->g = g;
	bc_hit_bundle->read_sorted_path = read_sorted_path;
	bc_hit_bundle->bc_pos_dict = bc_pos_dict;
	bc_hit_bundle->mm_edges = mm_edges;
}

void bc_hit_bundle_destroy(struct bc_hit_bundle_t *bc_hit_bundle)
{
	asm_graph_destroy(bc_hit_bundle->g);
	free(bc_hit_bundle->g);
	free(bc_hit_bundle->read_sorted_path);
	kh_destroy(bcpos, bc_hit_bundle->bc_pos_dict);
	//TODO: destroy everything
}

