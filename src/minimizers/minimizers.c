//
// Created by BioTuring on 2019-11-08.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <inttypes.h>

#include "../assembly_graph.h"
#include "attribute.h"
#include "../utils.h"
#include "minimizers.h"
#include "count_barcodes.h"
#include "../log.h"
#include "atomic.h"
#include "fastq_producer.h"


#ifdef DEBUG
#define DEBUG_MM printf
#else
#define DEBUG_MM
#endif

#define MOLECULE_MARGIN 6000
#define MIN_READS_TO_HITS 10
#define MAX_READS_TO_HITS 88
#define EMPTY_BX UINT64_MAX

struct minimizer_bundle_t {
    struct dqueue_t *q;
    struct mini_hash_t **bx_table;
    struct mini_hash_t **rp_table;
    struct asm_grapgh_t *g;
    struct mm_db_edge_t *mm_edges;
};

const char *bit_rep[16] = {
	[ 0] = "0000", [ 1] = "0001", [ 2] = "0010", [ 3] = "0011",
	[ 4] = "0100", [ 5] = "0101", [ 6] = "0110", [ 7] = "0111",
	[ 8] = "1000", [ 9] = "1001", [10] = "1010", [11] = "1011",
	[12] = "1100", [13] = "1101", [14] = "1110", [15] = "1111",
};
// shorthand way to get the key from hashtable or defVal if not found
#define kh_get_val(kname, hash, key, defVal) ({k=kh_get(kname, hash, key);(k!=kh_end(hash)?kh_val(hash,k):defVal);})

// shorthand way to set value in hash with single line command.  Returns value
// returns 0=replaced existing item, 1=bucket empty (new key), 2-adding element previously deleted
#define kh_set(kname, hash, key, val) ({int ret; k = kh_put(kname, hash,key,&ret); kh_value(hash,k) = val; ret;})

#define PRINT_BYTE(byte) printf("%s%s\n", bit_rep[(byte >> 4)& 0x0F], bit_rep[byte & 0x0F]);
#define PRINT_UINT64(byte) printf("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n", bit_rep[(byte >> 60) & 0x0F], bit_rep[(byte >> 56)& 0x0F], \
								bit_rep[(byte >> 52) & 0x0F], bit_rep[(byte >> 48) & 0x0F], \
								bit_rep[(byte >> 44) & 0x0F], bit_rep[(byte >> 40) & 0x0F], \
								bit_rep[(byte >> 36) & 0x0F], bit_rep[(byte >> 32) & 0x0F], \
								bit_rep[(byte >> 28) & 0x0F], bit_rep[(byte >> 24) & 0x0F], \
								bit_rep[(byte >> 20) & 0x0F], bit_rep[(byte >> 16) & 0x0F], \
								bit_rep[(byte >> 12) & 0x0F], bit_rep[(byte >> 8) & 0x0F],  \
								bit_rep[(byte >> 4) & 0x0F], bit_rep[(uint8_t)byte & 0x0F]);


#define BARCODES100M 100663320
#define BIG_CONSTANT(x) (x##LLU)

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

inline uint64_t rotl64(uint64_t x, int8_t r)
{
	return (x << r) | (x >> (64 - r));
}

#define	ROTL32(x,y)	rotl32(x,y)
#define ROTL64(x,y)	rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)

#endif // !defined(_MSC_VER)

static inline uint64_t fmix64(uint64_t k)
{
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xff51afd7ed558ccd);
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
	k ^= k >> 33;

	return k;
}

static inline uint64_t getblock64(const uint64_t * p, int i)
{
	return p[i];
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

static inline char *uin32t2seq(uint32_t *s, int l)
{
	char *seq;
	int i;
	seq = calloc(l, sizeof(char));
	for (i = 0; i < l; ++i)
		seq[i] = nt4_char[__binseq_get(s, i)];
	return seq;
}

void mm_print(uint64_t mm, int k)
{
	int i, j, pad;
	uint8_t c;
	for (j = 0 ; j < k; ++j) {
		pad = (62 - 2*j); /* 64 bits/number and 2bits/char */
		c = (mm & ((uint64_t)3 << pad)) >> pad;
		printf("%c", nt4_char[c]);
	}
	printf("\n");
}

void mm_db_insert(struct mm_db_t *db, uint64_t km, uint32_t p)
{
	if (db->n == db->size) {
		db->mm =  realloc(db->mm, (db->size << 1) * sizeof(uint64_t));
		db->p =  realloc(db->p, (db->size << 1) * sizeof(uint64_t));
		db->size <<= 1;
	}
	db->mm[db->n] = km;
	db->p[db->n++] = p;
}

struct mm_align_t *init_mm_align()
{
	struct mm_align_t *aln = calloc(1, sizeof(struct mm_align_t));
	aln->pos = aln->edge = aln->cnt = 0;
	return aln;
}

/**
 * @brief Insert one hit of minimizers to the hits database
 * @param hits hits database
 * @param mm encode minimizers
 * @param e edge ID
 * @param p mapped position
 */
void mm_hits_insert(struct mm_hits_t *hits, uint64_t mm, uint64_t e, uint32_t p) {
	khiter_t k;
	int miss;
	struct mm_align_t *aln;
	k = kh_get(mm_edges, hits->edges, e);
	// Increase the number of hits of each edge to 1
	if (k == kh_end(hits->edges)) {
		k = kh_put(mm_edges, hits->edges, e, &miss);
		kh_value(hits->edges, k) = 1;
	} else
		kh_value(hits->edges, k)++;
	k = kh_get(mm_align, hits->aln, mm);
	// Add the alignment object of each minimizers hit
	if (k == kh_end(hits->aln)) {
		k = kh_put(mm_align, hits->aln, mm, &miss);
		kh_value(hits->aln, k) = init_mm_align();
	}
	aln = kh_value(hits->aln, k);
	aln->cnt++;
	aln->edge = e;
	aln->pos  = p; // Single-ton minimizers so the pos and edge ID doens't change
	aln->mm = mm;

	hits->n++;
}

struct mm_db_t * mm_db_init()
{
	struct mm_db_t *db = calloc(1, sizeof(struct mm_db_t));
	db->n = 0;
	db->size = 8;
	db->mm = calloc(db->size, sizeof(struct mm_db_t));
	db->p = calloc(db->size, sizeof(struct mm_db_t));
	return db;
}

void mm_db_destroy(struct mm_db_t *db)
{
	free(db->mm);
	free(db->p);
	free(db);
}

struct mm_hits_t *mm_hits_init()
{
	struct mm_hits_t *hits = calloc(1, sizeof(struct mm_hits_t));
	hits->edges = kh_init(mm_edges);
	hits->aln = kh_init(mm_align);
	hits->n = 0;
	return hits;
}

void mm_hits_destroy(struct mm_hits_t *hits)
{
	for (khiter_t it = kh_begin(hits->aln); it != kh_end(hits->aln); ++it){
		if (!kh_exist(hits->aln, it))
			continue;
		free(kh_val(hits->aln, it));
	}
	kh_destroy(mm_edges, hits->edges);
	kh_destroy(mm_align, hits->aln);
	free(hits);
}

static inline uint64_t get_km_i_bin(uint32_t *s, int i, int k)
{
	uint64_t km, c;
	int j;
	int pad = (32 - k - 1)*2;
	km = 0;
	for (j = 0; j < k; ++j) {
		c = (uint64_t)__binseq_get(s, i + j);
		km |= c;
		km <<= 2;
	}
	km <<= pad;
	return km;
}

static inline uint64_t get_km_i_str(char *s, int i, int k)
{
	uint64_t km, c;
	int j;
	int pad = (32 - k - 1)*2;
	km = 0;
	for (j = 0; j < k; ++j) {
		c = (uint64_t)nt4_table[s[i + j]];
		km |= c;
		km <<= 2;
	}
	km <<= pad;
	return km;
}

static inline char *mm_dump_seq(uint64_t s, uint32_t l)
{
	char *seq = calloc(l, sizeof(char));
	for (int i = 0; i < l; ++i) {
		uint8_t c = (uint8_t)((s >> (64 - i*2 - 2)) & (uint64_t)3);
		seq[i] = nt4_char[c];
	}
	return seq;
}

#define HASH64(k) MurmurHash3_x64_64((uint8_t *)&k, sizeof(uint64_t)/sizeof(uint8_t))
#define DEBUG_PRINT

struct mm_db_t * mm_index_bin_str(uint32_t *s, int k, int w, int l)
{
	struct mm_db_t *db = mm_db_init();
	db->k = k;

	int i, j, p = -1;
	uint64_t km, mm, c;
	uint64_t km_h, mm_h = 0;
	int pad = (32 - k - 1)*2;

	for (i = 0; i < l - w -k + 1; ++i) {
		DEBUG_PRINT("[i = %d]\n", i);
		if (p < i) {
			km = mm = get_km_i_bin(s, i, k);
			mm_h = km_h = HASH64(mm);
			p = i;
			for (j = 0; j < w; ++j) {
				c = (uint64_t) __binseq_get(s, i + j + k - 1);
				km |= ((uint64_t) c << (pad + 2));
				km_h = HASH64(km);
				if (km_h < mm_h) {
					mm = km;
					mm_h = km_h;
					p = i + j;
					mm_db_insert(db, km, p);
				}
				km <<= 2;
			}
			DEBUG_PRINT("[1]minimizers at window %d: %d\n", i, p);

			continue;
		} else {
			c = (uint64_t) __binseq_get(s, i + w + k - 2);
			km |= ((uint64_t) c << (pad + 2));
			km_h = HASH64(km);
			if (km_h < mm_h){
				p = i + w - 1;
				mm = km;
				mm_h = km_h;
				mm_db_insert(db, km, p);
				DEBUG_PRINT("[2]minimizers at window %d: %d\n", i, p);
			}
			km <<= 2;
		}
	}
	//mm_print(db);
	return db;
}

struct mm_db_t * mm_index_char_str(char *s, int k, int w, int l)
{
	struct mm_db_t *db = mm_db_init();
	db->k = k;

	int i, j, p = -1;
	uint64_t km, mm, c;
	uint64_t km_h, mm_h = 0;
	uint64_t *km2 = calloc(l, 8);
	int pad = (32 - k - 1)*2;

	km2[0] = get_km_i_str(s, 0, k);
	for(int i =  1 ; i < l-k+1; i++) {
		c = (uint64_t) nt4_table[s[i + k - 1]];
		km2[i] = km2[i-1]| (c << (pad +2));
	}
	for(int i =  1 ; i < l-k+1; i++) {
		km2[i] = HASH64(km2[i]);
	}

	for (i = 0; i < l - w -k + 1; ++i) {
		DEBUG_PRINT("[i = %d]\n", i);
		if (p < i) {
			p = i;
			for (j = 0; j < w; ++j) {
				if (km2[i+j] < mm_h) {
					mm_h = km2[i+j];
					p = i + j;
					mm_db_insert(db, km2[i+j], p);
				}
			}
			DEBUG_PRINT("[1]minimizers at window %d: %d\n", i, p);
		} else {
			if (km2[i+w-1] < mm_h){
				p = i + w - 1;
				mm_h = km2[i+w-1];
				mm_db_insert(db, km2[i+w-1], p);
				DEBUG_PRINT("[2]minimizers at window %d: %d\n", i, p);
			}
		}
	}
	//mm_print(db);
	return db;
}

struct mm_db_edge_t *mm_db_edge_init()
{
	struct mm_db_edge_t *db = calloc(1, sizeof(struct mm_db_edge_t));
	db->h = kh_init(mm_hash);
	db->cnt = kh_init(mm_hash);
	db->p = kh_init(mm_hash);
	return db;
}

/**
 * @brief Insert one minimizer into the minimizer database of edges.
 * @param db edges minimizers database collection
 * @param mm encoded minimizer
 * @param e edge ID
 * @param p mapped position
 */
void mm_db_edge_insert(struct mm_db_edge_t *db, uint64_t mm, uint32_t e, uint32_t p)
{
	khiter_t k, ki, k_h;

	ki = kh_get(mm_hash, db->cnt, mm);
	if (ki == kh_end(db->cnt)) {
		kh_set(mm_hash, db->cnt, mm, 1);
		kh_set(mm_hash, db->h, mm, e);
		kh_set(mm_hash, db->p, mm, p);
	} else {
		k_h = kh_get(mm_hash, db->h, mm);
		//if (e != kh_value(db->h, k_h)) // Allow multiple mapped on one edge minimizer to be single-ton
		kh_value(db->cnt, ki)++;
	}

}

void mm_db_edge_destroy(struct mm_db_edge_t *db)
{
	kh_destroy(mm_hash, db->h);
	kh_destroy(mm_hash, db->cnt);
	kh_destroy(mm_hash, db->p);
	free(db);
}

static inline void mm_singleton_stats(struct mm_db_edge_t *db, int k_size)
{
	khiter_t k_cnt, k_h, k;
	uint64_t single_cnt = 0, n_mm = 0;
	for (k_cnt = kh_begin(db->cnt); k_cnt != kh_end(db->cnt); ++k_cnt) {
		if (kh_exist(db->cnt, k_cnt)) {
			n_mm++;
			char *s = mm_dump_seq(kh_key(db->cnt, k_cnt), k_size);
			uint32_t p = kh_get_val(mm_hash, db->p, kh_key(db->cnt, k_cnt), -1) ;
			if (kh_value(db->cnt, k_cnt) == 1) {
				single_cnt++;
				DEBUG_MM("Singleton %s: contig: %d, pos: %d\n", s, kh_get_val(mm_hash, db->h, kh_key(db->cnt, k_cnt), -1), p);
			} else {
				DEBUG_MM("Redundant minimizers %s: contig: %d, pos: %d\n", s, kh_get_val(mm_hash, db->h, kh_key(db->cnt, k_cnt), -1), p);

			}
		}
	}

}

struct mm_db_edge_t *mm_index_edges(struct asm_graph_t *g, int k, int w) {
	int j, e;
	struct mm_db_edge_t *mm_db_e = mm_db_edge_init();
	struct mm_db_t *mm_db;

	for (e = 0; e < g->n_e; ++e) {
		if (!(e % 10000))
			log_info("%d edges", e);
		mm_db = mm_index_bin_str(g->edges[e].seq, k, w, g->edges[e].seq_len);
		for (j = 0; j < mm_db->n; ++j) {
			mm_db_edge_insert(mm_db_e, mm_db->mm[j], e, mm_db->p[j]);
		}
		mm_db_destroy(mm_db);
	}
	mm_singleton_stats(mm_db_e, k);
	return mm_db_e;
}

/**
 * @brief Map the minimizers from query to database
 * @param db minimizers collection from query
 * @param db_e minimizers collection from database
 * @param hits object to store the results
 * @return
 */
void *mm_hits_cmp(struct mm_db_t *db, struct mm_db_edge_t *db_e, struct mm_hits_t *hits, struct asm_graph_t *g)
{
	khiter_t k;
	uint32_t i, p;
	uint64_t e;
	int single_ton = 0;
	for (i = 0; i < db->n; ++i) {
		k = kh_get(mm_hash, db_e->cnt, db->mm[i]);
		if (k != kh_end(db_e->cnt) && kh_get_val(mm_hash, db_e->cnt, db->mm[i], -1) == 1) {
			single_ton++;
			p = kh_get_val(mm_hash, db_e->p, db->mm[i], -1);
			e = kh_get_val(mm_hash, db_e->h, db->mm[i], -1);
			if (p > MOLECULE_MARGIN && abs(g->edges[e].seq_len - p) > MOLECULE_MARGIN)
				continue;
			mm_hits_insert(hits, db->mm[i], e, p);
		}
	}
	return hits;
}

void mm_hits_print(struct mm_hits_t *hits, const char *file_path)
{
	FILE *f = fopen(file_path, "w");
	khiter_t k;
	uint32_t even, key;
	fprintf(f, "edge,Colour,hits\n");
	for (k = kh_begin(hits->edges); k != kh_end(hits->edges); ++k) {
		if (kh_exist(hits->edges, k)) {
			key = kh_key(hits->edges, k);
			even = key % 2 ? key - 1 : key;
			fprintf(f, "%d_%d,red,%d\n", even, even + 1, kh_value(hits->edges, k));
		}
	}
	fclose(f);
}

/*void mm_align(struct read_t r1, struct read_t r2, uint64_t bx, struct minimizer_bundle_t *bundle)
{
	int pos1 = 0, pos2 = 0;
	int n_reads = 0;
	struct mm_db_t *db1, *db2;
	struct mm_hits_t *hits1, *hits2, *hits;
	struct mini_hash_t *bx_table1 = (struct mini_hash_t *)()
	struct asm_graph_t *g = bundle->g;
	struct mm_db_edge_t *mm_edges_db = bundle->mm_edges;
	struct mini_hash_t *rp_table = *(bundle->rp_table);

	kh_mm_pw_t *kh_pw = kh_init(mm_pw);
	khiter_t ki, k;
	hits = mm_hits_init();

	n_reads++;
	db1 = mm_index_char_str(r1.seq, MINIMIZERS_KMER, MINIMIZERS_WINDOW, r1.len);
	db2 = mm_index_char_str(.seq, MINIMIZERS_KMER, MINIMIZERS_WINDOW, r2.len);
	khiter_t k1, k2;
	hits1 = mm_hits_init();
	hits2 = mm_hits_init();
	mm_hits_cmp(db1, mm_edges_db, hits1, g);
	mm_hits_cmp(db2, mm_edges_db, hits2, g);

	uint64_t max = 0, er1 = UINT64_MAX, er2 = UINT64_MAX;
	if (hits1->n > 0) {
		for (k1 = kh_begin(hits1->edges); k1 != kh_end(hits1->edges); ++k1) {
			if (kh_exist(hits1->edges, k1)) {
				if (kh_val(hits1->edges, k1) > max) {
					max = kh_val(hits1->edges, k1);
					er1 = kh_key(hits1->edges, k1);
				}
			}
		}
	}

	if (max < RATIO_OF_CONFIDENT * hits1->n && hits1->n > MIN_NUMBER_SINGLETON)
		er1 = UINT64_MAX;
	max = 0;
	if (hits2->n > 0) {
		for (k2 = kh_begin(hits2->edges); k2 != kh_end(hits2->edges); ++k2) {
			if (kh_exist(hits2->edges, k2))
				if (kh_val(hits2->edges, k2) > max) {
					max = kh_val(hits2->edges, k2);
					er2 = kh_key(hits2->edges, k2);
				}
		}
	}
	if (max < RATIO_OF_CONFIDENT * hits2->n && hits2->n > MIN_NUMBER_SINGLETON)
		er2 = UINT64_MAX;

	if (er1 != UINT64_MAX)
		kh_set(mm_edges, hits->edges, er1, 1);
	if (er2 != UINT64_MAX)
		kh_set(mm_edges, hits->edges, er2, 1);
	if (er1 != er2  && er1 != UINT64_MAX && er2 != UINT64_MAX && g->edges[er1].rc_id != er2) {
		printf("%llu %llu\n", er1, er2);
	}
	mm_hits_destroy(hits1);
	mm_hits_destroy(hits2);
	mm_db_destroy(db1);
	mm_db_destroy(db2);
}

*//*
* @brief Worker for counting barcode
* @param data
* @return
*//*
static inline void *minimizer_iterator(void *data)
{
	struct minimizer_bundle_t *bundle = (struct minimizer_bundle_t *) data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2, readbc;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_trip_buffer();

	char *R1_buf, *R2_buf;
	int pos1, pos2, rc1, rc2, input_format;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf)
			break;
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		R1_buf = ext_buf->R1_buf;
		R2_buf = ext_buf->R2_buf;
		input_format = ext_buf->input_format;

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
			      get_read_from_fq(&read1, R1_buf, &pos1) :
			      get_read_from_fa(&read1, R1_buf, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
			      get_read_from_fq(&read2, R2_buf, &pos2) :
			      get_read_from_fa(&read2, R2_buf, &pos2);

			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				log_error("Wrong format file");

			*//* read_name + \t + BX:Z: + barcode + \t + QB:Z: + barcode_quality + \n *//*
			//TODO: this assumes bx only came from one type of library
			uint64_t barcode = get_barcode_biot(read1.info, &readbc);
			if (barcode != (uint64_t) -1) {
				// any main stuff goes here
				uint64_t *slot = mini_put(bundle->bx_table, barcode);
			} else {
				//read doesn't have barcode
			}
			if (rc1 == READ_END)
				break;
		}
	}
	return NULL;
}

void mm_hit_all_barcodes(struct opt_proc_t *opt)
{
	uint64_t i;
	struct mini_hash_t *bx_table, *rp_table;
	bx_table = count_bx_freq(opt); //count barcode freq
	for (i = 0; i < bx_table->size; ++i) {
		if (bx_table->key[i] != 0 && bx_table->h[i] < MAX_READS_TO_HITS && bx_table->h[i] > MIN_READS_TO_HITS) {
			struct mini_hash_t *h;
			init_mini_hash(&h, 0);
			bx_table->h[i] = h;
		} else {
			bx_table->h[i] = EMPTY_BX;
		}
	}

	struct asm_graph_t *g = calloc(1, sizeof(struct asm_graph_t));
	assert(opt->in_file != NULL);
	load_asm_graph(g, opt->in_file);
	struct mm_db_edge_t *mm_edges = mm_index_edges(g, MINIMIZERS_KMER,
	                                               MINIMIZERS_WINDOW);
	init_mini_hash(&rp_table, 16);

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	void *(*buffer_iterator)(void *) = minimizer_iterator;
	struct producer_bundle_t *producer_bundles = NULL;
	//TODO: implement for sorted library
	producer_bundles = init_fastq_pair(opt->n_threads, opt->n_files,
	                                   opt->files_1, opt->files_2);

	struct minimizer_bundle_t *worker_bundles; //use an arbitrary structure for worker bundle
	worker_bundles = malloc(opt->n_threads * sizeof(struct minimizer_bundle_t));

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].bx_table = &bx_table;
		worker_bundles[i].mm_edges = mm_edges;
		worker_bundles[i].rp_table = &rp_table;
		worker_bundles[i].g = g;
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_producer,
		               producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, buffer_iterator,
		               worker_bundles + i);

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	free_fastq_pair(producer_bundles, opt->n_files);
}*/
