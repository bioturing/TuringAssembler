/*
 * =====================================================================================
 *
 *       Filename:  count_barcodes.c
 *
 *    Description:  Count the frequency of each barcode
 *
 *        Version:  1.0
 *        Created:  09/22/2019 10:55:15 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Phan Thai Duy Tan (ptdtan), tan@bioturing.com
 *   Organization:  BioTuring Dev Team
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <stdbool.h>

#include "fastq_producer.h"
#include "attribute.h"
#include "dqueue.h"
#include "get_buffer.h"
#include "utils.h"
#include "verbose.h"
#include "atomic.h"
#include "log.h"

#include "count_barcodes.h"


/*
 Credit for primes table: Aaron Krowne
 http://br.endernet.org/~akrowne/
 http://planetmath.org/encyclopedia/GoodHashTablePrimes.html
 */
static const uint64_t primes[] = { 53, 97, 193, 389, 769, 1543, 3079, 6151,
                                   12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739,
                                   6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189,
                                   805306457, 1610612741 };
#define N_PRIMES_NUMBER 26

#define MAX_LOAD_FACTOR 0.8
#define FATAL_LOAD_FACTOR 0.9

pthread_mutex_t h_table_mut;

struct h_bundle_t {
    uint64_t *slot_ptr;
    void *arg;
};

void destroy_worker_bundles(struct readsort_bundle1_t *bundles, int n)
{
	destroy_mini_hash(*bundles[0].h_table_ptr);
	free(bundles);
}

/**
 * @brief Init a hash table with a pre-define size in the prime number table
 * @param p_index index of the prime number table
 * @return
 */
void init_mini_hash(struct mini_hash_t **h_table, uint32_t p_index)
{
	struct mini_hash_t *table = calloc(1, sizeof(struct mini_hash_t));
	table->prime_index = p_index;
	uint32_t h_size = primes[table->prime_index];
	table->h = calloc(h_size, sizeof(uint64_t));
	table->key = calloc(h_size, sizeof(uint64_t));
	memset(table->key, 255, sizeof(uint64_t) * h_size);
	for (int i = 0; i < h_size; i++) {
		assert(table->key[i] == -1);
	}
	table->size = h_size;
	table->max_cnt = (uint64_t) (table->size * MAX_LOAD_FACTOR);
	table->count = 0;
	*h_table = table;
}

void destroy_mini_hash(struct mini_hash_t *h_table)
{
	/** for (int i = 0; i < h_table->size; ++i){
		if (h_table->h[i] != EMPTY_BX)
			free((struct mini_hash_t *) h_table->h[i]);
	}
	 **/
	free(h_table->h);
	free(h_table->key);
}

uint64_t barcode_hash_mini(char *s)
{
	uint64_t ret = 0;
	for (int i = 0; i < 18; ++i) {
		ret = ret * 5 + nt4_table[(int) s[i]];
	}
	return ret;
}

uint64_t get_barcode_biot(char *s, struct read_t *r)
{
	if (s == NULL) {
		r->seq = r->qual = NULL;
		return (uint64_t) -1;
	}
	int i, len = 0;
	uint64_t ret = 0;
	char *p;
	p = strstr(s, "BX:Z:");
	if (p == NULL) {
		r->seq = r->qual = NULL;
		return (uint64_t) -1;
	}
	r->seq = p += 5;
	for (i = 0; p[i] && !__is_sep(p[i]); ++i) {
		ret = ret * 5 + nt4_table[(int) p[i]];
		++len;
	}
	p = strstr(s, "QB:Z:");
	if (p == NULL)
		r->qual = NULL;
	else
		r->qual = p + 5;
	r->len = len;
	return ret;
}

/*
 * Thomas Wang 64 bit mix hash function
 */

uint64_t twang_mix64(uint64_t key) {
	key = (~key) + (key << 21); // key *= (1 << 21) - 1; key -= 1;
	key = key ^ (key >> 24);
	key = key + (key << 3) + (key << 8); // key *= 1 + (1 << 3) + (1 << 8)
	key = key ^ (key >> 14);
	key = key + (key << 2) + (key << 4); // key *= 1 + (1 << 2) + (1 << 4)
	key = key ^ (key >> 28);
	key = key + (key << 31); // key *= 1 + (1 << 31)
	return key;
}

/*
 * Inverse of twang_mix64
 *
 * Note that twang_unmix64 is significantly slower than twang_mix64.
 */

uint64_t twang_unmix64(uint64_t key) {
	// See the comments in jenkins_rev_unmix32 for an explanation as to how this
	// was generated
	key *= 4611686016279904257U;
	key ^= (key >> 28) ^ (key >> 56);
	key *= 14933078535860113213U;
	key ^= (key >> 14) ^ (key >> 28) ^ (key >> 42) ^ (key >> 56);
	key *= 15244667743933553977U;
	key ^= (key >> 24) ^ (key >> 48);
	key = (key + 1) * 9223367638806167551U;
	return key;
}

/**
 * @brief Expand the hash table by re-hash and doubling size
 * when the load factor reach MAX_LOAD_FACTOR (0.65 by default)
 * key must be hashed by MurMurHash64
 */
void mini_expand(struct mini_hash_t **new_table_ptr)
{
	uint32_t i;
	uint64_t *slot;
	uint64_t key, data, val;
	uint32_t log_size = 1<<16;
	struct mini_hash_t *new_table;
	struct mini_hash_t *h_table = *new_table_ptr;
	init_mini_hash(&new_table, h_table->prime_index + 1);
	assert(h_table->h != NULL);
	for (i = 0; i < h_table->size; ++i) {
		if (!((i + 1) % log_size)) {
			log_info("Copied %d elements", log_size);
			log_size <<= 1;
		}
		if (__mini_empty(h_table, i))
			continue;
		val = h_table->h[i];
		data = h_table->key[i];
		key = twang_mix64(data);
		slot = mini_put_by_key(new_table, data, key);
		*slot = val;
	}
	destroy_mini_hash(h_table);
	*new_table_ptr = new_table;
}

/**
 * @brief Try expanding the hash table by doubling in size
 */
inline void try_expanding(struct mini_hash_t **h_table)
{
	struct mini_hash_t *table = *h_table;
	if (table->count >= table->max_cnt) {
		if (table->prime_index < N_PRIMES_NUMBER - 1) {
			log_info("Doubling hash table...");
			mini_expand(h_table);
		} else {
			if ((double)table->count > FATAL_LOAD_FACTOR * (table->size)) {
				log_warn("Hash table size reached limit!");
			}
		}
	}
}

/**
 * Return slot of data
 * @param data  barcode encoded as an uint64_t number
 * @param key   hash(data)
 */
uint64_t *mini_put_by_key(struct mini_hash_t *h_table, uint64_t data, uint64_t key)
{
	uint64_t i;
	uint64_t mask = h_table->size;
	uint64_t slot = key % mask;
	uint64_t is_empty = atomic_bool_CAS64(h_table->key + slot, EMPTY_SLOT, data);
	int go_around = 0;
	if (is_empty) { // slot is empty -> fill in
		atomic_add_and_fetch64(&(h_table->count), 1);
	} else if (!atomic_bool_CAS64(h_table->key + slot, data, data)) { // slot is reserved
		//linear probing
		for (i = slot + 1; i < h_table->size && !atomic_bool_CAS64(h_table->key + i, data, data); ++i) {
			is_empty = atomic_bool_CAS64(h_table->key + i, EMPTY_SLOT, data);
			if (is_empty)
				break;
		}
		if (i == h_table->size) {
			go_around = 1;
			for (i = 0; i < slot && !atomic_bool_CAS64(h_table->key + i, data, data); ++i) {
				is_empty = atomic_bool_CAS64(h_table->key + i, EMPTY_SLOT, data);
				if (is_empty)
					break;
			}
		}
		log_info("%lld", h_table->key[0]);
		if (i == slot) {
			for (int i = 0; i < h_table->size; i+= 1000) {
				log_info("%lld", h_table->key[i]);
			}
			log_info("size %d count %d maxcount %d key0 %lld EMPTY %lld goaround %d",
				h_table->size, h_table->count, h_table->max_cnt, h_table->key[0], EMPTY_SLOT, go_around);
			log_error("i == slot == %d", i);
		}
		assert(!atomic_bool_CAS64(&i, slot, slot));
		if (is_empty) //room at probe is empty -> fill in
			atomic_add_and_fetch64(&(h_table->count), 1);
		slot = i;
	}
	return h_table->h + slot;
}

/**
 * @param data  barcode encoded as an uint64_t number
 * @param key   hash(data)
 */
uint64_t *mini_get_by_key(struct mini_hash_t *h_table, uint64_t data, uint64_t key)
{
	uint64_t i;
	uint64_t mask = h_table->size;
	uint64_t slot = key % mask;
	uint64_t is_empty = atomic_bool_CAS64(h_table->key + slot, EMPTY_SLOT, EMPTY_SLOT);
	if (is_empty) { // slot is empty -> fill in
		return (uint64_t *)EMPTY_BX;
	} else if (!atomic_bool_CAS64(h_table->key + slot, data, data)) { // slot is reserved
		//linear probing
		for (i = slot + 1; i < h_table->size && !atomic_bool_CAS64(h_table->key + i, data, data); ++i) {
			is_empty = atomic_bool_CAS64(h_table->key + i, EMPTY_SLOT, EMPTY_SLOT);
			if (is_empty)
				break;
		}
		if (i == h_table->size) {
			for (i = 0; i < slot && !atomic_bool_CAS64(h_table->key + i, data, data); ++i) {
				is_empty = atomic_bool_CAS64(h_table->key + i, EMPTY_SLOT, EMPTY_SLOT);
				if (is_empty)
					break;
			}
		}
		assert(!atomic_bool_CAS64(&i, slot, slot));
		if (is_empty) //room at probe is empty -> fill in
			return (uint64_t *)EMPTY_BX;
		slot = i;
	}
	return h_table->h + slot;
}

/**
 * @brief
 * @param data
 * @param key
 * @return count of data
 */

uint64_t *mini_get(struct mini_hash_t *h_table, uint64_t data)
{
	uint64_t key = twang_mix64(data);
	uint64_t *slot = mini_get_by_key(h_table, data, key);
	return slot;
}

/**
 * @brief Return slot of the data
 * @param data byte array of data
 * @param len length in byte of data
 */
uint64_t *mini_put(struct mini_hash_t **h_table, uint64_t data)
{
	struct mini_hash_t *table = *h_table;
	uint64_t key = twang_mix64(data);
//	if (table->count > table->max_cnt) {
	if (atomic_bool_CAS64(&table->count, table->max_cnt, table->max_cnt)){
		pthread_mutex_lock(&h_table_mut);
		try_expanding(h_table);
		pthread_mutex_unlock(&h_table_mut);
	}
	return mini_put_by_key(*h_table, data, key);
}

/**
 * @brief Print the barcode frequencies to a file
 * @param bx_size
 * @param out_dir
 */
void mini_print(struct mini_hash_t *h_table, size_t bx_size, char *out_dir)
{
	char count_file[1024];
	sprintf(count_file, "%s/barcode_frequencies.txt", out_dir);
	FILE *fp = fopen(count_file, "w");
	uint32_t i, j, c;
	char bx[bx_size + 1];
	char nt5[5] = "ACGTN";

	memset(bx, '\0', bx_size + 1);
	for (i = 0; i < h_table->size; ++i) {
		if (h_table->h[i] != 0) {
			j = bx_size;
			uint64_t ret = h_table->key[i];
			assert(ret != EMPTY_SLOT);
			while (j) {
				c = ret % 5;
				bx[--j] = nt5[c];
				ret = (ret - c) / 5;
			}
			fprintf(fp, "%s\t%lu\n", bx, h_table->h[i]);
		}
	}
	fclose(fp);
}

static inline void *biot_buffer_iterator_simple(void *data);

/**
 * @brief main function that count barcode frequencies
 * @param opt
 * @param r_path
 */
struct mini_hash_t *count_bx_freq(struct opt_proc_t *opt)
{
	struct mini_hash_t *h_table;
	init_mini_hash(&h_table, INIT_PRIME_INDEX);

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	void *(*buffer_iterator)(void *) = biot_buffer_iterator_simple;
	struct producer_bundle_t *producer_bundles = NULL;
	producer_bundles = init_fastq_pair(opt->n_threads, opt->n_files,
	                                   opt->files_1, opt->files_2);
	buffer_iterator = biot_buffer_iterator_simple;
	struct readsort_bundle1_t *worker_bundles; //use an arbitrary structure for worker bundle
	worker_bundles = malloc(opt->n_threads * sizeof(struct readsort_bundle1_t));
	int i;
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].h_table_ptr = &h_table;
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

	free(producer_threads);
	free(worker_threads);
	free(worker_bundles);

	free_fastq_pair(producer_bundles, opt->n_files);
	mini_print(h_table, 18, opt->out_dir); // 18 is the barcode length from TELL-Seq technology
	return h_table;
}

/**
 * @brief Worker for counting barcode
 * @param data
 * @return
 */
static inline void *biot_buffer_iterator_simple(void *data)
{
	struct readsort_bundle1_t *bundle = (struct readsort_bundle1_t *) data;
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
				log_error("\nWrong format file\n");

			/* read_name + \t + BX:Z: + barcode + \t + QB:Z: + barcode_quality + \n */
			uint64_t barcode = get_barcode_biot(read1.info, &readbc);
			if (barcode != (uint64_t) -1) {
				// any main stuff goes here
				uint64_t *slot = mini_put(bundle->h_table_ptr, barcode);
				atomic_add_and_fetch64(slot, 1);
			} else {
				//read doesn't have barcode
			}
			if (rc1 == READ_END)
				break;
		}
	}
	return NULL;
}

uint64_t get_min_code(struct asm_graph_t *g, int u, int v)
{
	if (g->edges[u].rc_id < u)
		u = g->edges[u].rc_id;
	if (g->edges[v].rc_id < v)
		v = g->edges[v].rc_id;
	if (u > v) {
		int tmp = u;
		u = v;
		v = tmp;
	}
	return ((uint64_t)u << 32) + v;
}

/**
 * @brief: Create the hash table of share-barcode count of every edge pair u-v
 * @param: edge hits of all barcodes. See @function mm_hit_all_barcodes
 * @return: khash, key = (uint64_t)(e1|e2), e1 < e2, value: share barcode count
 */
khash_t(long_int) *count_edge_link_shared_bc(struct asm_graph_t *g,
		struct mini_hash_t *bc)
{
	khash_t(long_int) *res = kh_init(long_int);
	for (int i = 0; i < bc->size; ++i){
		if (bc->h[i] == EMPTY_BX)
			continue;
		struct mini_hash_t *wrapper = (struct mini_hash_t *)bc->h[i]; // ptr to mini_hash_t instead of count
		uint64_t *edges = calloc(wrapper->size, sizeof(uint64_t));
		int n_e = 0;
		for (int j = 0; j < wrapper->size; ++j){
			if (wrapper->h[j] == 0)
				continue;
			edges[n_e++] = wrapper->key[j];
		}

		for (int j = 0; j < n_e; ++j){
			for (int k = j + 1; k < n_e; ++k){
				int e1 = edges[j];
				int e2 = edges[k];
				if (e1 >= g->n_e || e2 >= g->n_e) {
					log_error("Wrong edges number of one barcode hit, e1 %d, e2 %d", e1, e2);
				}
				uint64_t code = get_min_code(g, e1, e2);

				int cur_count = kh_long_int_try_get(res, code, 0);
				kh_long_int_set(res, code, cur_count + 1);
			}
		}
		free(edges);
	}

	return res;
}

/**
 * @brief Query share barcode of a pair e1-e2
 * @param all_count
 * @param e1
 * @param e2
 * @return share barcode count
 */
int get_edge_link_bc_count(khash_t(long_int) *all_count, int e1, int e2)
{
	uint64_t code = (e1 <= e2) ? GET_CODE(e1, e2) : GET_CODE(e2, e1);
	return kh_long_int_try_get(all_count, code, 0);
}
