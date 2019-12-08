//
// Created by BioTuring on 2019-12-08.
//

#include "../assembly_graph.h"
#include "../minimizers/minimizers.h"
#include "../minimizers/count_barcodes.h"

#include "../log.h"

#define KMER_SIZE_COVERAGE 31


/**
 * Get the i kmer of the encoded DNA sequence s (from minimizers/minimizers.c)
 * @param s DNA sequence s
 * @param i pos
 * @param k k size
 * @return encoded kmer
 */
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

//Only for consecutive edge
int index_bin_edge(struct mini_hash_t **table, struct asm_edge_t e)
{
	uint32_t i;
	uint64_t c, cnt = 0;
	uint64_t *slot;
	int pad = (32 - KMER_SIZE_COVERAGE - 1)*2;
	uint64_t km = get_km_i_bin(e.seq, 0, KMER_SIZE_COVERAGE);
	for (i = 0 ; i < e.seq_len - KMER_SIZE_COVERAGE + 1; ++i) {
		c = (uint64_t)__binseq_get(e.seq, i + KMER_SIZE_COVERAGE - 1);
		km |= ((uint64_t) c << (pad + 2));
		slot = mini_put(table, km);
		if (*slot == 0)
			cnt++;
		km <<= 2;
	}
	return cnt;
}

struct mini_hash_t * construct_edges_hash(struct asm_graph_t *g)
{
	int64_t i, kmer_cnt = 0;
	struct mini_hash_t *table;
	init_mini_hash(&table, INIT_PRIME_INDEX);
	for (i = 0; i < g->n_e; ++i) {
		if (g->edges[i].seq_len < (KMER_SIZE_COVERAGE + 1))
			continue;
		kmer_cnt += index_bin_edge(&table, g->edges[i]);
	}
	log_info("Total number of kmer %d", kmer_cnt);
	return table;
}
