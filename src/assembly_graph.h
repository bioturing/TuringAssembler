#ifndef __ASSEMBLY_GRAPH_H__
#define __ASSEMBLY_GRAPH_H__

#include <stdint.h>

#include "attribute.h"
#include "k31hash.h"
#include "k63hash.h"

struct asm_node_t {
	gint_t rc_id;
	gint_t deg;
	gint_t *adj;
};

struct asm_edge_t {
	uint64_t count;
	uint32_t *seq;
	gint_t seq_len;
	gint_t source;
	gint_t target;
	gint_t rc_id;
};

struct asm_graph_t {
	int ksize;
	gint_t n_v, n_e;

	struct asm_node_t *nodes;
	struct asm_edge_t *edges;
};

void k31_process(struct opt_count_t *opt);

void k63_process(struct opt_count_t *opt);

void test_asm_graph(struct asm_graph_t *g);

void build_asm_graph_from_k31(struct opt_count_t *opt, int ksize,
			struct k31hash_t *kmer_hash, struct asm_graph_t *ret_g);

void build_asm_graph_from_k63(struct opt_count_t *opt, int ksize,
			struct k63hash_t *kmer_hash, struct asm_graph_t *ret_g);

int asm_is_edge_rc(uint32_t *seq1, gint_t l1, uint32_t *seq2, gint_t l2);

#endif  /* __ASSEMBLY_GRAPH_H__ */
