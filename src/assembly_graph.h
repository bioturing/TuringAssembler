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
	struct barcode_hash_t *bucks;
	pthread_mutex_t lock;
};

struct asm_graph_t {
	int ksize;
	gint_t n_v, n_e;

	struct asm_node_t *nodes;
	struct asm_edge_t *edges;
};

void assembly_process(struct opt_count_t *opt);

/* deprecated */
void k31_process(struct opt_count_t *opt);

/* deprecated */
void k63_process(struct opt_count_t *opt);

void test_asm_graph(struct asm_graph_t *g);

/* deprecated */
void build_asm_graph_from_k31(struct opt_count_t *opt, int ksize,
			struct k31hash_t *kmer_hash, struct asm_graph_t *ret_g);

/* deprecated */
void build_asm_graph_from_k63(struct opt_count_t *opt, int ksize,
			struct k63hash_t *kmer_hash, struct asm_graph_t *ret_g);

/* should not put here */
int asm_is_edge_rc(uint32_t *seq1, gint_t l1, uint32_t *seq2, gint_t l2);

void save_asm_graph(struct asm_graph_t *g, const char *path);

void load_asm_graph(struct asm_graph_t *g, const char *path);

void build0_1_process(struct opt_build_t *opt);

void build1_2_process(struct opt_build_t *opt);

void build2_3_process(struct opt_build_t *opt);

void build0_process(struct opt_count_t *opt);

void graph_convert_process(struct opt_build_t *opt);

void construct_barcode_map(struct asm_graph_t *g, struct opt_build_t *opt);

#endif  /* __ASSEMBLY_GRAPH_H__ */
