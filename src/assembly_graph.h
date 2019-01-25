#ifndef __ASSEMBLY_GRAPH_H__
#define __ASSEMBLY_GRAPH_H__

#include <stdint.h>

#include "kmhash.h"

typedef int64_t gint_t;

struct asm_node_t {
	/* node represents a 31-mer */
	kmkey_t seq;
	gint_t forward_adj[4];
	gint_t reverse_adj[4];
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

void test_graph_build(struct asm_graph_t *graph);
#endif  /* __ASSEMBLY_GRAPH_H__ */
