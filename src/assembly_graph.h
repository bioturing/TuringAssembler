#ifndef __ASSEMBLY_GRAPH_H__
#define __ASSEMBLY_GRAPH_H__

#include <stdint.h>

#include "kmhash.h"

typedef int64_t gint_t;

struct asm_node0_t {
	kmkey_t seq;
	gint_t forward_adj[4];
	gint_t reverse_adj[4];
};

struct asm_edge0_t {
	uint64_t count;
	uint32_t *seq;
	gint_t seq_len;
	gint_t source;
	gint_t target;
	gint_t rc_id;
};

struct asm_graph0_t {
	int ksize;
	gint_t n_v, n_e;

	struct asm_node0_t *nodes;
	struct asm_edge0_t *edges;
};

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
	// struct kmhash2_t bc_set;
};

struct asm_graph_t {
	int ksize;
	gint_t n_v, n_e;

	struct asm_node_t *nodes;
	struct asm_edge_t  *edges;
};

#endif  /* __ASSEMBLY_GRAPH_H__ */
