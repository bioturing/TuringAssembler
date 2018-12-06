#ifndef __GRAPH_ASSEMBLY_H__
#define __GRAPH_ASSEMBLY_H__

#include "attribute.h"
#include "kmhash.h"

typedef int64_t gint_t;

struct scrap_graph_t {
	gint_t *kmer_chain_id;

	// 1-based node id, 0-based stored
	gint_t n_v;

	gint_t *kmer_count;
	kmkey_t *kmer_beg;
	kmkey_t *kmer_end;

	gint_t **fadj, **radj;
	uint32_t *bin_fdeg, *bin_rdeg;
};


void assembly_process(struct opt_count_t *opt);

#endif /* __GRAPH_ASSEMBLY_H__ */
