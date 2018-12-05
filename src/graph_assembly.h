#ifndef __GRAPH_ASSEMBLY_H__
#define __GRAPH_ASSEMBLY_H__

#include "attribute.h"

typedef int64_t gint_t;

struct scrap_graph_t {
	// 1-based node id, 0-based stored
	gint_t n_v;

	gint_t *kmer_count;
	kmkey_t **chains;

	gint_t *fhead, *rhead;
	gint_t *fadj, *radj;
};


void assembly_process(struct opt_count_t *opt);

#endif /* __GRAPH_ASSEMBLY_H__ */
