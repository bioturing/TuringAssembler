#ifndef __contig_graph__
#define __contig_graph__

struct contig_edge {
	uint32_t src, des;
	uint32_t rv_src, rv_des;
	float score0;
};

#endif

