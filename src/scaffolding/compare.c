#include "assembly_graph.h"
#include "scaffolding/edge.h"

int ascending_unint32(const void *e0,const void *e1)
{
	uint32_t a = *(uint32_t*)(e0);
	uint32_t b = *(uint32_t*)(e1);
	return (a > b) - (a < b);
}

int decending_uint32(const void *e0,const void *e1)
{
	uint32_t a = *(uint32_t*)(e0);
	uint32_t b = *(uint32_t*)(e1);
	return (a < b) - (a > b);
}

int decending_edge_score(const void *d0, const void *d1)
{
	float a =  ((struct contig_edge *) d0)->score0;
	float b =  ((struct contig_edge *) d1)->score0;
	struct contig_edge *e1 =  (struct contig_edge *) d1;
	return (a < b) - (a > b);
}



int decending_candidate_edge(const void *d0, const void *d1)
{
	float a =  ((struct candidate_edge *) d0)->score;
	float b =  ((struct candidate_edge *) d1)->score;
	return (a < b) - (a > b);
}

