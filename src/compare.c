#include "assembly_graph.h"
#include "contig_graph.h"
int less_contig_edge(const void *e0,const void *e1)
{
	struct contig_edge * v0 = (struct contig_edge *) e0;
	struct contig_edge * v1 = (struct contig_edge *) e1;
	return (v0->src > v1->src || (v0->src == v1->src && v0->des > v1->des));
}

int less_uint32(const void *e0,const void *e1)
{
	return *(uint32_t*)(e0) < *(uint32_t*)(e1);
}

int greater_uint32(const void *e0,const void *e1)
{
	return *(uint32_t*)(e0) > *(uint32_t*)(e1);
}

uint32_t equal_contig_edge(struct contig_edge *e0,struct contig_edge *e1)
{
	return (e0->src == e1->src && e0->des == e1->des);
}

uint32_t better_contig_edge(struct contig_edge *e0, struct contig_edge *e1)
{
	return (e0->score0 > e1->score0);
}

