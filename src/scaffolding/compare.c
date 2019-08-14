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
	float a =  ((struct scaffold_edge *) d0)->score.bc_score;
	float b =  ((struct scaffold_edge *) d1)->score.bc_score;
	struct scaffold_edge *e1 =  (struct scaffold_edge *) d1;
	return (a < b) - (a > b);
}

int ascending_edge(const void *e0, const void *e1)
{
	struct scaffold_edge *a = (struct scaffold_edge *) e0;
	struct scaffold_edge *b = (struct scaffold_edge *) e1;
	int t0 = (a->src > b->src) - (a->src < b->src);
	int t1 = (a->des > b->des) - (a->des < b->des);
	if (t0 != 0)
		return t0;
	return t1;
}

int decending_candidate_edge(const void *d0, const void *d1)
{
	float a =  ((struct candidate_edge *) d0)->score.bc_score;
	float b =  ((struct candidate_edge *) d1)->score.bc_score;
	return (a < b) - (a > b);
}

int ascending_index_edge(const void *e0, const void *e1)
{
	struct candidate_edge *a = (struct candidate_edge *) e0;
	struct candidate_edge *b = (struct candidate_edge *) e1;
	int t0 = (a->src > b->src) - (a->src < b->src);
	int t1 = (a->des > b->des) - (a->des < b->des);
	if (t0 != 0)
		return t0;
	return t1;
}

int ascending_scaffold_edge_index(const void *e0, const void *e1)
{
	struct scaffold_edge *a = (struct scaffold_edge *) e0;
	struct scaffold_edge *b = (struct scaffold_edge *) e1;
	int t0 = (a->src > b->src) - (a->src < b->src);
	int t1 = (a->des > b->des) - (a->des < b->des);
	if (t0 != 0)
		return t0;
	return t1;
}

