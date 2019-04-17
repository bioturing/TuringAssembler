#ifndef __COMPARE_H__
#define __COMPARE_H__
#include "assembly_graph.h"
#include "contig_graph.h"
int less_contig_edge(const void *e0,const void *e1);
int less_uint32(const void *e0,const void *e1);
uint32_t equal_contig_edge(struct contig_edge *e0,struct contig_edge *e1);
uint32_t better_contig_edge(struct contig_edge *e0, struct contig_edge *e1);
int greater_uint32(const void *e0,const void *e1);

#endif
