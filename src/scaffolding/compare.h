#ifndef SCAFFOLDING_COMPARE_H
#define SCAFFOLDING_COMPARE_H
#include "assembly_graph.h"

int ascending_unint32(const void *e0,const void *e1);
int decending_uint32(const void *e0,const void *e1);
int decending_edge_score(const void *d0, const void *d1);
int decending_candidate_edge(const void *d0, const void *d1);
#endif
