#include "assembly_graph.h"

int less_uint32(const void *e0,const void *e1)
{
	return *(uint32_t*)(e0) < *(uint32_t*)(e1);
}

int greater_uint32(const void *e0,const void *e1)
{
	return *(uint32_t*)(e0) > *(uint32_t*)(e1);
}

