#ifndef __SIMPLE_FOREST_H__
#define __SIMPLE_FOREST_H__

#include <stdint.h>

#include "attribute.h"
#include "assembly_graph.h"
#include "k31hash.h"
#include "k63hash.h"

/* @abstract : compute the genome coverage using only contig > 7k
 * @param g0: level3 assembly graph struct
 * @return : genome coverage
 */
uint32_t get_genome_cov(struct asm_graph_t *g0);

void find_forest(struct asm_graph_t *g0);

#endif //__SIMPLE_FOREST_H__


