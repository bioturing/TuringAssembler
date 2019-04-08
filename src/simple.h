#ifndef __SIMPLE_FOREST_H__
#define __SIMPLE_FOREST_H__

#include <stdint.h>

#include "attribute.h"
#include "assembly_graph.h"
#include "k31hash.h"
#include "k63hash.h"
#include "khash.h"

/* @abstract : compute the genome coverage using only contig > 7k
 * @param g0: level3 assembly graph struct
 * @return : genome coverage
 */
uint32_t get_genome_cov(struct asm_graph_t *g0);

/* @abstract: whether an edge is the bridge of four sequence
 */
int is_bridge(struct asm_edge_t *e, struct asm_node_t *v, gint_t i);


/* @abstract: whether an edge is an element of a loop (the one with 
 * double coverage)
 */
int is_trivial_loop(struct asm_edge_t *e, struct asm_node_t *v, gint_t i);


/* @abstract: whether an edge is the bridge of four sequence
 */
int is_simple_tandem(struct asm_edge_t *e, struct asm_node_t *v, gint_t e_i,
			gint_t *k, uint32_t *comp_sz);

#endif //__SIMPLE_FOREST_H__


