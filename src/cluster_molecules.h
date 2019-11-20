#ifndef __CLUSTER_MOLECULES__
#define __CLUSTER_MOLECULES__
#include "assembly_graph.h"
#include "complex_resolve.h"
#include "minimizers/minimizers.h"
#include "sort_read.h"
#include "get_buffer.h"
KHASH_MAP_INIT_INT64(long_int, int);

struct dijkstra_node_t{
	int vertex;
	int len;
};

int get_shortest_path(struct asm_graph_t *g, int source, int target);
void count_edge_links_bc(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct read_path_t *read_sorted_path, khash_t(bcpos) *bx_pos_dict,
		struct mm_db_edge_t *mm_edges, char **bc_list, int n_bc);
void get_sub_graph(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct mm_hits_t *hits, khash_t(long_int) *pair_count);
#endif
