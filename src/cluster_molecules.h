#ifndef __CLUSTER_MOLECULES__
#define __CLUSTER_MOLECULES__
#include "assembly_graph.h"
#include "complex_resolve.h"
#include "minimizers/minimizers.h"
#include "sort_read.h"
#include "get_buffer.h"
#include "complex_resolve.h"

KHASH_MAP_INIT_INT64(long_int, int);
struct dijkstra_node_t{
	int vertex;
	int len;
};

struct bc_edges_path_t{
	char bc[19];
	int n_e;
	int *edges;
};

struct barcode_list_t{
	int n_bc;
	char **bc_list;
};

struct simple_node_t{
	int deg;
	int *adj;
};
KHASH_MAP_INIT_INT(int_node, struct simple_node_t *);

struct simple_graph_t{
	khash_t(set_int) *is_loop;
	khash_t(int_node) *nodes;
};

void init_simple_graph(struct simple_graph_t *sg);

int get_shortest_path(struct asm_graph_t *g, int source, int target);
void count_edge_links_bc(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct read_path_t *read_sorted_path, khash_t(bcpos) *bx_pos_dict,
		struct mm_db_edge_t *mm_edges, char **bc_list, int n_bc);
void get_sub_graph(struct asm_graph_t *g, struct mm_hits_t *hits,
		khash_t(long_int) *pair_count);
void print_barcode_graph(struct opt_proc_t *opt);
void get_barcode_edges_path(struct opt_proc_t *opt);
void get_barcode_list(char *bc_count_path, struct barcode_list_t *blist);
void get_all_pair_edge_count(char *file_path, khash_t(long_int) *pair_count);
void add_simple_node(struct simple_graph_t *sg, int u);
void add_simple_edge(struct simple_graph_t *sg, int u, int v);
void build_simple_graph(khash_t(long_int) *one_bc, khash_t(long_int) *all_bc,
		struct simple_graph_t *sg);
void simple_graph_destroy(struct simple_graph_t *sg);
void check_loop_dfs(struct simple_graph_t *sg, int u, khash_t(set_int) *visited,
		khash_t(set_int) *in_dfs);
void find_DAG(struct simple_graph_t *sg, struct asm_graph_t *g);
#endif
