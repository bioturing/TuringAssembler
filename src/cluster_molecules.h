#ifndef __CLUSTER_MOLECULES__
#define __CLUSTER_MOLECULES__
#include "assembly_graph.h"
#include "minimizers/minimizers.h"
#include "sort_read.h"
#include "get_buffer.h"
#include "simple_queue.h"
#include "khash_operations.h" 
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
	int *read_count;
};

struct simple_node_t{
	int deg;
	int *adj;
};
KHASH_MAP_INIT_INT(int_node, struct simple_node_t *);

struct simple_graph_t{
	khash_t(set_int) *is_loop;
	khash_t(int_int) *path_len;
	khash_t(int_int) *next;
	khash_t(int_node) *nodes;
};

void init_simple_graph(struct simple_graph_t *sg);

int get_shortest_path(struct asm_graph_t *g, int source, int target);
void dijkstra(struct asm_graph_t *g, int source, khash_t(int_int) *distance);
void get_all_shortest_paths(struct asm_graph_t *g, khash_t(long_int) *distance);
int get_pair_distance(int v, int u, khash_t(long_int) *distance);
void get_edge_links_by_distance(struct asm_graph_t *g, int *edges, int n_e,
		khash_t(long_int) *distance, khash_t(long_int) *count_link);
void count_edge_links_bc(struct opt_proc_t *opt);
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
void get_longest_path_dfs(struct simple_graph_t *sg, int u,
		khash_t(set_int) *done_dfs);
void get_longest_path(struct simple_graph_t *sg);
int cmp_dijkstra(void *node1, void *node2);
#endif
