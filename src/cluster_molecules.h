#ifndef __CLUSTER_MOLECULES__
#define __CLUSTER_MOLECULES__
#include "assembly_graph.h"
#include "minimizers/minimizers.h"
#include "sort_read.h"
#include "get_buffer.h"
#include "simple_queue.h"
#include "khash_operations.h"
KHASH_MAP_INIT_INT64(long_int, int);
struct shortest_path_info_t{
	int len;
	int trace;
};

KHASH_MAP_INIT_INT64(long_spath, struct shortest_path_info_t *);
KHASH_SET_INIT_INT64(set_long);

KHASH_MAP_OPERATIONS(long_spath, uint64_t, struct shortest_path_info_t *);
KHASH_SET_OPERATIONS(set_int, int);
KHASH_MAP_OPERATIONS(long_int, uint64_t, int);
KHASH_MAP_OPERATIONS(int_int, int, int);
KHASH_SET_OPERATIONS(set_long, uint64_t);

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
	struct asm_graph_t *g;
	khash_t(set_int) *is_loop;
	khash_t(int_int) *path_len;
	khash_t(int_int) *next;
	khash_t(int_node) *nodes;
};

void init_simple_graph(struct asm_graph_t *g, struct simple_graph_t *sg);

void get_all_shortest_paths_dp(struct asm_graph_t *g, khash_t(long_spath) *spath_info);

int get_pair_distance(int v, int u, khash_t(long_spath) *spath_info);

void get_edge_links_by_distance(struct asm_graph_t *g, int *edges, int n_e,
		khash_t(long_spath) *spath_info, khash_t(long_int) *is_connected,
		khash_t(long_int) *count_link);

int check_connected(struct asm_graph_t *g, int v, int u,
		khash_t(long_spath) *spath_info);

void count_edge_links_bc(struct opt_proc_t *opt);

void print_barcode_graph(struct opt_proc_t *opt);

void get_barcode_edges_path(struct opt_proc_t *opt);

void get_simple_components(struct opt_proc_t *opt);

void get_barcode_list(char *bc_count_path, struct barcode_list_t *blist);

void barcode_list_destroy(struct barcode_list_t *blist);

void get_all_pair_edge_count(char *file_path, khash_t(long_int) *pair_count);

void add_simple_node(struct simple_graph_t *sg, int u);

void add_simple_edge(struct simple_graph_t *sg, int u, int v);

void build_simple_graph(struct mm_hits_t *hits, khash_t(long_int) *all_bc,
		struct simple_graph_t *sg);

void build_simple_bigraph(struct mm_hits_t *hits, khash_t(long_int) *all_bc,
		struct simple_graph_t *sg);
void simple_graph_destroy(struct simple_graph_t *sg);
void check_loop_dfs(struct simple_graph_t *sg, int u, khash_t(set_int) *visited,
		khash_t(set_int) *in_dfs);
void find_DAG(struct simple_graph_t *sg, struct asm_graph_t *g);
void get_longest_path_dfs(struct simple_graph_t *sg, int u,
		khash_t(set_int) *done_dfs);
void get_longest_path(struct simple_graph_t *sg);
int cmp_dijkstra(void *node1, void *node2);
int is_repeat(struct asm_graph_t *g, int e);

void hits_to_edges(struct asm_graph_t *g, struct mm_hits_t *hits, int **edges,
		int *n_e);

void bfs_nearby(struct asm_graph_t *g, int source, int radius, int **edges, int *n_e);
void print_graph_component(struct simple_graph_t *sg, char *bc, FILE *f);
#endif
