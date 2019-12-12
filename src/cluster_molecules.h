#ifndef __CLUSTER_MOLECULES__
#define __CLUSTER_MOLECULES__
#include "assembly_graph.h"
#define MIN_EDGE_LEN 500
#define MAX_RADIUS 4000
#define MAX_PATH_LEN 30
#define MIN_BC_READ_COUNT 10
#define MAX_BC_READ_COUNT 88
#define MIN_BARCODE_EDGE_COUNT 100
#define MIN_COVERAGE_TO_BE_IGNORE 0.25
#define COVERAGE_RATIO_TO_BE_REPEAT 1.75
#include "minimizers/minimizers.h"
#include "sort_read.h"
#include "get_buffer.h"
#include "simple_queue.h"
#include "khash_operations.h"
#define GET_CODE(a, b) ((((uint64_t) (a)) << 32) | (b))
KHASH_MAP_INIT_INT64(long_int, int);
KHASH_MAP_OPERATIONS(long_int, uint64_t, int);

KHASH_SET_INIT_INT64(set_long);
KHASH_SET_OPERATIONS(set_long, uint64_t);


struct barcode_list_t{
	int n_bc;
	char **bc_list;
	int *read_count;
};

struct simple_node_t{
	int deg;
	int *adj;
	int rv_deg;
	int *rv_adj;
};
KHASH_MAP_INIT_INT(int_node, struct simple_node_t *);
KHASH_MAP_OPERATIONS(int_node, int, struct simple_node_t *);

struct simple_graph_t{
	struct asm_graph_t *g;
	khash_t(set_int) *is_loop;
	khash_t(set_int) *is_complex;
	khash_t(int_int) *path_len;
	khash_t(int_int) *next;
	khash_t(int_node) *nodes;
};

struct simple_path_t{
	int *edges;
	int n_e;
};

struct paths_bundle_t{
	struct simple_path_t *paths;
	int n_paths;
};

struct shortest_path_info_t{
	int sum_seq;
	int n_e;
	int *path;
};

struct len_info_t{
	int v;
	int sum_seq;
	int len;
};

KHASH_MAP_INIT_INT64(long_spath, struct shortest_path_info_t *);
KHASH_MAP_OPERATIONS(long_spath, uint64_t, struct shortest_path_info_t *);

void init_simple_graph(struct asm_graph_t *g, struct simple_graph_t *sg);

//void get_all_shortest_paths_dp(struct asm_graph_t *g, khash_t(long_spath) *spath_info);
/*int extract_shortest_path(struct asm_graph_t *g, khash_t(long_spath) *spath,
		int v, int u, int **path, int *n_v);*/

void trace_shortest_path(int s, int t, khash_t(long_int) *best_P,
		khash_t(int_int) *best_D, int **path, int *n_e);

struct shortest_path_info_t *get_shortest_path(struct asm_graph_t *g, int s,
		int t, khash_t(long_spath) *stored);

int get_pair_distance(int v, int u, khash_t(long_spath) *spath_info);

void get_edge_links_by_distance(struct asm_graph_t *g, int *edges, int n_e,
		khash_t(long_spath) *spath_info, khash_t(long_int) *is_connected,
		khash_t(long_int) *count_link);

int check_connected(struct asm_graph_t *g, int v, int u,
		khash_t(long_spath) *spath_info);

//void count_edge_links_bc(struct opt_proc_t *opt);

void print_barcode_graph(struct opt_proc_t *opt);

void get_barcode_edges_path(struct opt_proc_t *opt);

void get_barcode_list(char *bc_count_path, struct barcode_list_t *blist);

void barcode_list_destroy(struct barcode_list_t *blist);

void add_simple_node(struct simple_graph_t *sg, int u);

void add_simple_edge(struct simple_graph_t *sg, int u, int v);

void build_simple_graph(int *edges, int n_e, khash_t(long_int) *all_bc,
		struct simple_graph_t *sg);

void build_graph_from_edges_list(int *edges, int n_e, struct asm_graph_t *g,
		struct simple_graph_t *sg);

void simple_graph_destroy(struct simple_graph_t *sg);
void check_loop_dfs(struct simple_graph_t *sg, int u, khash_t(set_int) *visited,
		khash_t(set_int) *in_dfs);
void find_DAG(struct simple_graph_t *sg);
void filter_complex_regions(struct simple_graph_t *bi_sg);
void get_longest_path_dfs(struct simple_graph_t *sg, int u,
		khash_t(set_int) *done_dfs);
void get_longest_path(struct simple_graph_t *sg);

void create_barcode_molecules(struct opt_proc_t *opt, int *edges, int n_e,
		struct asm_graph_t *g);

int is_repeat(struct asm_graph_t *g, int e);

void hits_to_edges(struct asm_graph_t *g, struct mm_hits_t *hits, int **edges,
		int *n_e);

void bfs_nearby(struct asm_graph_t *g, int source, int radius, int **edges, int *n_e);
void print_graph_component(struct simple_graph_t *sg, char *bc, FILE *f);

void load_pair_edge_count(char *path, khash_t(long_int) *h_all);
void print_simple_graph(struct simple_graph_t *sg, int *edges, int n_e, FILE *f);
void fill_gap(struct asm_graph_t *g, int v, int u, khash_t(long_spath) *spath,
		struct simple_graph_t *sg, char **seq);
void get_all_pair_edges(struct asm_graph_t *g, khash_t(long_int) *pair_edges);
void get_all_longest_paths(int *edges, int n_e, struct asm_graph_t *g,
		struct paths_bundle_t *path_bundle);
int check_adj_edges(struct asm_graph_t *g, int v, int u);
#endif
