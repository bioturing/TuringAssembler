#ifndef __GRAPH_SEARCH__
#define __GRAPH_SEARCH__
#define MIN_DEPTH_RATIO 0.2
#define MIN_PATH_LENGTH 10000
#define MIN_RELATIVE_COV_RATIO 0.2
#define SIMPLE_PATH 0
#define COMPLEX_PATH 1
#define MAX_PATH_COUNT 100
#include "khash.h"
#include "assembly_graph.h"
#include <stdio.h>
#define PATH_NOT_FOUND -1
KHASH_MAP_INIT_INT64(gint_int, int);
struct graph_info_t{
	struct asm_graph_t *g;
	int start_edge;
	int end_edge;
	int *is_edge_trash;
	int *is_edge_vst;
	khash_t(gint_int) *is_link_trash;
	khash_t(gint_int) *is_link_vst;
	khash_t(gint_int) *link_max_vst;
};

struct path_info_t{
	int n_paths;
	int m_paths;
	int *path_lens;
	int **paths;
};

void graph_info_init(struct asm_graph_t *lg, struct graph_info_t *ginfo,
			int start_edge, int end_edge);
void graph_info_init_max_vst(struct graph_info_t *ginfo);
int check_edge_trash(struct graph_info_t *ginfo, int e);
int check_link_trash(struct graph_info_t *ginfo, int e1, int e2);
int check_edge_visted(struct graph_info_t *ginfo, int e);
void graph_info_destroy(struct graph_info_t *ginfo);
gint_t get_edge_code(gint_t u, gint_t v);
void mark_edge_trash(struct graph_info_t *ginfo, int e);
void mark_edge_visited(struct graph_info_t *ginfo, int e);
void mark_link_trash(struct graph_info_t *ginfo, int e1, int e2);
void mark_link_visited(struct graph_info_t *ginfo, int e1, int e2);
void unmark_link_visited(struct graph_info_t *ginfo, int e1, int e2);
void unmark_edge_visited(struct graph_info_t *ginfo, int e);
void copy_static_info(struct graph_info_t *dest, struct graph_info_t *source);

int check_key_exist(khash_t(gint_int) *h, gint_t key);
void insert_key(khash_t(gint_int) *h, gint_t key);
void remove_key(khash_t(gint_int) *h, gint_t key);


int get_path(struct asm_graph_t *lg, int start_edge, int end_edge,
		int *middle_edge, int **path, int *path_len);
int get_best_middle_edge(struct asm_graph_t *lg, struct graph_info_t *ginfo);
void filter_edges(struct asm_graph_t *lg, struct graph_info_t *ginfo);
void cov_filter(struct asm_graph_t *lg, struct graph_info_t *ginfo);
void link_filter(struct asm_graph_t *lg, struct graph_info_t *ginfo);
void connection_filter(struct asm_graph_t *lg, struct graph_info_t *ginfo);
void bfs(struct asm_graph_t *lg, struct graph_info_t *ginfo, int start_edge,
		int **bfs_len);
int find_path_hao(struct asm_graph_t *lg, struct graph_info_t *ginfo, int u,
		int depth, int **path, int *path_len);
int check_simple_path(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int *old_path, int old_path_len, int **new_path,
		int *new_path_len);
void print_path(int *path, int path_len);
void find_middle_edge_candidates(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int u, int *path, int depth, int *mark);
void print_graph(struct asm_graph_t *lg, struct graph_info_t *ginfo);
void get_all_paths(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		struct path_info_t *pinfo);
void find_all_paths(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int u, int depth, int *cur_path, struct path_info_t *pinfo);
void path_info_init(struct path_info_t *path);
void path_info_push(struct path_info_t *pinfo, int *path, int len);
void path_info_destroy(struct path_info_t *pinfo);
#endif
