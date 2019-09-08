#ifndef __GRAPH_SEARCH__
#define __GRAPH_SEARCH__
#define MIN_RELATIVE_COV_RATIO 0.2
#define SIMPLE_PATH 0
#define COMPLEX_PATH 1
#define MAX_PATH_COUNT 100
#include "khash.h"
#include "assembly_graph.h"
#include "resolve.h"
#include "map_contig.h"
#include "kmer_hash.h"
#include <stdio.h>
#define PATH_NOT_FOUND -1
KHASH_MAP_INIT_INT64(gint_int, int);
struct graph_info_t{
	struct asm_graph_t *g;
	int lc_e1;
	int lc_e2;
	int *is_edge_trash;
	int *edge_vst_count;
	int *edge_max_vst;
	khash_t(gint_int) *is_link_trash;
	khash_t(gint_int) *is_link_vst;
};

struct path_info_t{
	int n_paths;
	int m_paths;
	int *path_lens;
	int **paths;
};

void graph_info_init(struct asm_graph_t *lg, struct graph_info_t *ginfo,
			int lc_e1, int lc_e2);
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


int get_path(struct asm_graph_t *lg, int lc_e1, int lc_e2,
		int *middle_edge, int **path, int *path_len);
int get_best_middle_edge(struct asm_graph_t *lg, struct graph_info_t *ginfo);
void bfs(struct asm_graph_t *lg, struct graph_info_t *ginfo, int lc_e1,
		int **bfs_len);
int find_path_hao(struct asm_graph_t *lg, struct graph_info_t *ginfo, int u,
		int depth, int **path, int *path_len);
int check_simple_path(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int *old_path, int old_path_len, int **new_path,
		int *new_path_len);
void print_path(int *path, int path_len);
void find_middle_edge_candidates(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int u, int *path, int depth, int *mark);
void print_graph(struct asm_graph_t *lg, int lc_e1, int lc_e2);
void get_all_paths(struct asm_graph_t *g, struct asm_graph_t *lg,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2,
		struct path_info_t *pinfo);
void get_all_paths_kmer_check(struct asm_graph_t *g, struct asm_graph_t *lg,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2,
		struct path_info_t *pinfo, int ksize, khash_t(kmer_int) *h);
void find_all_paths(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int u, int depth, int *cur_path, struct path_info_t *pinfo);
void find_all_paths_kmer_check(struct asm_graph_t *lg, struct graph_info_t *ginfo,
		int u, int depth, int *cur_path, struct path_info_t *pinfo,
		int ksize, khash_t(kmer_int) *h);
void path_info_init(struct path_info_t *path);
void path_info_push(struct path_info_t *pinfo, int *path, int len);
void path_info_destroy(struct path_info_t *pinfo);
#endif
