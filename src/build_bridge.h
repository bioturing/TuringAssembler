#ifndef __BUILD_BRIDGE__
#define __BUILD_BRIDGE__
#define MIN_MATCH_LENG 500
#define MIN_UNMATCHED_RATIO 0.005
#define MIN_RADIUS 2
#define MIN_EDGE_LENGTH 10000
#define MIN_EDGE_LENGTH_RATIO 0.8
#define MIN_RELATIVE_DEPTH_RATIO 0.2
#define MIN_DEPTH_RATIO 0.05
#define MIN_PATH_LENGTH 20
#define READ_GAP 500
#define BRIDGE_GAP 1000
#define MIN_READ_MAP_RATIO 0.95
#include <stdlib.h>
#include "assembly_graph.h"
#include "verbose.h"
#include "khash.h"
KHASH_MAP_INIT_INT64(gint_t_int, int);
KHASH_MAP_INIT_STR(str_int, int);
int match_head(struct asm_edge_t P, struct asm_edge_t T);
int match_tail(struct asm_edge_t P, struct asm_edge_t T);
void get_local_edge_id_head(struct asm_graph_t g, struct asm_edge_t e,
				int *edge_id, int *pos);
void get_local_edge_id_tail(struct asm_graph_t g, struct asm_edge_t e,
				int *edge_id, int *pos);
int check_simple_path(struct asm_graph_t lg, int start_edge, int end_edge);
/*int find_path(struct asm_graph_t lg, int u, int end_edge, int *visited,
		int *trace, int ban_edge1, int ban_edge2);*/
/*int find_path(struct asm_graph_t lg, gint_t u, int start_edge, int end_edge,
		khash_t(gint_t_int) *visited, int depth, int **path,
		int *path_leng);*/
int find_path(struct asm_graph_t lg, gint_t u, int start_edge, int end_edge,
		khash_t(gint_t_int) *visited, int *is_disable, int depth,
		int **path, int *path_leng);
void find_all_paths(struct asm_graph_t g, gint_t u, int start_edge, int end_edge,
		khash_t(gint_t_int) *vistied, int *is_disable, int depth,
		int *path, FILE *record, int *count_paths);
void get_path(struct asm_graph_t lg, int start_edge, int end_edge,
		int **path, int *path_leng);
/*void print_path(struct asm_graph_t lg, int start_edge, int end_node,
		int *trace);*/
void print_path(int *path, int path_leng);

int get_bridge(struct asm_graph_t *g, struct asm_graph_t *lg, int e1, int e2,
		uint32_t **ret_seq, uint32_t *seq_len);

void combine_edges(struct asm_graph_t lg, uint32_t **seq, uint32_t *leng,
			int head_leng, int head_pos, int tail_leng,
			int tail_pos, int *path, int path_leng);
float get_cov(struct asm_graph_t g, int edge);
gint_t get_edge_code(gint_t u, gint_t v);

void cov_filter(struct asm_graph_t g, int start_edge, int end_edge,
		int *is_disable);
int dfs(struct asm_graph_t g, int u, int start_edge, int end_edge, int *visited,
		int *is_disable);
void bfs(struct asm_graph_t g, int start_edge, int *is_disable, int **path_leng);
int reachable(struct asm_graph_t g, int start_edge, int end_edge,
		int *is_disable);
void connection_filter(struct asm_graph_t g, int start_edge, int end_edge,
		int *is_disable);
void filter_edges(struct asm_graph_t g, int start_edge, int end_edge,
		int **is_disable);
void print_graph(struct asm_graph_t g, int *is_disable, char *path);
char int_to_base(int x);

int is_read_on_edge(struct read_t r, struct asm_edge_t e, int left, int right,
		int ksize);
void copy_seq_to_str(char *str, uint32_t *seq, int leng);

void get_midair_bridge(struct asm_graph_t lg, int start_edge, int end_edge,
		int *is_disable, int *midair_edge);
#endif
