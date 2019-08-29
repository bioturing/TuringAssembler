#ifndef __BUILD_BRIDGE__
#define __BUILD_BRIDGE__
#define MIN_OUTPUT_CONTIG_LEN 1000
#define MIN_MATCH_LENG 4000
#define MATCH_THRESH 8000
#define MIN_UNMATCHED_RATIO 0.005
#define MIN_RADIUS 0
#define MIN_EDGE_LENGTH 10000
#define MIN_EDGE_LENGTH_RATIO 0.8
#define READ_GAP 500
#define BRIDGE_GAP 1000
#define MIN_READ_MAP_RATIO 0.95
#define POS_IGNORED 1e9
#define NO_PATH_FOUND 0
#define TRIVIAL_BRIDGE 1
#define SINGLE_PATH 2
#define MULTIPLE_PATH 3
#define N_BRIDGE_TYPE 4
#define SYNC_KEEP_GLOBAL 0
#define SYNC_KEEP_LOCAL 1
#include <stdlib.h>
#include "assembly_graph.h"
#include "verbose.h"
#include "khash.h"
#include "map_contig.h"
#include "graph_search.h"
KHASH_MAP_INIT_STR(str_int, int);

void get_local_edge_head(struct asm_graph_t g, struct asm_graph_t lg,
		struct asm_edge_t e, int *edge_id, struct subseq_pos_t *gpos,
		struct subseq_pos_t *lpos);
void get_local_edge_tail(struct asm_graph_t g, struct asm_graph_t lg,
		struct asm_edge_t e, int *edge_id, struct subseq_pos_t *gpos,
		struct subseq_pos_t *lpos);
int get_bridge(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, int e1, int e2, int pre_e1, int next_e2,
		char **res_seq, int *seq_len);

void combine_edges(struct asm_graph_t lg, int *path, int path_len, char **seq);
gint_t get_edge_code(gint_t u, gint_t v);

void join_seq(char **dest, char *source);
void sync_global_local_edge(struct asm_edge_t global, struct asm_edge_t local,
		struct subseq_pos_t global_pos, struct subseq_pos_t local_pos,
		int sync_type, char **res_seq);
void join_trivial_bridge(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_graph_t lg, int local_edge,
		struct subseq_pos_t gpos1, struct subseq_pos_t lpos1,
		struct subseq_pos_t gpos2, struct subseq_pos_t lpos2,
		char **res_seq);
void join_bridge_by_path(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_graph_t lg, int *path, int path_len, 
		struct subseq_pos_t gpos1, struct subseq_pos_t lpos1,
		struct subseq_pos_t gpos2, struct subseq_pos_t lpos2,
		char **res_seq);
void get_contig_from_scaffold_path(struct opt_proc_t *opt, struct asm_graph_t *g,
		int *path, int path_len, char **contig);
int try_bridging(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, int e1, int e2, int pre_e1, int next_e2,
		int lc_e1, int lc_e2, struct subseq_pos_t gpos1,
		struct subseq_pos_t lpos1, struct subseq_pos_t gpos2,
		struct subseq_pos_t lpos2, char **res_seq, int *seq_len);
void join_complex_path(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_edge_t lc_e1, struct asm_edge_t lc_e2,
		struct subseq_pos_t gpos1, struct subseq_pos_t lpos1,
		struct subseq_pos_t gpos2, struct subseq_pos_t lpos2,
		char **res_seq);
void join_middle_edge(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_edge_t lc_e1, struct asm_edge_t lc_e2,
		struct subseq_pos_t gpos1, struct subseq_pos_t lpos1,
		struct subseq_pos_t gpos2, struct subseq_pos_t lpos2,
		struct asm_edge_t middle, char **res_seq);
void get_path_scores(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, struct path_info_t *pinfo,
		int e1, int e2, float **scores);
void join_bridge_center_by_path(struct asm_graph_t *lg, int *path, int path_len,
		char **seq);
void unrelated_filter(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_edge_t pre_e1, struct asm_edge_t next_e2,
		struct asm_graph_t *lg, struct graph_info_t *ginfo);
#endif
