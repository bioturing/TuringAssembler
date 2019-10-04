#ifndef __BUILD_BRIDGE__
#define __BUILD_BRIDGE__
#define KSIZE_CHECK (lg->ksize + 6)
#define MIN_PATH_LENGTH 100
#define MIN_DEPTH_RATIO 0.2
#define MIN_OUTPUT_CONTIG_LEN 500
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
#define BRIDGE_LOCAL_NOT_FOUND 0
#define BRIDGE_TRIVIAL_BRIDGE 1
#define BRIDGE_MULTIPLE_PATH 2
#define BRIDGE_PATH_NOT_FOUND 3
#define N_BRIDGE_TYPE 4
#include <stdlib.h>
#include "assembly_graph.h"
#include "verbose.h"
#include "khash.h"
#include "map_contig.h"
#include "graph_search.h"

KHASH_MAP_INIT_STR(str_int, int);

static char *local_asm_result[4] = {"BRIDGE_LOCAL_NOT_FOUND", "BRIDGE_TRIVIAL_BRIDGE", "BRIDGE_MULTIPLE_PATH", "BRIDGE_PATH_NOT_FOUND"};

struct scaffold_record_t{
	int n_paths;
	int *path_lens;
	int **paths;
};

struct query_record_t{
	int *e1;
	int *e2;
	int *path_id;
	int process_pos;
	int n_process;
};

struct build_bridge_bundle_t{
	struct opt_proc_t *opt;
	struct asm_graph_t *g;
	struct query_record_t *query_record;
	struct scaffold_record_t *scaffold_record;
	pthread_mutex_t *query_lock;
	char **bridges;
	pthread_mutex_t *bridge_lock;
};

void get_local_edge_head(struct asm_graph_t g, struct asm_graph_t lg,
		int e, struct edge_map_info_t *emap);
void get_local_edge_tail(struct asm_graph_t g, struct asm_graph_t lg,
		int e, struct edge_map_info_t *emap);
int get_bridge(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, int e1, int e2, int *scaffolds,
		int n_scaff, char **res_seq, int *seq_len);

void combine_edges(struct asm_graph_t lg, int *path, int path_len, char **seq);
gint_t get_edge_code(gint_t u, gint_t v);

void join_seq(char **dest, char *source);
void sync_global_local_edge(struct asm_edge_t global, struct asm_edge_t local,
		struct subseq_pos_t global_pos, struct subseq_pos_t local_pos,
		int sync_type, char **res_seq);
void join_trivial_bridge(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_graph_t lg, struct edge_map_info_t *emap1,
		struct edge_map_info_t *emap2, char **res_seq);
void join_bridge_by_path(struct asm_edge_t e1, struct asm_edge_t e2,
		struct asm_graph_t lg, int *path, int path_len,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2,
		char **res_seq);
int try_bridging(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, int *scaffolds, int n_scaff,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2,
		char **res_seq, int *seq_len);
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
		int e1, int e2, float **scores, float **error);
void join_bridge_center_by_path(struct asm_graph_t *lg, int *path, int path_len,
		char **seq);
void unrelated_filter(struct asm_graph_t *g, struct edge_map_info_t *emap1,
		struct edge_map_info_t *emap2, int *scaffolds, int n_scaff,
		struct asm_graph_t *lg);
void join_bridge_dump(struct asm_edge_t e1, struct asm_edge_t e2,
		char **res_seq);
void join_bridge_no_path(struct asm_graph_t *g, struct asm_graph_t *lg,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2,
		char **res_seq);
void get_best_path(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, struct edge_map_info_t *emap1,
		struct edge_map_info_t *emap2, int *scaffolds, int n_scaff,
		int **path, int *path_len);
void cov_filter(struct asm_graph_t *g, struct asm_graph_t *lg,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2);
void connection_filter(struct asm_graph_t *g, struct asm_graph_t *lg,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2);
void print_log_edge_map(struct edge_map_info_t *emap1, struct edge_map_info_t *emap2);
void link_filter(struct opt_proc_t *opt, struct asm_graph_t *g, struct asm_graph_t *lg,
		struct edge_map_info_t *emap1, struct edge_map_info_t *emap2);
int check_degenerate_graph(struct asm_graph_t *g, struct asm_graph_t *lg,
		int e1, int e2);
void build_bridge(struct opt_proc_t *opt, FILE *f);
void *build_bridge_iterator(void *data);
void get_all_local_graphs(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct query_record_t *query);
void cleanup(struct opt_proc_t *opt);
#endif
