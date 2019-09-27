#ifndef __MAP_CONTIG__
#define __MAP_CONTIG__
#include "assembly_graph.h"
#include "khash.h"
#include "helper.h"
#include "verbose.h"
#include "kmer_hash.h"
#define KSIZE 100
#define WINDOW_SIZE 1000
#define POINT_HIGH_THRESH 0.9

struct map_contig_t{
	struct asm_edge_t global_edge;
	struct asm_graph_t local_graph;
	khash_t(kmer_int) **kmers;
	int n_candidates;
	int pos;
	int best_match;
	int *is_match;
};

struct subseq_pos_t{
	int start;
	int end;
};

struct edge_map_info_t{
	int gl_e;
	int lc_e;
	struct subseq_pos_t gpos;
	struct subseq_pos_t lpos;
};

void init_map_contig(struct map_contig_t *mct, struct asm_edge_t global_edge,
		struct asm_graph_t local_graph);
void init_local_kmers(struct map_contig_t *mct);
void get_all_seq_kmers(char *seq, int len, khash_t(kmer_int) **kmers);
void get_all_edge_kmers(struct asm_edge_t edge, khash_t(kmer_int) **kmers);
int find_match_from_pos(struct map_contig_t *mct);
int count_match_kmer(khash_t(kmer_int) *first, khash_t(kmer_int) *second);
int find_match(struct map_contig_t *mct);
void map_contig_destroy(struct map_contig_t *mct);
void advance_pos(struct map_contig_t *mct);
void get_match_pos(struct map_contig_t *mct, struct subseq_pos_t *global,
		struct subseq_pos_t *local);
void get_global_match_pos(struct map_contig_t *mct, struct subseq_pos_t *pos);
void get_local_match_pos(struct map_contig_t *mct, struct subseq_pos_t *global,
		struct subseq_pos_t *local);
int get_next_len_global(struct map_contig_t *mct, int pos);
int get_next_len_local(struct map_contig_t *mct, int pos);
uint64_t get_one_seq_kmer_hash(char *seq);
int check_good_match(int point, float thresh);
void add_kmer(uint64_t hash, khash_t(kmer_int) *kmers);
void remove_kmer(uint64_t hash, khash_t(kmer_int) *kmers);
int check_stop(struct map_contig_t *mct);
#endif
