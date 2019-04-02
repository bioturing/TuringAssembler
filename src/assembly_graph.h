#ifndef __ASSEMBLY_GRAPH_H__
#define __ASSEMBLY_GRAPH_H__

#include <stdint.h>

#include "attribute.h"
#include "k31hash.h"
#include "k63hash.h"

struct asm_node_t {
	gint_t rc_id;
	gint_t deg;
	gint_t *adj;
};

struct asm_edge_t {
	uint64_t count;
	uint32_t *seq;		/* only store contigs */
	uint32_t seq_len;	/* length of "seq" */

	uint32_t n_holes;	/* number of gaps */
	/* the i-th gap starts right after the position p_holes[i]-th of seq and
	 * has length equal to l_holes[i]
	 */
	uint32_t *p_holes;
	uint32_t *l_holes;

	gint_t source;
	gint_t target;
	gint_t rc_id;
	struct barcode_hash_t *bucks;
	pthread_mutex_t lock;
};

struct asm_graph_t {
	int ksize;
	int bin_size;
	gint_t n_v, n_e;

	struct asm_node_t *nodes;
	struct asm_edge_t *edges;
};

#define MIN_CONNECT_SIZE				250
#define MIN_CONTIG_LEN					100
#define MAX_TIPS_LEN					250
#define TIPS_THRESHOLD			5.0
#define TIPS_RATIO_THRESHOLD		0.1

/************************* Build graph ultilities *****************************/
/******************************************************************************/

/* Build graph using kmer with length smaller than 32 */
void k31_build0(struct opt_count_t *opt, int ksize, struct asm_graph_t *g0);
void build_asm_graph_from_k31(struct opt_count_t *opt, int ksize,
			struct k31hash_t *kmer_hash, struct asm_graph_t *ret_g);
/* Build graph using kmer with length from 33 to 63 */
void k63_build0(struct opt_count_t *opt, int ksize, struct asm_graph_t *g0);
void build_asm_graph_from_k63(struct opt_count_t *opt, int ksize,
			struct k63hash_t *kmer_hash, struct asm_graph_t *ret_g);

/*************************** Barcoding ultilities *****************************/
/******************************************************************************/

/* Test if 2 edges share enough barcode */
int test_edge_barcode(struct asm_graph_t *g, gint_t e1, gint_t e2);
/* print matrix of shared barcodes */
void print_test_barcode_edge(struct asm_graph_t *g, gint_t e1, gint_t e2);
double get_barcode_ratio(struct asm_graph_t *g, gint_t e1, gint_t e2);

/********************* Utilities for edges manipulating ***********************/
/******************************************************************************/

/* Function signature:
 * int __int_cov_ratio(float cov1, float cov2);
 */
#define __int_cov_ratio(cov1, cov2) ((int)((cov1) / (cov2) + 0.49999999))
/* Function signature:
 * uint32_t __binseq_get(uint32_t *seq, uint32_t k);
 * Extract the numberic nucleotide at position k-th on seq
 */
#define __binseq_get(seq, k) (((seq)[(uint32_t)(k) >> 4] >>		\
		(((uint32_t)(k) & 15) << 1)) & (uint32_t)3)
/* Function signature:
 * float __get_edge_cov(struct asm_edge_t *e, int ksize);
 * */
#define __get_edge_cov(e, ksize) ((e)->count * 1.0 /			\
				((e)->seq_len - ((e)->n_holes + 1) * (ksize)))
/* Write plain nucleotide sequence to buffer seq */
gint_t dump_edge_seq(char **seq, uint32_t *m_seq, struct asm_edge_t *e);
/* Estimate coverage of 1 walk on genome */
float get_genome_coverage(struct asm_graph_t *g);
/* Copy sequence, gap and kmer count information from src to dst */
void asm_clone_edge(struct asm_edge_t *dst, struct asm_edge_t *src);
void asm_clone_edge2(struct asm_graph_t *g, gint_t dst, gint_t src);
/* Copy edge and do reversing complement */
void asm_clone_reverse(struct asm_edge_t *dst, struct asm_edge_t *src);
/* append edge, check for integrity */;
void asm_append_edge(struct asm_edge_t *dst, struct asm_edge_t *src,
							uint32_t overlap);
/* append sequence and gap (please note that this function does not
 * check for consistence, including source and target node)
 */
void asm_append_edge_seq(struct asm_edge_t *dst, struct asm_edge_t *src,
							uint32_t overlap);
void asm_append_edge_seq2(struct asm_graph_t *g, gint_t e1, gint_t e2);
void asm_join_edge_with_gap(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
		gint_t e2, gint_t e_rc2, uint32_t gap_size, uint64_t gap_count);
void asm_join_edge3(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
	gint_t e2, gint_t e_rc2, gint_t e3, gint_t e_rc3, uint64_t added_count);
void asm_append_seq_with_gap(struct asm_edge_t *dst, struct asm_edge_t *src,
							uint32_t gap_size);
/* clean local data */
void asm_clean_edge_seq(struct asm_edge_t *e);
/* only set the link to the edge to -1 */
void asm_remove_edge(struct asm_graph_t *g, gint_t e);
/* get the real length of the edge */
uint32_t get_edge_len(struct asm_edge_t *e);
/* get the index of the longest edge */
gint_t get_longest_edge(struct asm_graph_t *g);

/********************** Utilities for graph manipulating **********************/
/******************************************************************************/

/* Write confidence contigs/scaffolds to fasta file */
void write_fasta(struct asm_graph_t *g, const char *path);
/* Write graph as gfa format */
void write_gfa(struct asm_graph_t *g, const char *path);
/* Save graph topology as binary format */
void save_asm_graph_simple(struct asm_graph_t *g, const char *path);
/* Save graph topology and barcode info as binary format */
void save_asm_graph_barcode(struct asm_graph_t *g, const char *path);
/* Load graph, auto detect saved type */
void load_asm_graph(struct asm_graph_t *g, const char *path);
/* Load graph, only load graph topology */
void load_asm_graph_simple(struct asm_graph_t *g, const char *path);
/* Load graph and barcode info */
void load_asm_graph_complex(struct asm_graph_t *g, const char *path);
/* Construct barcode information from reads */
void construct_barcode_map(struct asm_graph_t *g, struct opt_build_t *opt);

/************************** Testing functions *********************************/
/******************************************************************************/

void test_asm_graph(struct asm_graph_t *g);
void print_test_barcode_edge(struct asm_graph_t *g, gint_t e1, gint_t e2);

#endif  /* __ASSEMBLY_GRAPH_H__ */
