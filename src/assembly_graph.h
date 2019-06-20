#ifndef __ASSEMBLY_GRAPH_H__
#define __ASSEMBLY_GRAPH_H__

#include <stdint.h>

#include "attribute.h"
#include "barcode_hash.h"

#define ASM_HAVE_BARCODE		0x1
#define ASM_HAVE_READPAIR		0x2

#define ASM_BUILD_BARCODE		0x1
#define ASM_BUILD_READPAIR		0x2
#define ASM_BUILD_COVERAGE		0x4

struct asm_node_t {
	gint_t rc_id;		/* reverse complement link */
	gint_t deg;		/* out degree */
	gint_t *adj;		/* list of out edges */
};

struct asm_edge_t {
	uint64_t count;
	uint32_t *seq;		/* only store contigs */
	uint32_t seq_len;	/* length of "seq" */

	uint32_t n_holes;	/* number of gaps */
	/* the i-th gap starts right after the position p_holes[i]-th of seq and
	 * has length equal to l_holes[i]
	 */
	uint32_t *p_holes;	/* positions of holes */
	uint32_t *l_holes;	/* length of holes */

	gint_t source;		/* start node */
	gint_t target;		/* end node */
	gint_t rc_id;		/* reverse complement link */
	pthread_mutex_t lock;	/* lock for build/mapping process */
	struct barcode_hash_t barcodes;		/* mapped barcode */
	struct barcode_hash_t *mate_barcodes;
	gint_t *mate_contigs;
	int n_mate_contigs;
	// struct barcode_hash_t mate_contigs;	/* list of mate contigs (build process only) */
	// gint_t best_mate_contigs;		/* best mate contigs picker */
};

struct asm_graph_t {
	int ksize;			/* ksize of nodes */
	int bin_size;			/* deprecated */
	uint32_t aux_flag;		/* aux flag, marker for storing BARCODE of READPAIR */
	gint_t n_v, n_e;		/* number of nodes, number of edges */

	struct asm_node_t *nodes;	/* list of nodes */
	struct asm_edge_t *edges;	/* list of edges */
};

#define MIN_FRAG_LEN			200
#define MIN_UNIQUE_BARCODE		100
#define MIN_NOTICE_LEN			200
#define MIN_CONNECT_SIZE		500
#define MAX_TIPS_LEN			250
#define TIPS_THRESHOLD			5.0
#define TIPS_RATIO_THRESHOLD		0.1

/* Add barcode upto prefix length */
#define MIN_CONTIG_BARCODE		3000
/* Only add and use barcode for contig with length minimum */
#define MIN_LONG_CONTIG			1000
#define MIN_CONTIG_READPAIR		500

#define MAX_PAIR_LEN			700
#define MAX_MOLECULE_LEN		30000

#define MIN_BARCODE_COUNT		150
#define MIN_READPAIR_COUNT		20

#define MIN_BARCODE_RATIO		0.1
#define MIN_SUB_BARCODE_RATIO		0.04

/************************* Build graph ultilities *****************************/
/******************************************************************************/

/* Build graph using KMC module for kmer counting */
void graph_build_KMC(struct opt_proc_t *opt, int ksize, struct asm_graph_t *g);

/*************************** Barcoding ultilities *****************************/
/******************************************************************************/

/* get the barcode share ratio between two contigs */
double get_barcode_ratio(struct asm_graph_t *g, gint_t e1, gint_t e2);
double get_barcode_ratio_unique(struct asm_graph_t *g, gint_t e1, gint_t e2);
/* construct the barcode map */
void construct_aux_information(struct opt_proc_t *opt, struct asm_graph_t *g, uint32_t aux_build);

/********************* Utilities for edges manipulating ***********************/
/******************************************************************************/

struct cov_range_t {
	int lo;
	int hi;
};

#define __ratio_greater(fcov1, fcov2) ((fcov1) > 2.0 * (fcov2))
#define __coverage_range_intersect(rcov1, rcov2) ((rcov1).lo <= (rcov2).hi && (rcov2).lo <= (rcov1).hi)
#define __check_coverage(fcov1, fcov2, rcov1, rcov2)			\
	(__coverage_range_intersect(rcov1, rcov2) && __abs((fcov1) - (fcov2)) < 0.3)

/* Function signature:
 * int __int_cov_ratio(float cov1, float cov2);
 */
#define __int_cov_ratio(cov1, cov2) ((int)((cov1) / (cov2) + 0.49999999))
#define __binseq_set(seq, k, c) ((seq)[(k) >> 4] |= (uint32_t)(c) << (((k) & 15) << 1))
/* Function signature:
 * uint32_t __binseq_get(uint32_t *seq, uint32_t k);
 * Extract the numberic nucleotide at position k-th on seq
 */
#define __binseq_get(seq, k) (((seq)[(k) >> 4] >> (((k) & 15) << 1)) & (uint32_t)3)
/* Function signature:
 * float __get_edge_cov(struct asm_edge_t *e, int ksize);
 * */
#define __get_edge_cov(e, ksize) ((e)->count * 1.0 /			\
				((e)->seq_len - ((e)->n_holes + 1) * (ksize)))
static inline struct cov_range_t get_edge_cov_range(struct asm_graph_t *g, gint_t e, double uni_cov)
{
	double fcov = __get_edge_cov(g->edges + e, g->ksize) / uni_cov;
	int icov = (int)fcov;
	if (fcov + EPS < icov + 0.25)
		return (struct cov_range_t){icov, icov};
	else if (fcov + EPS > icov + 0.25 && fcov + EPS < icov + 0.75)
		return (struct cov_range_t){icov, icov + 1};
	else
		return (struct cov_range_t){icov + 1, icov + 1};
}

static inline struct cov_range_t convert_cov_range(double fcov)
{
	int icov = (int)fcov;
	if (fcov + EPS < icov + 0.25)
		return (struct cov_range_t){icov, icov};
	else if (fcov + EPS > icov + 0.25 && fcov + EPS < icov + 0.75)
		return (struct cov_range_t){icov, icov + 1};
	else
		return (struct cov_range_t){icov + 1, icov + 1};
}
/* Estimate coverage of 1 walk on genome */
double get_genome_coverage(struct asm_graph_t *g);
/* Copy sequence, gap and kmer count information from src to dst */
void asm_clone_edge(struct asm_graph_t *g, gint_t dst, gint_t src);
gint_t asm_create_node(struct asm_graph_t *g);
gint_t asm_create_clone_edge(struct asm_graph_t *g, gint_t src);
void asm_append_seq(struct asm_edge_t *dst, struct asm_edge_t *src, uint32_t overlap);
void asm_clean_edge(struct asm_graph_t *g, gint_t e);
void asm_clone_seq(struct asm_edge_t *dst, struct asm_edge_t *src);
void asm_clone_seq_reverse(struct asm_edge_t *dst, struct asm_edge_t *src);
void asm_unroll_loop_forward(struct asm_graph_t *g, gint_t e1, gint_t e2);

void asm_join_edge(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
					gint_t e2, gint_t e_rc2);
void asm_join_edge_with_gap(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
				gint_t e2, gint_t e_rc2, uint32_t gap_size);
void asm_join_edge3(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
	gint_t e2, gint_t e_rc2, gint_t e3, gint_t e_rc3, uint64_t added_count);
void asm_join_edge_loop(struct asm_graph_t *g, gint_t e1, gint_t e_rc1,
			gint_t e2, gint_t e_rc2, uint64_t added_count);
void asm_graph_destroy(struct asm_graph_t *g);
/* only set the link to the edge to -1 */
void asm_remove_edge(struct asm_graph_t *g, gint_t e);
/* get the real length of the edge */
uint32_t get_edge_len(struct asm_edge_t *e);
/* get the index of the longest edge */
gint_t get_longest_edge(struct asm_graph_t *g);
void asm_duplicate_edge_seq(struct asm_graph_t *g, gint_t e, int cov);
void asm_duplicate_edge_seq2(struct asm_graph_t *g, gint_t e1, gint_t e2, int cov);
void asm_join_edge_loop_reverse(struct asm_graph_t *g, gint_t e1, gint_t e2,
				gint_t e_rc2, gint_t e_rc1);

/********************** Utilities for graph manipulating **********************/
/******************************************************************************/

/* Write confidence contigs/scaffolds to fasta file */
void write_fasta(struct asm_graph_t *g, const char *path);
/* Write graph as gfa format */
void write_gfa(struct asm_graph_t *g, const char *path);
/* Save graph topology as binary format */
void save_asm_graph(struct asm_graph_t *g, const char *path);
/* Load graph, auto detect saved type */
void load_asm_graph(struct asm_graph_t *g, const char *path);
/* Load fasta file and build graph base on reported contigs/scaffold */
void load_asm_graph_fasta(struct asm_graph_t *g, const char *path, int ksize);

/************************** Testing functions *********************************/
/******************************************************************************/

void test_asm_graph(struct asm_graph_t *g);
void print_test_barcode_edge(struct asm_graph_t *g, gint_t e1, gint_t e2);
void print_test_pair_end(struct asm_graph_t *g, gint_t e);
void print_test_barcode_superior(struct asm_graph_t *g, gint_t e1,
						gint_t e2, gint_t e2a);

#endif  /* __ASSEMBLY_GRAPH_H__ */
