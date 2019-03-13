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

/******************** Build graph ultilities **********************************/
/* Build graph using kmer with length smaller than 32 */
void k31_build0(struct opt_count_t *opt, int ksize, struct asm_graph_t *g0);
void build_asm_graph_from_k31(struct opt_count_t *opt, int ksize,
			struct k31hash_t *kmer_hash, struct asm_graph_t *ret_g);
/* Build graph using kmer with length from 33 to 63 */
void k63_build0(struct opt_count_t *opt, int ksize, struct asm_graph_t *g0);
void build_asm_graph_from_k63(struct opt_count_t *opt, int ksize,
			struct k63hash_t *kmer_hash, struct asm_graph_t *ret_g);

/********************* Utilities for edges manipulating ***********************/

/* Function signature: float __get_edge_cov(struct asm_edge_t *e, int ksize) */
#define __get_edge_cov(e, ksize) ((e)->count * 1.0 /			\
						(((e)->n_holes + 1) * (ksize)))
/* Write plain nucleotide sequence to buffer seq */
void dump_edge_seq(char **seq, uint32_t *m_seq, struct asm_edge_t *e);
/* Estimate coverage of 1 walk on genome */
float get_genome_coverage(struct asm_graph_t *g);
/* Copy sequence, gap and kmer count information from src to dst */
void asm_clone_edge(struct asm_edge_t *dst, struct asm_edge_t *src);
/* Copy edge and do reversing complement */
void asm_clone_reverse(struct asm_edge_t *dst, struct asm_edge_t *src);
/* append sequence, gap and kmer count (please note that this function does not
 * check for consistence, including source and target node)
 */
void asm_append_edge(struct asm_edge_t *dst, struct asm_edge_t *src,
							uint32_t overlap);

/********************** Utilities for graph manipulating **********************/

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

void test_asm_graph(struct asm_graph_t *g);

#endif  /* __ASSEMBLY_GRAPH_H__ */
