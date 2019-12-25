#include <string.h>

#include "assembly_graph.h"
#include "io_utils.h"
#include "kseq.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"
#include "khash.h"
#include "kmhash.h"

#include <zlib.h>

#if defined(_MSC_VER)

#define FORCE_INLINE	__forceinline

#include <stdlib.h>

#define ROTL32(x,y)	_rotl(x,y)
#define ROTL64(x,y)	_rotl64(x,y)

#define BIG_CONSTANT(x) (x)

// Other compilers

#else	// defined(_MSC_VER)

#define	FORCE_INLINE inline __attribute__((always_inline))

inline uint32_t rotl32(uint32_t x, int8_t r)
{
	return (x << r) | (x >> (32 - r));
}

static inline uint64_t rotl64(uint64_t x, int8_t r)
{
	return (x << r) | (x >> (64 - r));
}

#define	ROTL32(x,y)	rotl32(x,y)
#define ROTL64(x,y)	rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)

#endif // !defined(_MSC_VER)

KSEQ_INIT(gzFile, gzread);

static inline uint64_t getblock64(const uint64_t * p, int i)
{
	return p[i];
}

static inline uint64_t fmix64(uint64_t k)
{
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xff51afd7ed558ccd);
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
	k ^= k >> 33;

	return k;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  asm_fasta_edge_convert
 *  Description:  Dump the nucleotide sequence into the binary format
 * =====================================================================================
 */
static int asm_fasta_edge_convert(struct asm_graph_t *g, gint_t e, kseq_t *seq)
{
	uint32_t *p_holes, *l_holes, *bseq;
	uint32_t n_holes, i, j, k, c, last_c, m_seq;
	p_holes = NULL;
	l_holes = NULL;
	m_seq = 0x100;
	bseq = calloc(m_seq, sizeof(uint32_t));
	last_c = 0;
	n_holes = 0;
	for (i = k = 0; i < seq->seq.l; ++i) {
		c = nt4_table[(int)seq->seq.s[i]];
		if (c >= 4) {
			if (last_c >= 4) {
				++l_holes[j];
			} else {
				if (k == 0) {
					free(bseq);
					return 0;
				}
				p_holes = realloc(p_holes, (n_holes + 1) * sizeof(uint32_t));
				l_holes = realloc(l_holes, (n_holes + 1) * sizeof(uint32_t));
				j = n_holes;
				++n_holes;
				p_holes[j] = k - 1;
				l_holes[j] = 1;
			}
		} else {
			/* append new char */
			if (((k + 15) >> 4) >= m_seq) {
				uint32_t new_m = m_seq << 1;
				bseq = realloc(bseq, new_m * sizeof(uint32_t));
				memset(bseq + m_seq, 0, m_seq * sizeof(uint32_t));
				m_seq = new_m;
			}
			__binseq_set(bseq, k, c);
			++k;
		}
		last_c = c;
	}
	g->edges[e].seq_len = k;
	g->edges[e].seq = realloc(bseq, ((k + 15) >> 4) * sizeof(uint32_t));
	g->edges[e].count = 0;
	g->edges[e].n_holes = n_holes;
	g->edges[e].l_holes = l_holes;
	g->edges[e].p_holes = p_holes;
	return 1;
}

KHASH_MAP_INIT_INT64(khInt, uint64_t);

// shorthand way to get the key from hashtable or defVal if not found
#define kh_get_val(kname, hash, key, defVal) ({k=kh_get(kname, hash, key);(k!=kh_end(hash)?kh_val(hash,k):defVal);})
// // returns 0=replaced existing item, 1=bucket empty (new key), 2-adding element previously deleted
#define kh_set(kname, hash, key, val) ({int ret; k = kh_put(kname, hash,key,&ret); kh_value(hash,k) = val; ret;})

static unsigned char seq_nt6_table[256] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_adj_idx
 *  Description:  Find the index of an edge in the adjacencies list of one node
 *  Return     :  Return -1 is the edge doesn't exists
 * =====================================================================================
 */
static inline gint_t find_adj_idx(gint_t *adj, gint_t deg, gint_t id)
{
	gint_t i, ret;
	ret = -1;
	for (i = 0; i < deg; ++i) {
		if (adj[i] == id)
			ret = i;
	}
	return ret;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  fastg_hash_kmer1
 *  Description:  Hash the prefix length k of the sequence seq
 * =====================================================================================
 */
static uint64_t fastg_hash_kmer1(char *seq, int k, int l)
{
	uint64_t hash = MurmurHash3_x64_64((const uint8_t *) seq, k);
	return hash;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fastg_hash_kmer2
 *  Description:  Hash the suffix length k of sequence seq
 * =====================================================================================
 */
static uint64_t fastg_hash_kmer2(char *seq, int k, int l)
{
	int i, j = 0;
	char *s = malloc(k * sizeof(char));
	for (i = l - k; i < l; s[j++] = seq[i++]);
	uint64_t hash = MurmurHash3_x64_64((const uint8_t *) s, k);
	free(s);
	return hash;
}

static void insert_one_edge(struct asm_graph_t *g, gint_t u, gint_t v, gint_t e, kseq_t *seq)
{
	g->edges[e].source = u;
	g->edges[e].target = v;
	g->edges[e].count = 1;
	assert(asm_fasta_edge_convert(g, e, seq) && "Can not serialize sequence into edge!\n");
	g->nodes[u].adj = realloc(g->nodes[u].adj, sizeof(gint_t) * (g->nodes[u].deg + 1));
	g->nodes[u].adj[g->nodes[u].deg++] = e;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  load_asm_graph_fastg
 *  Description:  Load the assembly graph as fastg format into the skipping graph structure. 
 *  Haven't test for graph from other assembler except for MEGAHIT
 * =====================================================================================
 */
void load_asm_graph_fastg(struct asm_graph_t *g, const char *path, int ksize )
{
	g->ksize = ksize;
	g->aux_flag = 0;
	g->edges = NULL;
	g->nodes = NULL;
	g->n_v = g->n_e = 0;
	gint_t p_id = 0, q_id = 0;
	int c;
	khint_t k;
	int absent;
	gzFile fp = gzopen(path, "r");
	gzFile fp_edge = gzopen(path, "r");
	int v_cnt = 0;
	int edges_cnt = 0;
	
	// Keeping inserted kmer node
	khash_t(khInt) *map_v = kh_init(khInt);

	if (!fp)
		log_error("Load graph fastg: Unable to open file [%s] to read", path);
	kseq_t *seq = kseq_init(fp);

	while (kseq_read(seq) >= 0) {
		g->edges = realloc(g->edges, sizeof(struct asm_edge_t) * (g->n_e + 1));
		// Hash kmer1
		uint64_t hash1 = fastg_hash_kmer1(seq->seq.s, ksize, seq->seq.l);
		uint64_t hash2 = fastg_hash_kmer2(seq->seq.s, ksize, seq->seq.l);

		//log_info("Hashed kmer1: %u, kmer2: %u", hash1, hash2);
		char *p, *s = seq->name.s;
		int is_comp;
		for (p = s; *p && *p != ':' && *p != ';'; ++p);
		c = *p, *p = 0;
		is_comp = (p > s && *(p-1) == '\''); // if we are looking at a complement segment
		if (is_comp) *(p-1) = 0;

		if (kh_get(khInt, map_v, hash1) == kh_end(map_v)) {
			g->nodes = realloc(g->nodes, sizeof(struct asm_node_t) * (g->n_v + 1));
			k = kh_put(khInt, map_v, hash1, &absent);
			kh_value(map_v, k) = g->n_v;
			g->nodes[g->n_v].deg = 0;
			g->nodes[g->n_v].adj = NULL;
			++g->n_v;
		}
		if (kh_get(khInt, map_v, hash2) == kh_end(map_v)) {
			g->nodes = realloc(g->nodes, sizeof(struct asm_node_t) * (g->n_v + 1));
			k = kh_put(khInt, map_v, hash2, &absent);
			kh_value(map_v, k) = g->n_v;
			g->nodes[g->n_v].deg = 0;
			g->nodes[g->n_v].adj = NULL;
			++g->n_v;
		}

		gint_t u = kh_get_val(khInt, map_v, hash1, (uint64_t)-1 );
		gint_t v = kh_get_val(khInt, map_v, hash2, (uint64_t)-1 );

		if (find_adj_idx(g->nodes[u].adj, g->nodes[u].deg, g->n_e) == -1) {
			insert_one_edge(g, u, v, g->n_e, seq);
		}
		++g->n_e;
	}	
	
	log_info("Number of nodes: %d", (int) g->n_v);
	log_info("Number of edges: %d", (int) g->n_e);

	// Set the reverse complement of node and edge
	int i;
	for (i = 0; i < g->n_e; ++i) {
		if (!(i & 1)) {
			gint_t e_rc = i + 1; // 1 and 2 is the reverse complement of each other
			gint_t src = g->edges[i].source; 
			gint_t dst = g->edges[i].target;
			g->nodes[src].rc_id = g->edges[e_rc].target; // Set the reverse complement of node
			g->nodes[g->edges[e_rc].target].rc_id = src; 
			g->nodes[dst].rc_id = g->edges[e_rc].source;
			g->nodes[g->edges[e_rc].source].rc_id = dst;
			g->edges[i].rc_id = e_rc; // Set the reverse complement of edge
			g->edges[e_rc].rc_id = i;
		}
	}

	kseq_destroy(seq);
	gzclose(fp);
}
