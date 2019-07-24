#include <stdlib.h>
#include <string.h>

#include "assembly_graph.h"
#include "barcode_hash.h"
#include "fastq_producer.h"
#include "io_utils.h"
#include "kmer_build.h"
#include "kseq.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"
#include "../include/kmc_skipping.h"

KSEQ_INIT(gzFile, gzread);

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

static int add_one_edge(struct asm_graph_t *g, gint_t id, kseq_t *seq, int added)
{
	if (id > g->n_e) {
		// +1 is used for the reverse complement of the new edge
		g->edges = realloc(g->edges, (id + 1) * sizeof(struct asm_edge_t));
		g->n_e = id + 1;
	} 

	//add an edge, not a neighbor
	if (seq != NULL) {
		if (!asm_fasta_edge_convert(g, id - 1, seq))
			return 1;
		asm_clone_seq_reverse(g->edges + id, g->edges + id - 1);
	}

	if (seq != NULL && added)
		return 0;

	g->nodes = realloc(g->nodes, (g->n_v + 4) * sizeof(struct asm_node_t));
	// edge id starts from zero
	g->edges[id - 1].rc_id = id;
	g->edges[id].rc_id = id - 1;

	g->edges[id - 1].source = g->n_v;
	g->edges[id - 1].target = g->n_v + 1;
	g->edges[id - 1].count = 0;
	g->nodes[g->n_v].adj = malloc(sizeof(gint_t));
	g->nodes[g->n_v].adj[0] = id - 1;
	g->nodes[g->n_v].deg = 1;
	g->nodes[g->n_v + 1].adj = NULL;
	g->nodes[g->n_v + 1].deg = 0;

	g->edges[id].source = g->n_v + 2;
	g->edges[id].target = g->n_v + 3;
	g->edges[id].count = 0;
	g->nodes[g->n_v + 2].adj = malloc(sizeof(gint_t));
	g->nodes[g->n_v + 2].adj[0] = id;
	g->nodes[g->n_v + 2].deg = 1;
	g->nodes[g->n_v + 3].adj = NULL;
	g->nodes[g->n_v + 3].deg = 0;

	g->nodes[g->n_v].rc_id = g->n_v + 3;
	g->nodes[g->n_v + 3].rc_id = g->n_v;
	g->nodes[g->n_v + 1].rc_id = g->n_v + 2;
	g->nodes[g->n_v + 2].rc_id = g->n_v + 1;
	g->n_v += 4;
	return 0;
}

static void add_one_connection(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
	// Short-hand for variable
	gint_t e1_target = g->edges[e1].target;
	gint_t e1rc = g->edges[e1].rc_id;
	// To avoid re-add a reverse complement connection
	if (e1 == 5651 && e2 == 5173){
		__VERBOSE("Herer");
	}
	if (find_adj_idx(g->nodes[e1_target].adj, g->nodes[e1_target].deg, e2) == -1) {
		// allocate 1 more room for target of e1
		g->nodes[e1_target].adj = realloc(g->nodes[e1_target].adj, (g->nodes[e1_target].deg + 1) * sizeof(gint_t));
		// add e2 as an out-going edge of e1_target
		g->nodes[e1_target].adj[g->nodes[e1_target].deg++] = e2;
		// e2.source is useless now, we must decrease its degree by 1
		--g->nodes[g->edges[e2].source].deg;
		// Assign new source for e2
		g->edges[e2].source = e1_target;
		// Do the reverse complement thing
		g->edges[g->edges[e2].rc_id].target = g->nodes[e1_target].rc_id;
		//g->nodes[g->edges[e1].target].rc_id = g->edges[e1rc].source;
		//g->nodes[g->edges[e1rc].source].rc_id = g->edges[e1].target;
	} else{
		__VERBOSE("Reverse complement!\n");
	}
	gint_t e1rc_src = g->edges[g->edges[e1].rc_id].source;
	assert(find_adj_idx(g->nodes[e1_target].adj, g->nodes[e1_target].deg, e2) != -1 &&
			find_adj_idx(g->nodes[e1rc_src].adj, g->nodes[e1rc_src].deg, g->edges[e1].rc_id) != -1);
}

void load_asm_graph_fastg(struct asm_graph_t *g, const char *path, int ksize )
{
	g->ksize = ksize;
	g->aux_flag = 0;
	g->edges = NULL;
	g->nodes = NULL;
	g->n_v = g->n_e = 0;
	gint_t p_id = 0, q_id = 0;
	int c;
	char *p_i, *q_i, *p, *q;
	p_i = 0;
	q_i = 0;
	int *added = calloc(1, sizeof(int));
	gzFile fp = gzopen(path, "r");
	if (!fp)
		__ERROR("Unable to open file [%s] to read", path);
	kseq_t *seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		/* add new edge */
		if ((int)seq->seq.l < ksize)
			continue;
		char *s = seq->name.s;
		int is_comp;
		for (p = s; *p && *p != ':' && *p != ';'; ++p);
		c = *p, *p = 0;
		is_comp = (p > s && *(p-1) == '\''); // if we are looking at a complement segment
		if (is_comp) *(p-1) = 0;
		for(p_i = p; *p_i != '_'; --p_i);
		p_id = atoi(p_i + 1);
		if (!is_comp) {
			add_one_edge(g, p_id, seq, g->n_e > 0 && added[p_id]);
		}
		if (p_id + 1 == g->n_e)
			added = realloc(added, sizeof(int) * g->n_e);
		added[p_id] = 1;
		added[p_id - 1] = 1;
		if (c == ':') { // have neighbors
			q = p + 1;
			do {
				int is_comp2 = 0;
				for (p = q; *p && *p != ',' && *p != ';'; ++p);
				c = *p, *p = 0;
				is_comp2 = (p > q && *(p-1) == '\'');
				if (is_comp2) *(p-1) = 0;
				for(q_i = p; *q_i != '_'; --q_i);
				q_id = atoi(q_i + 1);
				if (q_id > p_id) {
					if (q_id > g->n_e || !added[q_id]) {
						add_one_edge(g, q_id , NULL, 0);
						added = realloc(added, sizeof(int) * g->n_e);
					}
					added[q_id] = 1;
					added[q_id - 1] = 1;
				}
				add_one_connection(g, p_id - 1 + is_comp, q_id - 1 + is_comp2);
				q = p + 1;
			} while (c != 0 && c != ';');
		}
		__VERBOSE("EDGE %d\n", p_id);
		assert(g->edges[0].source == 0 && g->edges[0].target == 1);
	}
	kseq_destroy(seq);
	gzclose(fp);

	__VERBOSE("Number of nodes: %d\n", (int)g->n_v);
	__VERBOSE("Number of edges: %d\n", (int)g->n_e);
}
