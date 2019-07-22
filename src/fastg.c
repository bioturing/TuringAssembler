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
        g->nodes[g->n_v].adj = malloc(sizeof(gint_t));
        g->nodes[g->n_v].adj[0] = id - 1;
        g->nodes[g->n_v].deg = 1;
        g->nodes[g->n_v + 1].adj = NULL;
        g->nodes[g->n_v + 1].deg = 0;

        g->edges[id].source = g->n_v + 2;
        g->edges[id].target = g->n_v + 3;
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

static int add_one_connection(struct asm_graph_t *g, gint_t e1, gint_t e2)
{
        if (g->edges[e1].target == g->edges[e2].source) {
                return 0;
        }
        gint_t e1_target = g->edges[e1].target;
        g->nodes[e1_target].adj = realloc(g->nodes[e1_target].adj, (g->nodes[e1_target].deg + 1) * sizeof(gint_t));
        g->nodes[e1_target].adj[g->nodes[e1_target].deg++] = e2;
        g->nodes[g->edges[e2].source].deg--;
        g->edges[e2].source = e1_target;
        g->edges[g->edges[e2].rc_id].target = g->nodes[e1_target].rc_id;
}

void load_asm_graph_fastg(struct asm_graph_t *g, const char *path, int ksize )
{
        g->ksize = ksize;
        g->aux_flag = 0;
        g->edges = NULL;
        g->nodes = NULL;
        g->n_v = g->n_e = 0;
        int c;
        char *p_i;
        int *added;
        gzFile fp = gzopen(path, "r");
        if (!fp)
                __ERROR("Unable to open file [%s] to read", path);
        kseq_t *seq = kseq_init(fp);
        while (kseq_read(seq) >= 0) {
                /* add new edge */
                if ((int)seq->seq.l < ksize)
                        continue;
                char *p, *s = seq->name.s;
                int is_comp;
                for (p = s; *p && *p != ':' && *p != ';'; ++p);
                c = *p, *p = 0;
                is_comp = (p > s && *(p-1) == '\''); // if we are looking at a complement segment
                if (is_comp) *(p-1) = 0;
                if (!is_comp) {
                        for(p_i = p; *p_i != '_'; --p_i);
                        add_one_edge(g, atoi(p_i + 1), seq, g->n_e > 0 && added[atoi(p_i + 1)]);
                        //printf("S\t%s\t%s\tLN:i:%ld\n", s, seq->seq.s, (long)strlen(seq->seq.s));
                }
                if (atoi(p_i + 1) + 1 == g->n_e)
                        added = realloc(added, sizeof(int) * g->n_e);
                added[atoi(p_i + 1)] = 1;
                added[atoi(p_i + 1) - 1] = 1;
                if (c == ':') { // have neighbors
                        char *q = p + 1;
                        do {
                                int is_comp2 = 0;
                                for (p = q; *p && *p != ',' && *p != ';'; ++p);
                                c = *p, *p = 0;
                                is_comp2 = (p > q && *(p-1) == '\'');
                                if (is_comp2) *(p-1) = 0;
                                char *q_i;
                                for(q_i = p; *q_i != '_'; --q_i);
                                gint_t q_id = atoi(q_i + 1);
                                gint_t p_id = atoi(p_i + 1);
                                if (q_id > p_id) {
                                        if (q_id > g->n_e || !added[q_id]) {
                                                add_one_edge(g, q_id , NULL, 0);
                                                added = realloc(added, sizeof(int) * g->n_e);
                                        }
                                        added[q_id] = 1;
                                        added[q_id - 1] = 1;
                                }
                                add_one_connection(g, p_id - 1 + is_comp, q_id - 1 + is_comp2);
                                //printf("L\t%s\t%c\t%s\t%c\t0M\n", s, "+-"[!!is_comp], q, "+-"[!!is_comp2]);
                                q = p + 1;
                        } while (c != 0 && c != ';');
                }
                //__VERBOSE("ID %d\n", atoi(p_i + 1));
                //assert(g->n_v > g->n_e);
        }
        kseq_destroy(seq);
        gzclose(fp);

        __VERBOSE("Number of nodes: %d\n", g->n_v);
        __VERBOSE("Number of edges: %d\n", g->n_e);
}
