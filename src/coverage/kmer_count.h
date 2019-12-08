//
// Created by BioTuring on 2019-12-08.
//

#ifndef SKIPPING_KMER_COUNT_H
#define SKIPPING_KMER_COUNT_H

struct mini_hash_t * construct_edges_hash(struct asm_graph_t *g);
struct mini_hash_t *kmer_count_on_edges(struct opt_proc_t *opt);
void add_cnt_to_graph(struct asm_graph_t *g, struct mini_hash_t *kmer_table);
#endif //SKIPPING_KMER_COUNT_H
