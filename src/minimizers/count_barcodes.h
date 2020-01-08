//
// Created by BioTuring on 2019-09-23.
//

#ifndef SKIPPING_COUNT_BARCODES_H
#define SKIPPING_COUNT_BARCODES_H

#include <stdint.h>
#include "../cluster_molecules.h"

#define EMPTY_SLOT 0xffffffffffffffff
#define __mini_empty(table, i) (table->key[i] == EMPTY_SLOT)
#define INIT_PRIME_INDEX 16

struct mini_hash_t {
    uint64_t *h;
    uint64_t *key;
    uint64_t size;
    uint64_t count;
    uint64_t max_cnt;
    int prime_index;
};

struct readsort_bundle1_t {
    struct dqueue_t *q;
    char prefix[MAX_PATH];
    int64_t sm;
    struct mini_hash_t **h_table_ptr;
};

struct mini_hash_t *count_bx_freq(struct opt_proc_t *opt);
uint64_t barcode_hash_mini(char *s);
void init_mini_hash(struct mini_hash_t **h_table, uint32_t p_index);
void destroy_mini_hash(struct mini_hash_t *h_table);
void destroy_worker_bundles(struct readsort_bundle1_t *bundles, int n);
uint64_t *mini_put(struct mini_hash_t **h_table, uint64_t data);
uint64_t *mini_put_value(struct mini_hash_t **h_table, uint64_t data, int value);
uint64_t *mini_put_by_key(struct mini_hash_t *h_table, uint64_t data, uint64_t key);
uint64_t get_barcode_biot(char *s, struct read_t *r);
uint64_t *mini_get(struct mini_hash_t *h_table, uint64_t data);
khash_t(long_int) *count_edge_link_shared_bc(struct asm_graph_t *g,
		struct mini_hash_t *bc);
uint64_t get_min_code(struct asm_graph_t *g, int u, int v);
void destroy_bx_table(struct mini_hash_t *bx_table);
int get_edge_link_bc_count(khash_t(long_int) *all_count, int e1, int e2);

#endif //SKIPPING_COUNT_BARCODES_H
