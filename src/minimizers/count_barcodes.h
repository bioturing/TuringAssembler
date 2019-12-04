//
// Created by BioTuring on 2019-09-23.
//

#ifndef SKIPPING_COUNT_BARCODES_H
#define SKIPPING_COUNT_BARCODES_H

#include <stdint.h>

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
uint64_t *mini_put_by_key(struct mini_hash_t *h_table, uint64_t data, uint64_t key);
uint64_t get_barcode_biot(char *s, struct read_t *r);
uint64_t *mini_get(struct mini_hash_t *h_table, uint64_t data);

#endif //SKIPPING_COUNT_BARCODES_H
