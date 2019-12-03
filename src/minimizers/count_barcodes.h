//
// Created by BioTuring on 2019-09-23.
//

#ifndef SKIPPING_COUNT_BARCODES_H
#define SKIPPING_COUNT_BARCODES_H

#include <stdint.h>

struct mini_hash_t {
    uint32_t *h;
    uint64_t *key;
    uint64_t size;
    uint64_t count;
    uint64_t max_cnt;
    int prime_index;
};

void count_bx_freq(struct opt_proc_t *opt, struct read_path_t *r_path);
uint64_t barcode_hash_mini(char *s);
void init_mini_hash(struct mini_hash_t **h_table, uint32_t p_index);
void destroy_mini_hash(struct mini_hash_t *h_table);
int mini_inc_by_key(struct mini_hash_t *h_table, uint64_t data, uint64_t key);
uint64_t mini_get(struct mini_hash_t *h_table, uint64_t data, uint64_t key);

#endif //SKIPPING_COUNT_BARCODES_H
