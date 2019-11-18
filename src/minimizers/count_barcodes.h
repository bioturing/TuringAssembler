//
// Created by BioTuring on 2019-09-23.
//

#ifndef SKIPPING_COUNT_BARCODES_H
#define SKIPPING_COUNT_BARCODES_H


#include <stdint.h>

void count_bx_freq(struct opt_proc_t *opt, struct read_path_t *r_path);
uint64_t barcode_hash_mini(char *s);
struct mini_hash_t *init_mini_hash();
void mini_inc_by_key(uint64_t data, uint64_t key);
uint64_t mini_get(uint64_t data, uint64_t key);

#endif //SKIPPING_COUNT_BARCODES_H
