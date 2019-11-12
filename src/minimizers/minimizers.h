//
// Created by BioTuring on 2019-11-10.
//

#ifndef SKIPPING_MINIMIZERS_H
#define SKIPPING_MINIMIZERS_H

#include <stdint.h>
#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include "khash.h"
#include "assembly_graph.h"

KHASH_MAP_INIT_INT64(mm_hash, uint32_t);        /* Hash table structure of the minimizers */

struct mm_db_t {
    uint64_t *mm;
    uint32_t *p;
    size_t n;
    size_t size;
    int k;
};

struct mm_db_edge_t {
    kh_mm_hash_t *h;
    kh_mm_hash_t *cnt;
    kh_mm_hash_t *p;
};

struct mm_db_t * mm_index_bin_str(uint32_t *s, int k, int w, int l);
struct mm_db_t * mm_index_char_str(char *s, int k, int w, int l);
struct mm_db_edge_t *mm_index_edges(struct asm_graph_t *g, int k, int w);

#endif //SKIPPING_MINIMIZERS_H
