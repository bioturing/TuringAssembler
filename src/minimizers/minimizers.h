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
#include "../khash.h"
#include "../assembly_graph.h"
#include "minimizers/count_barcodes.h"

#define RATIO_OF_CONFIDENT 0.85
#define MIN_NUMBER_SINGLETON 2
#define EMPTY_BX UINT64_MAX

struct mm_align_t {
    uint64_t edge;
    uint32_t pos;
    uint32_t cnt;
    uint64_t mm;
};


struct mm_db_t {
    uint64_t *mm;
    uint32_t *p;
    size_t n;
    size_t size;
    int k;
};

struct mm_hits_t {
    struct mini_hash_t *edges;
    uint32_t n;
};

struct mm_db_edge_t {
    struct mini_hash_t *table;
    uint64_t cnt;
};

struct mm_bundle_t {
    struct mini_hash_t *bx_table; /* hash table for share barcode count */
    struct mini_hash_t *rp_table; /* hash table for read pair count */
};

struct mm_hits_t *mm_hits_init();
void mm_hits_print(struct mm_hits_t *hits, const char *file_path);
struct mm_db_t * mm_index_bin_str(uint32_t *s, int k, int w, int l);
struct mm_db_t * mm_index_char_str(char *s, int k, int w, int l);
struct mm_db_edge_t *mm_index_edges(struct asm_graph_t *g, int k, int w);
void *mm_hits_cmp(struct mm_db_t *db, struct mm_db_edge_t *db_e, struct mm_hits_t *hits, struct asm_graph_t *g);
void mm_db_destroy(struct mm_db_t *db);
void mm_hits_destroy(struct mm_hits_t *hits);
struct mm_bundle_t *mm_hit_all_barcodes(struct opt_proc_t *opt, struct asm_graph_t *g);

#endif //SKIPPING_MINIMIZERS_H
