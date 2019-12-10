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

#define RATIO_OF_CONFIDENT 0.85
#define MIN_NUMBER_SINGLETON 2
#define EMPTY_BX UINT64_MAX

KHASH_MAP_INIT_INT64(mm_hash, uint32_t);        /* Hash table structure of the minimizers */
KHASH_MAP_INIT_INT64(mm_edges, uint32_t);                 /* Hash table edge count */
KHASH_MAP_INIT_INT64(mm_pw, uint32_t);                 /* Hash table edge count */

struct mm_align_t {
    uint64_t edge;
    uint32_t pos;
    uint32_t cnt;
    uint64_t mm;
};

KHASH_MAP_INIT_INT64(mm_align, struct mm_align_t *);

struct mm_db_t {
    uint64_t *mm;
    uint32_t *p;
    size_t n;
    size_t size;
    int k;
};

struct mm_hits_t {
    kh_mm_edges_t *edges;
    kh_mm_align_t *aln;
    uint32_t n;
};

struct mm_db_edge_t {
    kh_mm_hash_t *h;
    kh_mm_hash_t *cnt;
    kh_mm_hash_t *p;
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
