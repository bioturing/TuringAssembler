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

struct mm_db_t {
    uint64_t *mm;
    uint32_t *p;
    size_t n;
    size_t size;
    int k;
};

struct mm_db_t * mm_index_str(uint32_t *s, int k, int w, int l);

#endif //SKIPPING_MINIMIZERS_H
