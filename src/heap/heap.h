//
// Created by BioTuring on 2019-11-24.
//

#ifndef SKIPPING_HEAP_H
#define SKIPPING_HEAP_H

#include <stdio.h>
#include <stdlib.h>

struct kheap_t {
	uint64_t *H;
	uint64_t *key;
	uint64_t *pos;
	uint32_t n;
	uint32_t s;
};


struct kheap_t *kheap_init(uint32_t n);
int kheap_delete(struct kheap_t *h, uint64_t v);
int kheap_insert(struct kheap_t *h, uint64_t v);

#endif //SKIPPING_HEAP_H
