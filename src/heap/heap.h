//
// Created by BioTuring on 2019-11-24.
//

#ifndef SKIPPING_HEAP_H
#define SKIPPING_HEAP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>


struct kheap_t {
	uint32_t *H;
	uint32_t *key;
	uint32_t *pos;
        uint32_t *v_cnt;
	uint32_t n;
	uint32_t s;
};


struct kheap_t *kheap_init(uint32_t n);
int kheap_delete(struct kheap_t *h, uint32_t v);
int kheap_insert(struct kheap_t *h, uint32_t v);
int heapify_up(struct kheap_t *h, uint32_t i);

#endif //SKIPPING_HEAP_H
