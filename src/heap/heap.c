#include <stdio.h>
#include <stdlib.h>

#define __parent(i) i >> 1
#define __left(i) (i << 1)
#define __right(i) ((i << 1) + 1)

#define __swap_h(H, i, j, t) \
	t = H[i]; \
	H[i] = H[j]; \
	H[j] = t;

struct kheap_t {
    uint64_t *H;
    uint64_t *key;
    uint64_t *pos;
    uint32_t n;
    uint32_t s;
};

struct kheap_t *kheap_init(uint32_t n)
{
	struct kheap_t *h = malloc(sizeof(struct kheap_t));
	h->H = calloc(n, sizeof(uint64_t));
	h->key = calloc(n, sizeof(uint64_t));
	h->pos = calloc(n, sizeof(uint64_t));
	h->n = 0; //1-based
	h->s = n;
	return h;
}

void print_heap(struct kheap_t *h)
{
	for (uint64_t i = 1 ; i < h->n + 1; ++i) {
		printf("%llu ", h->key[h->H[i]]);
	}
	printf("\n");
}

int heapify_up(struct kheap_t *h, uint64_t i)
{
	uint64_t t, j;
	if (i == 1)
		return 0; // nothing to do
	j = __parent(i);
	if (h->key[h->H[i]] < h->key[h->H[j]]) {
		h->pos[h->H[i]] = j;
		h->pos[h->H[j]] = i;
		__swap_h(h->H, i, j, t)
		heapify_up(h, j);
	}
	return 1;
}

void kheap_insert(struct kheap_t *h, uint64_t v)
{
	if (h->n + 1 > h->s) { //1-based
		printf("Heap is too large! exit.\n");
		exit(1);
	}
	h->H[++h->n] = v;
	heapify_up(h, h->n); //1-based
}

int heapify_down(struct kheap_t *h, uint64_t i)
{
	if (__left(i) >= (h->n)) //1-based. Reached the leave, do nothing
		return 0;
	if (h->key[h->H[__left(i)]] < h->key[h->H[i]] || h->key[h->H[__right(i)]] < h->key[h->H[i]]) {
		uint64_t t;
		uint64_t j = (h->key[h->H[__left(i)]] < h->key[h->H[__right(i)]]) ? __left(i):__right(i);
		__swap_h(h->H, i, j, t);
		heapify_down(h, j);
	}
	return 1;
}

void kheap_delete(struct kheap_t *h, uint64_t v)
{
	if (h->n == 0) {
		printf("Trying to remove an element in an empty heap! exit.\n");
		exit(1);
	}
	h->H[v] = h->H[h->n--];
	heapify_down(h, v);
}

int main()
{
	struct kheap_t *h = kheap_init(16);
	uint32_t i;
	int e[14] = {1,2,5,10,3,7,11,15,17,20,9,15,8,16};
	for (i = 0; i < 14; ++i){
		h->key[i] = e[i];
	}
	for (i = 7; i < 14; ++i){
		kheap_insert(h, i);
	}
	for (i = 0; i < 1; ++i){
		kheap_insert(h, i);
	}
	print_heap(h);
	return 0;
}
