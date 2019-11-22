#include <stdlib.h>
#include <string.h>
#include "simple_queue.h"

void *pointerize(void *data, int size)
{
	void *res = malloc(size);
	memcpy(res, data, size);
	return res;
}

void init_queue(struct queue_t *q, int cap)
{
	q->data = calloc(cap, sizeof(void *));
	q->cap = cap;
	q->front = 0;
	q->back = 0;
}

void push_queue(struct queue_t *q, void *ptr)
{
	if (q->back == q->cap){
		q->cap = q->cap == 0 ? 1 : (q->cap << 1);
		q->data = realloc(q->data, sizeof(void *) * q->cap);
	}
	q->data[q->back++] = ptr;
}

void *get_queue(struct queue_t *q)
{
	return q->data[q->front];
}

void pop_queue(struct queue_t *fq)
{
	fq->front++;
}

int is_queue_empty(struct queue_t *q)
{
	return q->front >= q->back;
}

void destroy_queue(struct queue_t *q)
{
	free(q->data);
}

void free_queue_content(struct queue_t *q)
{
	for (int i = q->front; i < q->back; ++i)
		free(q->data[i]);
}

void init_heap(struct heap_t *heap, int (*cmp)(void *, void *))
{
	heap->q = calloc(1, sizeof(struct queue_t));
	init_queue(heap->q, 1024);
	heap->cmp = cmp;
}

void up_heap(struct heap_t *heap)
{
	struct queue_t *q = heap->q;
	int (*cmp)(void *, void *) = heap->cmp;
	for (int i = q->back - 1; i / 2 >= q->front; i = i / 2){
		if (cmp(q->data[i], q->data[i / 2]) == -1){
			void *tmp = q->data[i];
			q->data[i] = q->data[i / 2];
			q->data[i / 2] = tmp;
		} else {
			break;
		}
	}
}

void down_heap(struct heap_t *heap)
{
	struct queue_t *q = heap->q;
	if (is_queue_empty(q))
		return;
	int (*cmp)(void *, void *) = heap->cmp;
	for (int i = q->front; i * 2 < q->back; ){
		int j = i * 2;
		if (j + 1 < q->back && cmp(q->data[j + 1], q->data[j]) == -1)
			++j;
		if (cmp(q->data[j], q->data[i]) == -1){
			void *tmp = q->data[i];
			q->data[i] = q->data[j];
			q->data[j] = tmp;
			i = j;
		} else {
			break;
		}
	}
}

void push_heap(struct heap_t *heap, void *data)
{
	push_queue(heap->q, data);
	up_heap(heap);
}

void *get_heap(struct heap_t *heap)
{
	return get_queue(heap->q);
}

void pop_heap(struct heap_t *heap)
{
	pop_queue(heap->q);
	down_heap(heap);
}

int is_heap_empty(struct heap_t *heap)
{
	return is_queue_empty(heap->q);
}

void heap_destroy(struct heap_t *heap)
{
	free_queue_content(heap->q);
	free(heap->q);
}

