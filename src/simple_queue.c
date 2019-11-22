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

