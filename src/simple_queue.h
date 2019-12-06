#ifndef __SIMPLE_QUEUE__
#define __SIMPLE_QUEUE__
struct queue_t{
	void **data;
	int cap;
	int front;
	int back;
};

void init_queue(struct queue_t *q, int cap);
void push_queue(struct queue_t *q, void *ptr);
void *get_queue(struct queue_t *q);
void pop_queue(struct queue_t *q);
int is_queue_empty(struct queue_t *q);
void destroy_queue(struct queue_t *q);
void *pointerize(void *data, int size);
void free_queue_content(struct queue_t *q);

struct heap_t{
	struct queue_t *q;
	int (*cmp)(void *, void *);
};

void init_heap(struct heap_t *heap, int (*cmp)(void *, void *));
void up_heap(struct heap_t *heap);
void down_heap(struct heap_t *heap);
void *get_heap(struct heap_t *heap);
void push_heap(struct heap_t *heap, void *data);
void pop_heap(struct heap_t *heap);
int is_heap_empty(struct heap_t *heap);
void heap_destroy(struct heap_t *heap);
#endif
