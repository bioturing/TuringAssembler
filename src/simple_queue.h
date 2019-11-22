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

#endif
