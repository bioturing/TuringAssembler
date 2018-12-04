#include "fastq_producer.h"
#include "utils.h"

struct pair_buffer_t *init_pair_buffer()
{
	struct pair_buffer_t *ret = malloc(sizeof(struct pair_buffer_t));
	ret->buf1 = malloc(BUF_SIZE + 1);
	ret->buf2 = malloc(BUF_SIZE + 1);
	return ret;
}

void free_pair_buffer(struct pair_buffer_t *p)
{
	if (!p) return;
	free(p->buf1);
	free(p->buf2);
	free(p);
}

static struct dqueue_t *init_dqueue_PE(int cap)
{
	struct dqueue_t *ret = init_dqueue(cap);
	struct pair_buffer_t *p;
	int i;
	for (i = 0; i < cap; ++i) {
		p = init_pair_buffer();
		d_enqueue_out(ret, p);
	}
	return ret;
}

void *fastq_PE_producer(void *data)
{
	struct producer_bundle_t *bundle = (struct producer_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct gb_pair_data *input_stream = bundle->stream;
	struct pair_buffer_t *own_buf = init_pair_buffer();
	struct pair_buffer_t *external_buf;
	int64_t offset;
	int64_t n_frag = 0;

	while ((offset = gb_get_pair(input_stream, &own_buf->buf1, &own_buf->buf2)) != -1) {
		own_buf->input_format = input_stream->type;
		external_buf = d_dequeue_out(q);
		d_enqueue_in(q, own_buf);
		own_buf = external_buf;
		++n_frag;
	}
	free_pair_buffer(own_buf);

	int cur;
	pthread_barrier_wait(bundle->barrier);
	while (1) {
		pthread_mutex_lock(bundle->lock);
		cur = *(bundle->n_consumer);
		if (*(bundle->n_consumer) > 0)
			--*(bundle->n_consumer);
		pthread_mutex_unlock(bundle->lock);
		if (cur == 0)
			break;
		external_buf = d_dequeue_out(q);
		free_pair_buffer(external_buf);
		d_enqueue_in(q, NULL);
	}

	pthread_exit(NULL);
}

struct producer_bundle_t *init_fastq_PE(struct opt_count_t *opt)
{
	struct producer_bundle_t *ret;
	struct dqueue_t *q;
	int i;
	int *n_consumer;
	pthread_mutex_t *lock;
	pthread_barrier_t *barrier;

	ret = calloc(opt->n_files, sizeof(struct producer_bundle_t));
	q = init_dqueue_PE(__max(opt->n_threads, opt->n_files) * 2);
	n_consumer = malloc(sizeof(int));
	*n_consumer = opt->n_threads * 2;
	lock = malloc(sizeof(pthread_mutex_t));
	barrier = malloc(sizeof(pthread_barrier_t));

	for (i = 0; i < opt->n_files; ++i) {
		struct gb_pair_data *data = calloc(1, sizeof(struct gb_pair_data));
		gb_pair_init(data, opt->files_1[i], opt->files_2[i]);

		ret[i].n_consumer = n_consumer;
		ret[i].stream = (void *)data;
		ret[i].q = q;
		ret[i].barrier = barrier;
		ret[i].lock = lock;
	}
	return ret;
}

void free_fastq_PE(struct producer_bundle_t *bundles, int n)
{
	int i;
	for (i = 0; i < n; ++i) {
		gb_pair_destroy((struct gb_pair_data *)bundles[i].stream);
		free((struct gb_pair_data *)bundles[i].stream);
	}
	pthread_mutex_destroy(bundles->lock);
	pthread_barrier_destroy(bundles->barrier);
	dqueue_destroy(bundles->q);
	free(bundles->n_consumer);
	free(bundles->lock);
	free(bundles->barrier);
	free(bundles);
}
