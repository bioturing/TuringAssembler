#include <stdlib.h>

#include "atomic.h"
#include "fastq_producer.h"
#include "get_buffer.h"
#include "utils.h"
#include "verbose.h"

void *init_single_buffer()
{
	struct single_buffer_t *ret = malloc(sizeof(struct single_buffer_t));
	ret->R_buf = malloc(BUF_SIZE + 1);
	return ret;
}

void *init_pair_buffer()
{
	struct pair_buffer_t *ret = malloc(sizeof(struct pair_buffer_t));
	ret->R1_buf = malloc(BUF_SIZE + 1);
	ret->R2_buf = malloc(BUF_SIZE + 1);
	return ret;
}

void *init_trip_buffer()
{
	struct trip_buffer_t *ret = malloc(sizeof(struct trip_buffer_t));
	ret->R1_buf = malloc(BUF_SIZE + 1);
	ret->R2_buf = malloc(BUF_SIZE + 1);
	ret->I_buf = malloc(BUF_SIZE + 1);
	return ret;
}

void free_single_buffer(void *vp)
{
	if (!vp) return;
	struct single_buffer_t *p = (struct single_buffer_t *)vp;
	free(p->R_buf);
	free(p);
}

void free_pair_buffer(void *vp)
{
	if (!vp) return;
	struct pair_buffer_t *p = (struct pair_buffer_t *)vp;
	free(p->R1_buf);
	free(p->R2_buf);
	free(p);
}

void free_trip_buffer(void *vp)
{
	if (!vp) return;
	struct trip_buffer_t *p = (struct trip_buffer_t *)vp;
	free(p->R1_buf);
	free(p->R2_buf);
	free(p->I_buf);
	free(p);
}

static struct dqueue_t *init_dqueue_fastq(int cap, void *(*buffer_init)())
{
	struct dqueue_t *ret = init_dqueue(cap);
	struct pair_buffer_t *p;
	int i;
	for (i = 0; i < cap; ++i) {
		p = buffer_init();
		d_enqueue_out(ret, p);
	}
	return ret;
}

void *fastq_producer(void *data)
{
	struct producer_bundle_t *bundle = (struct producer_bundle_t *)data;
	void *(*buffer_init)() = bundle->buffer_init;
	void (*buffer_free)(void *) = bundle->buffer_free;
	int64_t (*gb_get_data)(void *, void *) = bundle->gb_get_data;
	struct dqueue_t *q = bundle->q;
	void *input_stream = bundle->stream;
	void *own_buf = buffer_init();
	void *external_buf;
	int64_t total_size = bundle->total_size;

	int64_t prev_processed = 0, cur_processed, global_processed, percentage;
	percentage = 0;
	while ((cur_processed = gb_get_data(input_stream, own_buf)) >= 0) {
		external_buf = d_dequeue_out(q);
		d_enqueue_in(q, own_buf);
		own_buf = external_buf;
		global_processed = atomic_add_and_fetch64(bundle->processed_size,
						cur_processed - prev_processed);
		prev_processed = cur_processed;
		percentage = global_processed * 100 / total_size;
		percentage = __min(percentage, 99);
		__VERBOSE("\rLoad %ld%%", percentage);
	}
	buffer_free(own_buf);
	__VERBOSE("\rLoad 100%%");

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
		buffer_free(external_buf);
		d_enqueue_in(q, NULL);
	}

	pthread_exit(NULL);
}

struct producer_bundle_t *init_fastq_single(int n_threads, int n_files,
								char **files)
{
	struct producer_bundle_t *ret;
	struct dqueue_t *q;
	int i;
	int *n_consumer;
	pthread_mutex_t *lock;
	pthread_barrier_t *barrier;

	ret = calloc(n_files, sizeof(struct producer_bundle_t));
	q = init_dqueue_fastq(__max(n_threads, n_files) * 2, init_single_buffer);
	n_consumer = malloc(sizeof(int));
	*n_consumer = n_threads * 2;
	lock = malloc(sizeof(pthread_mutex_t));
	barrier = malloc(sizeof(pthread_barrier_t));

	pthread_mutex_init(lock, NULL);
	pthread_barrier_init(barrier, NULL, n_files);

	int64_t total_size = 0;
	int64_t *processed_size = malloc(sizeof(int64_t));
	*processed_size = 0;

	for (i = 0; i < n_files; ++i) {
		struct gb_single_data *data = calloc(1, sizeof(struct gb_single_data));
		gb_single_init(data, files[i]);
		total_size += data->compressed_size;
		ret[i].processed_size = processed_size;

		ret[i].n_consumer = n_consumer;
		ret[i].q = q;
		ret[i].barrier = barrier;
		ret[i].lock = lock;
		ret[i].stream = (void *)data;
		ret[i].buffer_init = init_single_buffer;
		ret[i].buffer_free = free_single_buffer;
		ret[i].gb_get_data = gb_get_single;
	}

	for (i = 0; i < n_files; ++i)
		ret[i].total_size = total_size;
	return ret;
}

struct producer_bundle_t *init_fastq_pair(int n_threads, int n_files,
					char **R1_files, char **R2_files)
{
	struct producer_bundle_t *ret;
	struct dqueue_t *q;
	int i;
	int *n_consumer;
	pthread_mutex_t *lock;
	pthread_barrier_t *barrier;

	ret = calloc(n_files, sizeof(struct producer_bundle_t));
	q = init_dqueue_fastq(__max(n_threads, n_files) * 2, init_pair_buffer);
	n_consumer = malloc(sizeof(int));
	*n_consumer = n_threads * 2;
	lock = malloc(sizeof(pthread_mutex_t));
	barrier = malloc(sizeof(pthread_barrier_t));

	pthread_mutex_init(lock, NULL);
	pthread_barrier_init(barrier, NULL, n_files);

	int64_t total_size = 0;
	int64_t *processed_size = malloc(sizeof(int64_t));
	*processed_size = 0;

	for (i = 0; i < n_files; ++i) {
		struct gb_pair_data *data = calloc(1, sizeof(struct gb_pair_data));
		gb_pair_init(data, R1_files[i], R2_files[i]);
		total_size += data->compressed_size;
		ret[i].processed_size = processed_size;

		ret[i].n_consumer = n_consumer;
		ret[i].q = q;
		ret[i].barrier = barrier;
		ret[i].lock = lock;
		ret[i].stream = (void *)data;
		ret[i].buffer_init = init_pair_buffer;
		ret[i].buffer_free = free_pair_buffer;
		ret[i].gb_get_data = gb_get_pair;
	}

	for (i = 0; i < n_files; ++i)
		ret[i].total_size = total_size;
	return ret;
}

struct producer_bundle_t *init_fastq_triple(int n_threads, int n_files,
			char **R1_files, char **R2_files, char **I_files)
{
	struct producer_bundle_t *ret;
	struct dqueue_t *q;
	int i;
	int *n_consumer;
	pthread_mutex_t *lock;
	pthread_barrier_t *barrier;

	ret = calloc(n_files, sizeof(struct producer_bundle_t));
	q = init_dqueue_fastq(__max(n_threads, n_files) * 2, init_trip_buffer);
	n_consumer = malloc(sizeof(int));
	*n_consumer = n_threads * 2;
	lock = malloc(sizeof(pthread_mutex_t));
	barrier = malloc(sizeof(pthread_barrier_t));

	pthread_mutex_init(lock, NULL);
	pthread_barrier_init(barrier, NULL, n_files);

	int64_t total_size = 0;
	int64_t *processed_size = malloc(sizeof(int64_t));
	*processed_size = 0;

	for (i = 0; i < n_files; ++i) {
		struct gb_trip_data *data = calloc(1, sizeof(struct gb_trip_data));
		gb_trip_init(data, R1_files[i], R2_files[i], I_files[i]);
		total_size += data->compressed_size;
		ret[i].processed_size = processed_size;

		ret[i].n_consumer = n_consumer;
		ret[i].q = q;
		ret[i].barrier = barrier;
		ret[i].lock = lock;
		ret[i].stream = (void *)data;
		ret[i].buffer_init = init_trip_buffer;
		ret[i].buffer_free = free_trip_buffer;
		ret[i].gb_get_data = gb_get_trip;
	}

	for (i = 0; i < n_files; ++i)
		ret[i].total_size = total_size;
	return ret;
}

void free_fastq_single(struct producer_bundle_t *bundles, int n)
{
	int i;
	for (i = 0; i < n; ++i) {
		gb_single_destroy((struct gb_single_data *)bundles[i].stream);
		free((struct gb_single_data *)bundles[i].stream);
	}
	pthread_mutex_destroy(bundles->lock);
	pthread_barrier_destroy(bundles->barrier);
	dqueue_destroy(bundles->q);
	free(bundles->n_consumer);
	free(bundles->lock);
	free(bundles->barrier);
	free(bundles);
}

void free_fastq_pair(struct producer_bundle_t *bundles, int n)
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

void free_fastq_triple(struct producer_bundle_t *bundles, int n)
{
	int i;
	for (i = 0; i < n; ++i) {
		gb_trip_destroy((struct gb_trip_data *)bundles[i].stream);
		free((struct gb_trip_data *)bundles[i].stream);
	}
	pthread_mutex_destroy(bundles->lock);
	pthread_barrier_destroy(bundles->barrier);
	dqueue_destroy(bundles->q);
	free(bundles->n_consumer);
	free(bundles->lock);
	free(bundles->barrier);
	free(bundles);
}
