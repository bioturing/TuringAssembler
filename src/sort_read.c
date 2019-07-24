#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "attribute.h"
#include "fastq_producer.h"
#include "io_utils.h"
#include "radix_sort.h"
#include "sort_read.h"
#include "time_utils.h"
#include "utils.h"
#include "verbose.h"

struct readbc_t {
	uint64_t barcode;
	int64_t offset;
	int len1;
	int len2;
};

#define bc_get_block(p, s, mask) ((p).barcode >> (s) & (mask))
#define bc_less_than(x, y) ((x).barcode < (y).barcode)

RS_IMPL(read_sort, struct readbc_t, 64, 8, bc_less_than, bc_get_block)

struct readsort_bundle_t {
	struct dqueue_t *q;
	char prefix[MAX_PATH];
	int64_t sm;
};

static inline uint32_t unpack_int32(uint8_t *buf)
{
	uint32_t ret = 0;
	ret =	(uint32_t)buf[0]		|
		((uint32_t)buf[1] << 8)	|
		((uint32_t)buf[2] << 16)	|
		((uint32_t)buf[3] << 24);
	return ret;
}

static inline uint64_t unpack_int64(uint8_t *buf)
{
	uint64_t ret = 0;
	ret =	(uint64_t)buf[0]		|
		((uint64_t)buf[1] << 8)	|
		((uint64_t)buf[2] << 16)	|
		((uint64_t)buf[3] << 24)	|

		((uint64_t)buf[4] << 32)	|
		((uint64_t)buf[5] << 40)	|
		((uint64_t)buf[6] << 48)	|
		((uint64_t)buf[7] << 56);
	return ret;
}

static inline void pack_int32(uint8_t *buffer, uint32_t value)
{
	buffer[0] = value;
	buffer[1] = value >> 8;
	buffer[2] = value >> 16;
	buffer[3] = value >> 24;
}

static inline void pack_int64(uint8_t *buffer, uint64_t value)
{
	buffer[0] = value;
	buffer[1] = value >> 8;
	buffer[2] = value >> 16;
	buffer[3] = value >> 24;

	buffer[4] = value >> 32;
	buffer[5] = value >> 40;
	buffer[6] = value >> 48;
	buffer[7] = value >> 56;
}

static inline uint64_t get_barcode_ust_raw(struct read_t *I)
{
	uint64_t ret = 0;
	int i;
	for (i = 0; i < I->len; ++i)
		ret = ret * 5 + nt4_table[(int)I->seq[i]];
	return ret;
}

static inline void write_buffer_to_file(struct readbc_t *p, int n, char *buf,
		int64_t buf_len, char *fbuf, int64_t fbuf_size, const char *path)
{
	/* sort p */
	rs_sort(read_sort, p, p + n);
	/*******/
	FILE *fp = xfopen(path, "wb");
	setbuffer(fp, fbuf, fbuf_size);
	int i;
	for (i = 0; i < n; ++i) {
		fwrite(buf + p->offset, p->len1 + p->len2 + 16, 1, fp);
	}
	fclose(fp);
}

static inline int ust_add_record(struct read_t *r, struct read_t *rI,
						char *buf, int input_format)
{
	int len = 0, tlen, i;
	tlen = strlen(r->name);
	memcpy(buf, r->name, tlen);
	len = tlen;
	buf[len++] = '\t';
	memcpy(buf + len, "BX:Z:", 5);
	len += 5;
	memcpy(buf + len, rI->seq, rI->len);
	len += rI->len;
	buf[len++] = '\t';
	memcpy(buf + len, "QB:Z:", 5);
	len += 5;
	if (input_format == TYPE_FASTQ) {
		memcpy(buf + len, rI->qual, rI->len);
		len += rI->len;
	} else {
		for (i = 0; i < rI->len; ++i) {
			buf[len++] = nt4_table[(int)rI->seq[i]] < 4 ? 'F' : '#';
		}
	}
	buf[len++] = '\n';
	memcpy(buf + len, r->seq, r->len);
	len += r->len;
	buf[len++] = '\n';
	memcpy(buf + len, "+\n", 2);
	len += 2;
	if (input_format == TYPE_FASTQ) {
		memcpy(buf + len, r->qual, r->len);
		len += r->len;
	} else {
		for (i = 0; i < r->len; ++i) {
			buf[len++] = nt4_table[(int)r->seq[i]] < 4 ? 'F' : '#';
		}
	}
	buf[len++] = '\n';
	return len;
}

static inline void extract_read_barcode(FILE *fp, struct readbc_t *r, char *buf)
{
	size_t ret = fread(buf, 1, 16, fp);
	if (ret == 0) {
		r->barcode = (uint64_t)-1;
	} else if (ret < 16) {
		__ERROR("Corrupted temporary files");
	} else {
		r->barcode = unpack_int64((uint8_t *)buf);
		r->len1 = unpack_int32((uint8_t *)buf + 8);
		r->len2 = unpack_int32((uint8_t *)buf + 12);
	}
}

void merge_sorted_small(const char *prefix, int64_t sm, int n_file)
{
	int64_t sm_per_file;
	char **fbuf, *tmp_buf, *obuf, *buf, *path;
	FILE **fp, *fo;
	struct readbc_t *reads;
	int i, m_buf;
	sm_per_file = sm / n_file;
	__round_up_64(sm_per_file);
	if (sm_per_file > sm / n_file)
		sm_per_file >>= 1;
	fp = malloc(n_file * sizeof(FILE *));
	fbuf = malloc(n_file * sizeof(char *));
	reads = malloc(n_file * sizeof(struct readbc_t));
	tmp_buf = malloc(16);
	path = malloc(MAX_PATH);
	for (i = 0; i < n_file; ++i) {
		sprintf(path, "%s.tmp.%d", prefix, i);
		fp[i] = xfopen(path, "rb");
		fbuf[i] = malloc(sm_per_file);
		setbuffer(fp[i], fbuf[i], sm_per_file);
		extract_read_barcode(fp[i], reads + i, tmp_buf);
	}
	sprintf(path, "%s.tmp", prefix);
	fo = xfopen(path, "wb");
	obuf = malloc(SIZE_2MB);
	setbuffer(fo, obuf, SIZE_2MB);
	m_buf = 0x100;
	buf = malloc(m_buf);
	while (1) {
		uint64_t cur_bc = (uint64_t)-1;
		int idx = -1;
		for (i = 0; i < n_file; ++i) {
			if (reads[i].barcode != (uint64_t)-1 &&
				(cur_bc == (uint64_t)-1 || reads[i].barcode < cur_bc)) {
				cur_bc = reads[i].barcode;
				idx = i;
			}
		}
		if (cur_bc == (uint64_t)-1)
			break;
		if (reads[idx].len1 + reads[idx].len2 > m_buf) {
			m_buf = reads[idx].len1 + reads[idx].len2;
			buf = realloc(buf, m_buf);
		}
		pack_int64((uint8_t *)tmp_buf, reads[idx].barcode);
		pack_int32((uint8_t *)tmp_buf + 8, reads[idx].len1);
		pack_int32((uint8_t *)tmp_buf + 12, reads[idx].len2);
		fwrite(tmp_buf, 16, 1, fo);
		size_t byte_read = fread(buf, 1, reads[idx].len1 + reads[idx].len2, fp[idx]);
		if ((int64_t)byte_read < reads[idx].len1 + reads[idx].len2)
			__ERROR("Corrupted temporary files");
		fwrite(buf, reads[idx].len1 + reads[idx].len2, 1, fo);
		extract_read_barcode(fp[idx], reads + idx, tmp_buf);
	}
	for (i = 0; i < n_file; ++i) {
		sprintf(path, "%s.tmp.%d", prefix, i);
		fclose(fp[i]);
		free(fbuf[i]);
		remove(path);
	}
	free(buf);
	free(obuf);
	free(path);
	free(tmp_buf);
	free(reads);
	free(fbuf);
	free(fp);
	fclose(fo);
}

void *biot_buffer_iterator(void *data)
{
	return NULL;
}

void *x10_buffer_iterator(void *data)
{
	return NULL;
}

void *ust_buffer_iterator(void *data)
{
	struct readsort_bundle_t *bundle = (struct readsort_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2, readI;
	struct trip_buffer_t *own_buf, *ext_buf;
	own_buf = init_trip_buffer();
	char *prefix = bundle->prefix;
	int64_t sm = bundle->sm;
	char *path = malloc(MAX_PATH);
	fprintf(stderr, "thread info: prefix = %s; sm = %ld\n", prefix, sm);

	char *buf, *fbuf;
	int64_t buf_len, buf_size, fbuf_size;
	buf = malloc(sm);
	buf_len = 0;
	buf_size = sm;
	fbuf = malloc(SIZE_2MB);
	fbuf_size = SIZE_2MB;

	int m_read, n_read, n_file;
	struct readbc_t *rbc;
	m_read = sm / 600;
	__round_up_32(m_read);
	rbc = malloc(m_read * sizeof(struct readbc_t));
	n_read = 0;
	n_file = 0;

	char *R1_buf, *R2_buf, *I_buf;
	int pos1, pos2, posI, rc1, rc2, rcI, input_format;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf)
			break;
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = posI = 0;
		R1_buf = ext_buf->R1_buf;
		R2_buf = ext_buf->R2_buf;
		I_buf = ext_buf->I_buf;
		input_format = ext_buf->input_format;

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, R1_buf, &pos1) :
				get_read_from_fa(&read1, R1_buf, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, R2_buf, &pos2) :
				get_read_from_fa(&read2, R2_buf, &pos2);

			rcI = input_format == TYPE_FASTQ ?
				get_read_from_fq(&readI, I_buf, &posI) :
				get_read_from_fa(&readI, I_buf, &posI);

			if (rc1 == READ_FAIL || rc2 == READ_FAIL || rcI == READ_FAIL)
				__ERROR("\nWrong format file\n");

			/* read_name + \t + BX:Z: + barcode + \t + QB:Z: + barcode_quality + \n */
			uint64_t barcode = get_barcode_ust_raw(&readI);
			int record_len, len1, len2;
			record_len = strlen(read1.name) + strlen(read2.name) +
				2 * (1 + (readI.len + 6) * 2) +
				(read1.len + 1) * 2 + (read2.len + 1) * 2 + 2 * 2 +
				sizeof(barcode) + sizeof(len1) * 2;
			if (record_len > buf_size) {
				/* extend buf size */
				buf_size = record_len;
				buf = realloc(buf, buf_size);
			}
			if (record_len > buf_size - buf_len) {
				/* sort and flush buffer */;
				sprintf(path, "%s.tmp.%d", prefix, n_file);
				write_buffer_to_file(rbc, n_read, buf, buf_len,
							fbuf, fbuf_size, path);
				++n_file;
				n_read = 0;
				buf_len = 0;
			}
			/* append record */
			pack_int64((uint8_t *)buf + buf_len, barcode);
			len1 = ust_add_record(&read1, &readI,
					buf + buf_len + 16, input_format);
			len2 = ust_add_record(&read2, &readI,
				buf + buf_len + 16 + len1, input_format);
			if (record_len != len1 + len2 + 16) {
				fprintf(stderr, "record_len = %d; sum_len = %d\n",
					record_len, len1 + len2 + 16);
				assert(0);
			}
			pack_int32((uint8_t *)buf + buf_len + 8, len1);
			pack_int32((uint8_t *)buf + buf_len + 12, len2);
			if (n_read == m_read) {
				m_read <<= 1;
				rbc = realloc(rbc, m_read * sizeof(struct readbc_t));
			}
			rbc[n_read].barcode = barcode;
			rbc[n_read].offset = buf_len;
			rbc[n_read].len1 = len1;
			rbc[n_read].len2 = len2;
			buf_len += len1 + len2 + 16;
			++n_read;
			if (rc1 == READ_END)
				break;
		}
	}
	if (n_read) {
		sprintf(path, "%s.tmp.%d", prefix, n_file);
		write_buffer_to_file(rbc, n_read, buf, buf_len, fbuf, fbuf_size, path);
		++n_file;
		n_read = 0;
		buf_len = 0;
	}
	free(buf);
	free(fbuf);
	free(path);
	free(rbc);
	merge_sorted_small(prefix, sm, n_file);
	return NULL;
}

void merge_sorted_large(const char *prefix, int64_t sm, int n_file)
{
	int64_t sm_per_file;
	char **fbuf, *tmp_buf, *obuf1, *obuf2, *buf, *path;
	FILE **fp, *fo1, *fo2;
	struct readbc_t *reads;
	int i, m_buf;
	sm_per_file = sm / n_file;
	__round_up_64(sm_per_file);
	if (sm_per_file > sm / n_file)
		sm_per_file >>= 1;
	fp = malloc(n_file * sizeof(FILE *));
	fbuf = malloc(n_file * sizeof(char *));
	reads = malloc(n_file * sizeof(struct readbc_t));
	tmp_buf = malloc(16);
	path = malloc(MAX_PATH);
	for (i = 0; i < n_file; ++i) {
		sprintf(path, "%s/thread_%d.tmp", prefix, i);
		fp[i] = xfopen(path, "rb");
		fbuf[i] = malloc(sm_per_file);
		setbuffer(fp[i], fbuf[i], sm_per_file);
		extract_read_barcode(fp[i], reads + i, tmp_buf);
	}
	sprintf(path, "%s/R1.sorted.fq", prefix);
	fo1 = xfopen(path, "wb");
	sprintf(path, "%s/R2.sorted.fq", prefix);
	fo2 = xfopen(path, "wb");
	obuf1 = malloc(SIZE_2MB);
	obuf2 = malloc(SIZE_2MB);
	setbuffer(fo1, obuf1, SIZE_2MB);
	setbuffer(fo2, obuf2, SIZE_2MB);
	m_buf = 0x100;
	buf = malloc(m_buf);
	while (1) {
		uint64_t cur_bc = (uint64_t)-1;
		int idx = -1;
		for (i = 0; i < n_file; ++i) {
			if (reads[i].barcode != (uint64_t)-1 &&
				(cur_bc == (uint64_t)-1 || reads[i].barcode < cur_bc)) {
				cur_bc = reads[i].barcode;
				idx = i;
			}
		}
		if (cur_bc == (uint64_t)-1)
			break;
		if (reads[idx].len1 + reads[idx].len2 > m_buf) {
			m_buf = reads[idx].len1 + reads[idx].len2;
			buf = realloc(buf, m_buf);
		}
		// pack_int64((uint8_t *)tmp_buf, reads[idx].barcode);
		// pack_int32((uint8_t *)tmp_buf + 8, reads[idx].len1);
		// pack_int32((uint8_t *)tmp_buf + 12, reads[idx].len2);
		// fwrite(tmp_buf, 16, 1, fo);
		size_t byte_read = fread(buf, 1, reads[idx].len1 + reads[idx].len2, fp[idx]);
		if ((int64_t)byte_read < reads[idx].len1 + reads[idx].len2)
			__ERROR("Corrupted temporary files");
		fwrite(buf, reads[idx].len1, 1, fo1);
		fwrite(buf + reads[idx].len1, reads[idx].len2, 1, fo2);
		// fwrite(buf, reads[idx].len1 + reads[idx].len2, 1, fo);
		extract_read_barcode(fp[idx], reads + idx, tmp_buf);
	}
	for (i = 0; i < n_file; ++i) {
		sprintf(path, "%s/thread_%d.tmp", prefix, i);
		fclose(fp[i]);
		free(fbuf[i]);
		remove(path);
	}
	free(buf);
	free(path);
	free(tmp_buf);
	free(reads);
	free(fbuf);
	free(fp);
	fclose(fo1);
	fclose(fo2);
	free(obuf1);
	free(obuf2);
}

void sort_read(struct opt_proc_t *opt)
{
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	void *(*buffer_iterator)(void *) = biot_buffer_iterator;
	struct producer_bundle_t *producer_bundles = NULL;
	if (opt->lib_type == LIB_TYPE_BIOT || opt->lib_type == LIB_TYPE_10X) {
		producer_bundles = init_fastq_pair(opt->n_threads, opt->n_files,
						opt->files_1, opt->files_2);
		if (opt->lib_type == LIB_TYPE_BIOT)
			buffer_iterator = biot_buffer_iterator;
		else
			buffer_iterator = x10_buffer_iterator;
	} else if (opt->lib_type == LIB_TYPE_UST) {
		producer_bundles = init_fastq_triple(opt->n_threads, opt->n_files,
				opt->files_1, opt->files_2, opt->files_I);
		buffer_iterator = ust_buffer_iterator;
	}
	struct readsort_bundle_t *worker_bundles;
	worker_bundles = malloc(opt->n_threads * sizeof(struct readsort_bundle_t));
	pthread_mutex_t lock;
	pthread_mutex_init(&lock, NULL);
	int i;
	int64_t sm_in_byte, sm_per_thread;
	sm_in_byte = (int64_t)opt->mmem * 1024 * 1024 * 1024;
	sm_per_thread = sm_in_byte / opt->n_threads;
	__round_up_64(sm_per_thread);
	if (sm_per_thread > sm_in_byte / opt->n_threads)
		sm_per_thread >>= 1;
	char *read_dir = malloc(MAX_PATH);
	sprintf(read_dir, "%s/reads", opt->out_dir);
	mkdir(read_dir, 0755);
	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].sm = sm_per_thread;
		sprintf(worker_bundles[i].prefix, "%s/thread_%d", read_dir, i);
	}

	pthread_t *producer_threads, *worker_threads;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	for (i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_producer,
				producer_bundles + i);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, buffer_iterator,
				worker_bundles + i);

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	if (opt->lib_type == LIB_TYPE_BIOT || opt->lib_type == LIB_TYPE_10X)
		free_fastq_pair(producer_bundles, opt->n_files);
	else if (opt->lib_type == LIB_TYPE_UST)
		free_fastq_triple(producer_bundles, opt->n_files);
	free(worker_bundles);

	merge_sorted_large(read_dir, sm_in_byte, opt->n_threads);

	free(producer_threads);
	free(worker_threads);
}

