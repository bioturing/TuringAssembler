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

struct buffered_file_t {
	FILE *fp;
	char *buf;
	int64_t buf_len;
	int64_t buf_size;
	int64_t buf_pos;
	int mode;
};

#define BF_READ		0
#define BF_WRITE	1

void bf_open(struct buffered_file_t *bfp, const char *path, const char *mode, int64_t buf_size)
{
	bfp->fp = xfopen(path, mode);
	bfp->buf_size = buf_size;
	bfp->buf = malloc(buf_size);
	if (!strcmp(mode, "rb") || !strcmp(mode, "r")) {
		bfp->buf_len = fread(bfp->buf, 1, bfp->buf_size, bfp->fp);
		bfp->mode = BF_READ;
	} else {
		bfp->buf_len = 0;
		bfp->mode = BF_WRITE;
	}
	bfp->buf_pos = 0;
}

static inline int64_t bf_read(struct buffered_file_t *bfp, void *ptr, int64_t sz)
{
	if (bfp->buf_len - bfp->buf_pos >= sz) {
		memcpy(ptr, bfp->buf + bfp->buf_pos, sz);
		bfp->buf_pos += sz;
		return sz;
	}
	int64_t rem_fill, rem_buf, to_fill;
	rem_fill = sz;
	do {
		rem_buf = bfp->buf_len - bfp->buf_pos;
		if (rem_buf == 0 && bfp->buf_len < bfp->buf_size)
			return sz - rem_fill; /* EOF */
		to_fill = __min(rem_buf, rem_fill);
		memcpy(ptr + sz - rem_fill, bfp->buf + bfp->buf_pos, to_fill);
		rem_fill -= to_fill;
		rem_buf -= to_fill;
		bfp->buf_pos += to_fill;
		if (rem_buf == 0) {
			bfp->buf_len = fread(bfp->buf, 1, bfp->buf_size, bfp->fp);
			bfp->buf_pos = 0;
		}
	} while (rem_fill);
	return sz;
}

static inline int64_t bf_write(struct buffered_file_t *bfp, const void *ptr, int64_t sz)
{
	if (bfp->buf_size - bfp->buf_len >= sz) {
		memcpy(bfp->buf + bfp->buf_len, ptr, sz);
		bfp->buf_len += sz;
		return sz;
	}
	int64_t rem_buf, rem_fill, to_fill;
	rem_fill = sz;
	do {
		rem_buf = bfp->buf_size - bfp->buf_len;
		to_fill = __min(rem_buf, rem_fill);
		memcpy(bfp->buf + bfp->buf_len, ptr + sz - rem_fill, to_fill);
		rem_fill -= to_fill;
		rem_buf -= to_fill;
		bfp->buf_len += to_fill;
		if (rem_buf == 0) {
			xfwrite(bfp->buf, 1, bfp->buf_len, bfp->fp);
			bfp->buf_len = 0;
		}
	} while (rem_fill);
	return sz;
}

void bf_close(struct buffered_file_t *bfp)
{
	if (bfp->mode == BF_WRITE)
		xfwrite(bfp->buf, 1, bfp->buf_len, bfp->fp);
	free(bfp->buf);
	fclose(bfp->fp);
}

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

static inline uint64_t get_barcode_biot(char *s, struct read_t *r)
{
	if (s == NULL) {
		r->seq = r->qual = NULL;
		return (uint64_t)-1;
	}
	int i, len = 0;
	uint64_t ret = 0;
	char *p;
	p = strstr(s, "BX:Z:");
	if (p == NULL) {
		r->seq = r->qual = NULL;
		return (uint64_t)-1;
	}
	r->seq = p += 5;
	for (i = 0; p[i] && !__is_sep(p[i]); ++i) {
		ret = ret * 5 + nt4_table[(int)p[i]];
		++len;
	}
	p = strstr(s, "QB:Z:");
	if (p == NULL)
		r->qual = NULL;
	else
		r->qual = p + 5;
	r->len = len;
	return ret;
}

static inline uint64_t get_barcode_10x(struct read_t *r1, struct read_t *rbc)
{
	if (r1->len < 23) {
		rbc->seq = rbc->qual = NULL;
		return (uint64_t)-1;
	}
	uint64_t ret = 0;
	int i;
	char *s = r1->seq;
	for (i = 0; i < 16; ++i)
		ret = ret * 5 + nt4_table[(int)s[i]];
	rbc->seq = r1->seq;
	rbc->qual = r1->qual;
	rbc->len = 16;
	r1->seq += 23;
	r1->qual += 23;
	return ret;
}

static inline void write_buffer_to_file(struct readbc_t *p, int n, char *buf,
					int64_t file_buf_sz, const char *path)
{
	/* sort p */
	rs_sort(read_sort, p, p + n);
	/*******/
	struct buffered_file_t fp;
	bf_open(&fp, path, "wb", file_buf_sz);
	int i;
	for (i = 0; i < n; ++i) {
		bf_write(&fp, buf + p[i].offset, p[i].len1 + p[i].len2 + 16);
	}
	bf_close(&fp);
}

static inline int ust_add_record(struct read_t *r, struct read_t *rI,
						char *buf, int input_format)
{
	int len = 0, tlen, i;
	tlen = strlen(r->name);
	memcpy(buf, r->name, tlen);
	len = tlen;
	if (rI->seq != NULL) {
		buf[len++] = '\t';
		memcpy(buf + len, "BX:Z:", 5);
		len += 5;
		memcpy(buf + len, rI->seq, rI->len);
		len += rI->len;
		if (input_format == TYPE_FASTQ && rI->qual != NULL) {
			buf[len++] = '\t';
			memcpy(buf + len, "QB:Z:", 5);
			len += 5;
			memcpy(buf + len, rI->qual, rI->len);
			len += rI->len;
		} else {
			for (i = 0; i < rI->len; ++i) {
				buf[len++] = nt4_table[(int)rI->seq[i]] < 4 ? 'F' : '#';
			}
		}
	}
	buf[len++] = '\n';
	memcpy(buf + len, r->seq, r->len);
	len += r->len;
	buf[len++] = '\n';
	memcpy(buf + len, "+\n", 2);
	len += 2;
	if (input_format == TYPE_FASTQ && r->qual != NULL) {
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

static inline void extract_read_barcode(struct buffered_file_t *fp, struct readbc_t *r, char *buf)
{
	int64_t ret = bf_read(fp, buf, 16);
	if (ret == 0) {
		r->barcode = (uint64_t)-1;
	} else if (ret < 16) {
		__ERROR("Corrupted temporary files %lu", ret);
	} else {
		r->barcode = unpack_int64((uint8_t *)buf);
		r->len1 = unpack_int32((uint8_t *)buf + 8);
		r->len2 = unpack_int32((uint8_t *)buf + 12);
	}
}

void merge_sorted_small(const char *prefix, int64_t sm, int n_file)
{
	int64_t sm_per_file;
	char *tmp_buf, *buf, *path;
	struct buffered_file_t *fp, fo;
	fp = malloc(n_file * sizeof(struct buffered_file_t));
	struct readbc_t *reads;
	int i, m_buf;
	sm_per_file = sm / (n_file + 1);
	__round_up_64(sm_per_file);
	if (sm_per_file * (n_file + 1) > sm + SIZE_16MB)
		sm_per_file >>= 1;
	reads = malloc(n_file * sizeof(struct readbc_t));
	tmp_buf = malloc(16);
	path = malloc(MAX_PATH);
	for (i = 0; i < n_file; ++i) {
		sprintf(path, "%s.tmp.%d", prefix, i);
		bf_open(fp + i, path, "rb", sm_per_file);
		extract_read_barcode(fp + i, reads + i, tmp_buf);
	}
	sprintf(path, "%s.tmp", prefix);
	bf_open(&fo, path, "wb", sm_per_file);
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
		bf_write(&fo, tmp_buf, 16);
		int64_t byte_read = bf_read(fp + idx, buf, reads[idx].len1 + reads[idx].len2);
		if (byte_read < reads[idx].len1 + reads[idx].len2)
			__ERROR("small merge Corrupted temporary files");
		bf_write(&fo, buf, reads[idx].len1 + reads[idx].len2);
		extract_read_barcode(fp + idx, reads + idx, tmp_buf);
	}
	for (i = 0; i < n_file; ++i) {
		sprintf(path, "%s.tmp.%d", prefix, i);
		bf_close(fp + i);
		remove(path);
	}
	free(buf);
	free(path);
	free(tmp_buf);
	free(reads);
	free(fp);
	bf_close(&fo);
}

void *biot_buffer_iterator(void *data)
{
	struct readsort_bundle_t *bundle = (struct readsort_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2, readbc;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_trip_buffer();
	char *prefix = bundle->prefix;
	int64_t sm = bundle->sm;
	char *path = malloc(MAX_PATH);

	char *buf, *fbuf;
	int64_t buf_len, buf_size;
	buf = malloc(sm);
	buf_len = 0;
	buf_size = sm;

	int m_read, n_read, n_file;
	struct readbc_t *rbc;
	m_read = sm / 600;
	__round_up_32(m_read);
	rbc = malloc(m_read * sizeof(struct readbc_t));
	n_read = 0;
	n_file = 0;

	char *R1_buf, *R2_buf;
	int pos1, pos2, rc1, rc2, input_format;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf)
			break;
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		R1_buf = ext_buf->R1_buf;
		R2_buf = ext_buf->R2_buf;
		input_format = ext_buf->input_format;

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, R1_buf, &pos1) :
				get_read_from_fa(&read1, R1_buf, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, R2_buf, &pos2) :
				get_read_from_fa(&read2, R2_buf, &pos2);

			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file\n");

			/* read_name + \t + BX:Z: + barcode + \t + QB:Z: + barcode_quality + \n */
			uint64_t barcode = get_barcode_biot(read1.info, &readbc);
			int record_len, len1, len2;
			if (barcode != (uint64_t)-1) {
				record_len = strlen(read1.name) + strlen(read2.name) +
					2 * (1 + (readbc.len + 6) * 2) +
					(read1.len + 1) * 2 + (read2.len + 1) * 2 + 2 * 2 +
					8 + 4 * 2;
			} else {
				record_len = strlen(read1.name) + strlen(read2.name) +
					(read1.len + 1) * 2 + (read2.len + 1) * 2 + 2 * 2 +
					8 + 4 * 2;
			}
			if (record_len > buf_size) {
				/* extend buf size */
				buf_size = record_len;
				buf = realloc(buf, buf_size);
			}
			if (record_len > buf_size - buf_len) {
				/* sort and flush buffer */;
				sprintf(path, "%s.tmp.%d", prefix, n_file);
				write_buffer_to_file(rbc, n_read, buf, SIZE_16MB, path);
				++n_file;
				n_read = 0;
				buf_len = 0;
			}
			/* append record */
			pack_int64((uint8_t *)buf + buf_len, barcode);
			len1 = ust_add_record(&read1, &readbc,
					buf + buf_len + 16, input_format);
			len2 = ust_add_record(&read2, &readbc,
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
		write_buffer_to_file(rbc, n_read, buf, SIZE_16MB, path);
		++n_file;
		n_read = 0;
		buf_len = 0;
	}
	free(buf);
	free(path);
	free(rbc);
	merge_sorted_small(prefix, sm, n_file);
	return NULL;
}

void *x10_buffer_iterator(void *data)
{
	struct readsort_bundle_t *bundle = (struct readsort_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2, readbc;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_trip_buffer();
	char *prefix = bundle->prefix;
	int64_t sm = bundle->sm;
	char *path = malloc(MAX_PATH);

	char *buf, *fbuf;
	int64_t buf_len, buf_size;
	buf = malloc(sm);
	buf_len = 0;
	buf_size = sm;

	int m_read, n_read, n_file;
	struct readbc_t *rbc;
	m_read = sm / 600;
	__round_up_32(m_read);
	rbc = malloc(m_read * sizeof(struct readbc_t));
	n_read = 0;
	n_file = 0;

	char *R1_buf, *R2_buf;
	int pos1, pos2, rc1, rc2, input_format;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf)
			break;
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		R1_buf = ext_buf->R1_buf;
		R2_buf = ext_buf->R2_buf;
		input_format = ext_buf->input_format;

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, R1_buf, &pos1) :
				get_read_from_fa(&read1, R1_buf, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, R2_buf, &pos2) :
				get_read_from_fa(&read2, R2_buf, &pos2);

			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file\n");

			/* read_name + \t + BX:Z: + barcode + \t + QB:Z: + barcode_quality + \n */
			uint64_t barcode = get_barcode_10x(&read1, &readbc);
			int record_len, len1, len2;
			if (barcode != (uint64_t)-1) {
				record_len = strlen(read1.name) + strlen(read2.name) +
					2 * (1 + (readbc.len + 6) * 2) +
					(read1.len + 1) * 2 + (read2.len + 1) * 2 + 2 * 2 +
					8 + 4 * 2;
			} else {
				record_len = strlen(read1.name) + strlen(read2.name) +
					(read1.len + 1) * 2 + (read2.len + 1) * 2 + 2 * 2 +
					8 + 4 * 2;
			}
			if (record_len > buf_size) {
				/* extend buf size */
				buf_size = record_len;
				buf = realloc(buf, buf_size);
			}
			if (record_len > buf_size - buf_len) {
				/* sort and flush buffer */;
				sprintf(path, "%s.tmp.%d", prefix, n_file);
				write_buffer_to_file(rbc, n_read, buf, SIZE_16MB, path);
				++n_file;
				n_read = 0;
				buf_len = 0;
			}
			/* append record */
			pack_int64((uint8_t *)buf + buf_len, barcode);
			len1 = ust_add_record(&read1, &readbc,
					buf + buf_len + 16, input_format);
			len2 = ust_add_record(&read2, &readbc,
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
		write_buffer_to_file(rbc, n_read, buf, SIZE_16MB, path);
		++n_file;
		n_read = 0;
		buf_len = 0;
	}
	free(buf);
	free(path);
	free(rbc);
	merge_sorted_small(prefix, sm, n_file);
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

	char *buf, *fbuf;
	int64_t buf_len, buf_size;
	buf = malloc(sm);
	buf_len = 0;
	buf_size = sm;

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
				write_buffer_to_file(rbc, n_read, buf, SIZE_16MB, path);
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
		write_buffer_to_file(rbc, n_read, buf, SIZE_16MB, path);
		++n_file;
		n_read = 0;
		buf_len = 0;
	}
	free(buf);
	free(path);
	free(rbc);
	merge_sorted_small(prefix, sm, n_file);
	return NULL;
}

void merge_sorted_large(const char *prefix, int64_t sm, int n_file)
{
	int64_t sm_per_file;
	char *tmp_buf, *buf, *path;
	struct buffered_file_t *fp, fo1, fo2, f_idx;
	fp = malloc(n_file * sizeof(struct buffered_file_t));
	struct readbc_t *reads;
	int i, m_buf;
	sm_per_file = sm / (n_file + 2);
	__round_up_64(sm_per_file);
	if (sm_per_file * (n_file + 2) > sm + SIZE_16MB)
		sm_per_file >>= 1;
	reads = malloc(n_file * sizeof(struct readbc_t));
	tmp_buf = malloc(40);
	path = malloc(MAX_PATH);
	for (i = 0; i < n_file; ++i) {
		sprintf(path, "%s/thread_%d.tmp", prefix, i);
		bf_open(fp + i, path, "rb", sm_per_file);
		extract_read_barcode(fp + i, reads + i, tmp_buf);
	}
	sprintf(path, "%s/R1.sorted.fq", prefix);
	bf_open(&fo1, path, "wb", sm_per_file);
	sprintf(path, "%s/R2.sorted.fq", prefix);
	bf_open(&fo2, path, "wb", sm_per_file);
	sprintf(path, "%s/barcode.idx", prefix);
	bf_open(&f_idx, path, "wb", SIZE_2MB);
	m_buf = 0x100;
	buf = malloc(m_buf);
	int64_t offset_R1, offset_R2, poffset_R1, poffset_R2;
	offset_R1 = offset_R2 = poffset_R1 = poffset_R2 = 0;
	uint64_t pbarcode, cur_bc;
	pbarcode = (uint64_t)-1;
	while (1) {
		cur_bc = (uint64_t)-1;
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
		int64_t byte_read = bf_read(fp + idx, buf, reads[idx].len1 + reads[idx].len2);
		if (byte_read < reads[idx].len1 + reads[idx].len2)
			__ERROR("[large merge] Corrupted temporary files");
		bf_write(&fo1, buf, reads[idx].len1);
		bf_write(&fo2, buf + reads[idx].len1, reads[idx].len2);
		offset_R1 += reads[idx].len1;
		offset_R2 += reads[idx].len2;
		if (cur_bc != pbarcode && pbarcode != (uint64_t)-1) {
			pack_int64((uint8_t *)tmp_buf, pbarcode);
			pack_int64((uint8_t *)tmp_buf + 8, poffset_R1);
			pack_int64((uint8_t *)tmp_buf + 16, poffset_R2);
			pack_int64((uint8_t *)tmp_buf + 24, offset_R1 - poffset_R1);
			pack_int64((uint8_t *)tmp_buf + 32, offset_R2 - poffset_R2);
			bf_write(&f_idx, tmp_buf, 40);
			poffset_R1 = offset_R1;
			poffset_R2 = offset_R2;
		}
		pbarcode = cur_bc;
		extract_read_barcode(fp + idx, reads + idx, tmp_buf);
	}
	pack_int64((uint8_t *)tmp_buf, pbarcode);
	pack_int64((uint8_t *)tmp_buf + 8, poffset_R1);
	pack_int64((uint8_t *)tmp_buf + 16, poffset_R2);
	pack_int64((uint8_t *)tmp_buf + 24, offset_R1 - poffset_R1);
	pack_int64((uint8_t *)tmp_buf + 32, offset_R2 - poffset_R2);
	bf_write(&f_idx, tmp_buf, 40);
	for (i = 0; i < n_file; ++i) {
		sprintf(path, "%s/thread_%d.tmp", prefix, i);
		bf_close(fp + i);
		remove(path);
	}
	free(buf);
	free(path);
	free(tmp_buf);
	free(reads);
	free(fp);
	bf_close(&fo1);
	bf_close(&fo2);
	bf_close(&f_idx);
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

