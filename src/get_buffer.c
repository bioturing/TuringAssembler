#include <string.h>
#include "get_buffer.h"
#include "io_utils.h"
#include "utils.h"
#include "verbose.h"

static int get_format(const char *file_path)
{
	char buf[1024];
	gzFile file;
	int n, ret;

	file = gzopen(file_path, "r");
	if (!file)
		log_error("Could not open file: %s", file_path);

	n = gzread(file, buf, 1);
	if (!n)
		log_error("Could not read file: %s", file_path);

	if (buf[0] == '@')
		ret = TYPE_FASTQ;
	else if (buf[0] == '>')
		ret = TYPE_FASTA;
	else
		log_error("Unsupport format of file: %s", file_path);

	gzclose(file);
	return ret;
}

static int64_t gb_file_get_size(const char *path)
{
	FILE *fp = xfopen(path, "rb");
	fseek(fp, 0L, SEEK_END);
	int64_t ret = ftell(fp);
	fclose(fp);
	return ret;
}

static void gb_file_init(struct gb_file_inf *f, const char *path)
{
	f->name = path;
	f->fp = gzopen(path, "r");
	f->buf = malloc(BUF_SIZE + 1);
	f->buf_size = 0;
	f->is_eof = 0;
	f->processed = 0;
	f->compressed_size = gb_file_get_size(path);
}

void gb_trip_init(struct gb_trip_data *data, const char *R1_path,
				const char *R2_path, const char *I_path)
{
	if (!strcmp(R1_path, R2_path) || !strcmp(R1_path, I_path) || !strcmp(R2_path, I_path))
		log_error("Two identical read files");

	data->type = get_format(R1_path);
	if (get_format(R2_path) != data->type || get_format(I_path) != data->type)
		log_error("Error in three read files are not equal");
	gb_file_init(&(data->R1), R1_path);
	gb_file_init(&(data->R2), R2_path);
	gb_file_init(&(data->I), I_path);
	data->finish_flag = 0;
	data->warning_flag = 0;
	data->compressed_size = data->R1.compressed_size + data->R2.compressed_size + data->I.compressed_size;
}

void gb_trip_destroy(struct gb_trip_data *data)
{
	free(data->R1.buf);
	free(data->R2.buf);
	free(data->I.buf);
	gzclose(data->R1.fp);
	gzclose(data->R2.fp);
	gzclose(data->I.fp);
}

void gb_pair_init(struct gb_pair_data *data, const char *R1_path, const char *R2_path)
{
	if (!strcmp(R1_path, R2_path))
		log_error("Two identical read files");

	data->type = get_format(R1_path);
	if (get_format(R2_path) != data->type)
		log_error("Format in two read files are not equal");

	gb_file_init(&(data->R1), R1_path);
	gb_file_init(&(data->R2), R2_path);
	data->finish_flag = 0;
	data->warning_flag = 0;
	data->compressed_size = data->R1.compressed_size + data->R2.compressed_size;
}

void gb_pair_destroy(struct gb_pair_data *data)
{
	free(data->R1.buf);
	free(data->R2.buf);
	gzclose(data->R1.fp);
	gzclose(data->R2.fp);
}

void gb_single_init(struct gb_single_data *data, char *file_path)
{
	data->type = get_format(file_path);
	gb_file_init(&(data->R), file_path);
	data->finish_flag = 0;
	data->compressed_size = data->R.compressed_size;
}

void gb_single_destroy(struct gb_single_data *data)
{
	free(data->R.buf);
	gzclose(data->R.fp);
}

/*
 * return 0 if read is not found
 * otherwise return the next read's position
 */
static inline int get_next_pos(struct gb_file_inf *data, int prev, int type)
{
	int i = prev;
	char *buf = data->buf;
	int size = data->buf_size;
	int id = 0;
	int mod = (type == TYPE_FASTQ ? 3 : 1);

	while (1) {
		for (; buf[i] != '\n' && i < size; ++i);
		if (i == size) {
			if (!data->is_eof)
				break;
		} else {
			++i;
		}
		id = (id + 1) & mod;
		if (id == 0)
			return i;
		if (i == size)
			break;
	}

	if (size - prev >= BUF_SIZE) {
		__VERBOSE("\n");
		log_error("Read is too long from file: %s", data->name);
	}

	return 0;
}

static void load_buffer(struct gb_file_inf *data)
{
	int size = BUF_SIZE - data->buf_size;
	int n_byte = gzread(data->fp, data->buf + data->buf_size, size);
	data->buf_size += n_byte;
	data->processed = gzoffset(data->fp);
	if (n_byte < size)
		data->is_eof = 1;
}

static void split_buffer(struct gb_file_inf *data, char *old_buf, int prev)
{
	int padding = data->buf_size - prev;
	/* BUF_SIZE + 1 for case end of file is not '/n' */
	memcpy(data->buf, old_buf + prev, padding);
	data->buf_size = padding;
	data->buf[padding] = '\0';
	old_buf[prev] = '\0';
}

int64_t gb_get_trip(void *vdata, void *vp)
{
	struct gb_trip_data *data = (struct gb_trip_data *)vdata;
	struct trip_buffer_t *p = (struct trip_buffer_t *)vp;
	if (data->finish_flag)
		return -1;

	int prev1, prev2, prevI, new_pos1, new_pos2, new_posI;

	load_buffer(&data->R1);
	load_buffer(&data->R2);
	load_buffer(&data->I);
	data->processed = data->R1.processed + data->R2.processed + data->I.processed;
	prev1 = prev2 = prevI = 0;

	while (1) {
		new_pos1 = get_next_pos(&data->R1, prev1, data->type);
		new_pos2 = get_next_pos(&data->R2, prev2, data->type);
		new_posI = get_next_pos(&data->I, prevI, data->type);
		if (!new_pos1 || !new_pos2 || !new_posI)
			break;
		prev1 = new_pos1;
		prev2 = new_pos2;
		prevI = new_posI;
	}

	/* no more read to get */
	if (!prev1 && !new_pos1 && !new_pos2 && !new_posI) {
		data->finish_flag = 1;
		return -1;
	}

	/* read of two files are not equal */
	if (!prev1 || !prev2 || !prevI) {
		data->warning_flag = 1;
		data->finish_flag = 1;
		return -2;
	}

	/* split complete buffer */
	__SWAP(p->R1_buf, data->R1.buf);
	__SWAP(p->R2_buf, data->R2.buf);
	__SWAP(p->I_buf, data->I.buf);
	split_buffer(&data->R1, p->R1_buf, prev1);
	split_buffer(&data->R2, p->R2_buf, prev2);
	split_buffer(&data->I, p->I_buf, prevI);
	p->input_format = data->type;
	return data->processed;
}

int64_t gb_get_pair(void *vdata, void *vp)
{
	struct gb_pair_data *data = (struct gb_pair_data *)vdata;
	struct pair_buffer_t *p = (struct pair_buffer_t *)vp;
	if (data->finish_flag)
		return -1;

	int prev1, prev2, new_pos1, new_pos2;

	load_buffer(&data->R1);
	load_buffer(&data->R2);
	data->processed = data->R1.processed + data->R2.processed;
	prev1 = prev2 = 0;

	while (1) {
		new_pos1 = get_next_pos(&data->R1, prev1, data->type);
		new_pos2 = get_next_pos(&data->R2, prev2, data->type);
		if (!new_pos1 || !new_pos2)
			break;
		prev1 = new_pos1;
		prev2 = new_pos2;
	}

	/* no more read to get */
	if (!prev1 && !prev2 && !new_pos1 && !new_pos2) {
		data->finish_flag = 1;
		return -1;
	}

	/* read of two files are not equal */
	if (!prev1 || !prev2) {
		data->warning_flag = 1;
		data->finish_flag = 1;
		return -2;
	}

	/* split complete buffer */
	__SWAP(p->R1_buf, data->R1.buf);
	__SWAP(p->R2_buf, data->R2.buf);
	split_buffer(&data->R1, p->R1_buf, prev1);
	split_buffer(&data->R2, p->R2_buf, prev2);
	p->input_format = data->type;
	return data->processed;
}

int64_t gb_get_single(void *vdata, void *vp)
{
	struct gb_single_data *data = (struct gb_single_data *)vdata;
	struct single_buffer_t *p = (struct single_buffer_t *)vp;
	if (data->finish_flag)
		return -1;

	int prev, new_pos;

	load_buffer(&data->R);
	data->processed = data->R.processed;
	prev = 0;

	while (1) {
		new_pos = get_next_pos(&data->R, prev, data->type);
		if (!new_pos) 
			break;
		prev = new_pos;
	}

	/* no more read to get */
	if (!prev) {
		data->finish_flag = 1;
		return -1;
	}

	/* split complete buffer */
	__SWAP(p->R_buf, data->R.buf);
	split_buffer(&data->R, p->R_buf, prev);
	p->input_format = data->type;
	return data->processed;
}

void read_destroy(struct read_t *read, int is_allocated)
{
	if (!is_allocated) {
		free(read->seq);
		free(read->qual);
		free(read->name);
		free(read->note);
	}
}

int get_read_from_fq(struct read_t *read, char *buf, int *pos)
{
	int i = *pos, prev, k = -1;

	/* name part */
	prev = i;
	if (buf[i] != '@')
		return READ_FAIL;
	read->info = NULL;
	read->name = buf + prev + 1; /* skip @ character */
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i) {
		if (__is_sep(buf[i]) && read->info == NULL) {
			buf[i] = '\0';
			read->info = buf + i + 1;
			k = i - 2;
		}
	}
	if (buf[i] == '\0')
		return READ_FAIL;
	if (read->info == NULL)
		k = i - 2;

	if (k > 0 &&
	    (strncmp(buf + k, "/1", 2) == 0 || strncmp(buf + k, "/2", 2) == 0)) {
		buf[k] = '\0';
	}

	buf[i++] = '\0';

	/* seq part */
	prev = i;
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i);
	if (buf[i] == '\0')
		return READ_FAIL;
	read->seq = buf + prev;
	read->len = i - prev;
	buf[i++] = '\0';
	if (read->len == 0)
		return READ_FAIL;

	/* optionally part */
	prev = i;
	if (buf[i] != '+')
		return READ_FAIL;
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i);
	if (buf[i] == '\0')
		return READ_FAIL;
	read->note = buf + prev;
	buf[i++] = '\0';

	/* quality part */
	prev = i;
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i);
	if (i - prev != read->len)
		return READ_FAIL;
	read->qual = buf + prev;
	if (buf[i] == '\0')
		return READ_END;
	buf[i++] = '\0';
	if (buf[i] == '\0')
		return READ_END;

	*pos = i;
	return READ_SUCCESS;
}

int get_read_from_fa(struct read_t *read, char *buf, int *pos)
{
	int i = *pos, prev, k = -1;
	read->qual = read->note = NULL;

	/* name part */
	prev = i;
	if (buf[i] != '>')
		return READ_FAIL;
	read->info = NULL;
	read->name = buf + prev + 1; /* skip > character */
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i) {
		if (__is_sep(buf[i]) && read->info == NULL) {
			buf[i] = '\0';
			read->info = buf + i + 1;
			k = i - 2;
		}
	}
	if (buf[i] == '\0')
		return READ_FAIL;
	if (read->info == NULL)
		k = i - 2;

	if (k > 0 &&
	    (strncmp(buf + k, "/1", 2) == 0 || strncmp(buf + k, "/2", 2) == 0)) {
		buf[k] = '\0';
	}

	buf[i++] = '\0';

	/* seq part */
	prev = i;
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i);
	read->seq = buf + prev;
	read->len = i - prev;
	if (read->len == 0)
		return READ_FAIL;
	if (buf[i] == '\0')
		return READ_END;
	buf[i++] = '\0';
	if (buf[i] == '\0')
		return READ_END;

	*pos = i;
	return READ_SUCCESS;
}
