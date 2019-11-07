#ifndef __BUFFER_FILE_WRAPPER_H__
#define __BUFFER_FILE_WRAPPER_H__

#define BF_READ		0
#define BF_WRITE	1

struct buffered_file_t {
	FILE *fp;
	char *buf;
	int64_t buf_len;
	int64_t buf_size;
	int64_t buf_pos;
	int mode;
};

void bf_open(struct buffered_file_t *bfp, const char *path, const char *mode, int64_t buf_size);
void bf_close(struct buffered_file_t *bfp);

int64_t bf_read(struct buffered_file_t *bfp, void *ptr, int64_t sz);
int64_t bf_write(struct buffered_file_t *bfp, const void *ptr, int64_t sz);
int64_t check_data(const void *ptr, int64_t sz, int64_t bc);

#endif /* __BUFFER_FILE_WRAPPER__ */
