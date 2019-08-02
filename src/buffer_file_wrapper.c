#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "buffer_file_wrapper.h"
#include "io_utils.h"
#include "utils.h"

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

int64_t bf_read(struct buffered_file_t *bfp, void *ptr, int64_t sz)
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

int64_t bf_write(struct buffered_file_t *bfp, const void *ptr, int64_t sz)
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

