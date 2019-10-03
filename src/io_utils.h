#ifndef _IO_UTILS_H_
#define _IO_UTILS_H_

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

/* Reconstruct goal:
 * io_utils.h:	shared_fstream_t
 *		align_file_t
 *		
 */
// #if defined()
// #define RED			""
// #define GREEN			""
// #define YELLOW			""
// #define BLUE			""
// #define MAGENTA			""
// #define CYAN			""
// #define WHITE			""
// #define RESET			""
// #else
// #define RED			"\x1B[31m"
// #define GREEN			"\x1B[32m"
// #define YELLOW			"\x1B[33m"
// #define BLUE			"\x1B[34m"
// #define MAGENTA			"\x1B[35m"
// #define CYAN			"\x1B[36m"
// #define WHITE			"\x1B[37m"
// #define RESET			"\x1B[0m"
// #endif

/* shared stream struct
 * Usage: multiple threads write data to single file, provide each thread a
 * buffer
 */
#define SFS_BUF_SZ		2097152

struct shared_fstream_t {
	FILE *fp;
	pthread_mutex_t *lock;
	char *buf;
	int buf_len;
};


/* support openning file with UNICODE path */
FILE *xfopen(const char *file_path, const char *mode);

/* fflush before close file */
void xwfclose(FILE *f);

/* check fread function read enough nmemb */
size_t xfread(void *ptr, size_t size, size_t nmemb, FILE *stream);

/* check fwrite function write enough nmemb */
size_t xfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);

/* auto remove /n character if found */
ssize_t xgetline(char **str, size_t *size, FILE *stream);

/* remove redundant / character */
void normalize_dir(char *path);

/* make directory if is not exist */
void make_dir(const char *path);

/* get size of all file from file_path */
size_t fetch_size(char **file_path, int n_file);

/* ------------ shared_stream_t utils ------------ */

/* Init shared stream on n threads */
struct shared_fstream_t *init_shared_stream(const char *path, int n);

/* Destroy shared stream on n threads */
void destroy_shared_stream(struct shared_fstream_t *p, int n);

/* flush buffer to output */
void sfs_flush(struct shared_fstream_t *p);

int check_file_empty(char *file_path);
#endif
