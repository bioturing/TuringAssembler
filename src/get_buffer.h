#ifndef _GET_BUFFER_H_
#define _GET_BUFFER_H_

#include <zlib.h>

#include "attribute.h"

#define TYPE_FASTQ		0
#define TYPE_FASTA		1
#define TYPE_

// #define BUF_OK			0
// #define BUF_FAIL		1

#define BUF_SIZE		SIZE_1MB

struct gb_file_inf {
	gzFile fp;
	int buf_size;
	int is_eof;
	char *buf;
	const char *name;
	int64_t processed;
	int64_t compressed_size;
};

/* pair */

struct gb_pair_data {
	struct gb_file_inf R1;
	struct gb_file_inf R2;
	int64_t compressed_size;
	int finish_flag;
	int warning_flag;
	int type;
	int64_t processed;
};

void gb_pair_init(struct gb_pair_data *data, const char *R1_path, const char *R2_path);
void gb_pair_destroy(struct gb_pair_data *data);
int gb_get_pair(struct gb_pair_data *data, char **R1_buf, char **R2_buf);

/* triplet-barcoding */

struct gb_trip_data {
	struct gb_file_inf R1;
	struct gb_file_inf R2;
	struct gb_file_inf I;
	int64_t compressed_size;
	int finish_flag;
	int warning_flag;
	int type;
	int64_t processed;
};

void gb_trip_init(struct gb_trip_data *data, const char *R1_path,
				const char *R2_path, const char *I_path);
void gb_trip_destroy(struct gb_trip_data *data);
int gb_get_trip(struct gb_trip_data *data, char **R1_buf, char **R2_buf, char **I_buf);

/* single */

struct gb_single_data {
	struct gb_file_inf R;
	int64_t compressed_size;
	int finish_flag;
	int type;
	int64_t processed;
};

void gb_single_init(struct gb_single_data *data, char *R_path);
void gb_single_destroy(struct gb_single_data *data);
int gb_get_single(struct gb_single_data *data, char **buf);

/* get read */

#define READ_SUCCESS		0
#define READ_END		1
#define	READ_FAIL		2

void read_destroy(struct read_t *read, int is_buf);
int get_read_from_fq(struct read_t *read, char *buf, int *pos);
int get_read_from_fa(struct read_t *read, char *buf, int *pos);

#endif /* _GET_BUFFER_H_ */
