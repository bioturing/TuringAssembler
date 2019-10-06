//
// Created by BioTuring on 2019-09-25.
//

#include <zlib.h>

#include "smart_load.h"
#include "sort_read.h"
#include "verbose.h"
#include "radix_sort.h"
#include "buffer_file_wrapper.h"
#include "io_utils.h"
#include "kseq.h"
#include "utils.h"

/**
 * Construct the dictionary index of
 * barcode position in the barcode sorted file
 * @param rpath pointer to the read_path_t struct, which contains path to R1 and R2
 * @param h khash struct that whose barcodes indices. Will be used to look for position of any barcode
 */
void smart_construct_read_index(struct read_path_t *rpath, khash_t(bcpos) *h)
{
	FILE *fp = xfopen(rpath->idx_path, "rb");
	size_t byte_read;
	khint_t k;
	int ret;
	uint64_t barcode;
	char *buf = alloca(40);
	while ((byte_read = fread(buf, 1, 40, fp))) {
		if (byte_read != 40)
			__ERROR("Corrupted barcode in read index file");
		barcode = unpack_int64((uint8_t *)buf);
		k = kh_put(bcpos, h, barcode, &ret);
		if (ret != 1)
			__ERROR("Insert barcode failed");
		kh_value(h, k).r1_offset = unpack_int64((uint8_t *)buf + 8);
		kh_value(h, k).r2_offset = unpack_int64((uint8_t *)buf + 16);
		kh_value(h, k).r1_len = unpack_int64((uint8_t *)buf + 24);
		kh_value(h, k).r2_len = unpack_int64((uint8_t *)buf + 32);
	}
	fclose(fp);
}

/**
 * @brief Reads of one set of barcodes into string stream
 * @param ref read_path struct that contain path to original read1 and read2 (maybe barcode included)
 * @param dict khash_t(bcpost) hash table to lookup the offset of one barcode in the sorted data
 * @param shared array of barcodes pseudo hash values
 * @param n_shared number of barcodes in shared array
 * @param buf1 buffer stream for read 1
 * @param buf2 buffer stream for read 2
 * @param size1 size of the buffer stream for read 1
 * @param size2 size of the buffer stream for read 2
 */
void stream_filter_read(struct read_path_t *ref, khash_t(bcpos) *dict,
                 uint64_t *shared, int n_shared, char **buf1, char **buf2,
                 uint64_t *size1, uint64_t *size2)
{
	struct read_index_t *pos;
	pos = malloc(n_shared * sizeof(struct read_index_t));
	khiter_t k;
	int i;
	for (i = 0; i < n_shared; ++i) {
		k = kh_get(bcpos, dict, shared[i]);
		pos[i] = kh_value(dict, k);
	}
	rs_sort(read_index, pos, pos + n_shared); //TODO: get rid of this function
	FILE *fi1 = xfopen(ref->R1_path, "rb");
	FILE *fi2 = xfopen(ref->R2_path, "rb");
	uint64_t m_buf1 = 0, m_buf2 = 0;
	uint64_t prev1 = 0 , prev2 = 0;

	for (i = 0 ; i < n_shared; ++i) {
		m_buf1 += pos[i].r1_len;
		m_buf2 += pos[i].r2_len;
	}
	char *buf1_ = malloc(m_buf1);
	char *buf2_ = malloc(m_buf2);

	for (i = 0; i < n_shared; ++i) {
		fseek(fi1, pos[i].r1_offset, SEEK_SET);
		xfread(buf1_ + prev1, 1, pos[i].r1_len, fi1);
		fseek(fi2, pos[i].r2_offset, SEEK_SET);
		xfread(buf2_ + prev2, 1, pos[i].r2_len, fi2);
		prev1 += pos[i].r1_len ;
		prev2 += pos[i].r2_len ;
	}
	fclose(fi1);
	fclose(fi2);
	free(pos);
	*buf1 = buf1_;
	*buf2 = buf2_;
	*size1 = m_buf1;
	*size2 = m_buf2;
}
/**
 * @brief Seeking and loading the barcode sequences into the buffer stream
 * @param opt
 */
void smart_load_barcode(struct opt_proc_t *opt)
{
	struct read_path_t read_sorted_path;
	if (opt->lib_type == LIB_TYPE_SORTED) {
		read_sorted_path.R1_path = opt->files_1[0];
		read_sorted_path.R2_path = opt->files_2[0];
		read_sorted_path.idx_path = opt->files_I[0];
	} else {
		__VERBOSE("Reads are not sorted. Sort reads by barcode sequence...");
		sort_read(opt, &read_sorted_path);
	}
	uint64_t bx[2] = {42983706558, 106615444678}; //43 15 mock barcode pseudo hash id here
	khash_t(bcpos) *bx_pos_dict = kh_init(bcpos);
	smart_construct_read_index(&read_sorted_path, bx_pos_dict); //load the barcode indices

	char *buf1, *buf2;
	uint64_t m_buf1, m_buf2;
	stream_filter_read(&read_sorted_path, bx_pos_dict, bx, 2, &buf1, &buf2, &m_buf1, &m_buf2);

	// Test string stream
	FILE * buf1_stream = fmemopen((void *)buf1, m_buf1, "r");
	FILE * buf2_stream = fmemopen((void *)buf2, m_buf2, "r");

	char *line;
	size_t len, read;
	while ((read = getline(&line, &len, buf1_stream)) != -1) {
		printf("Retrieved line of length %zu :\n", read);
		printf("%s", line);
	}
}
