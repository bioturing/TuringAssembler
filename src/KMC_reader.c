#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "attribute.h"
#include "assembly_graph.h"
#include "io_utils.h"
#include "KMC_reader.h"
#include "utils.h"
#include "verbose.h"

void destroy_kmc_info(struct kmc_info_t *inf)
{
	free(inf->prefix_offset);
	free(inf->sig_map);
	inf->prefix_offset = NULL;
	inf->sig_map = NULL;
}

void KMC_read_prefix(const char *path, struct kmc_info_t *data)
{
	FILE *fp = xfopen(path, "rb");
	char sig[4];
	/* read signature to check whether not truncated file */
	fseek(fp, 0, SEEK_SET);
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, KMCP, 4))
		__ERROR("Wrong KMC prefix file");
	fseek(fp, -4, SEEK_END);
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, KMCP, 4))
		__ERROR("KMC prefix file is truncated");
	fseek(fp, 0, SEEK_END);
	uint64_t file_size = ftell(fp);

	/* read header position */
	int header_position;
	fseek(fp, -8, SEEK_END);
	xfread(&header_position, 4, 1, fp);

	/* read header */
	fseek(fp, -header_position - 8, SEEK_END);
	xfread(&(data->header), sizeof(struct kmc_header_t), 1, fp);

	/* calculate some auxiliary information */
	data->sig_map_size = ((uint64_t)1 << (2 * data->header.signature_length)) + 1;
	uint64_t sig_map_size_byte = data->sig_map_size * sizeof(uint32_t);
	uint64_t prefix_offset_size_byte = file_size - sig_map_size_byte - header_position - 12;
	data->prefix_offset_size = (prefix_offset_size_byte - 8) / sizeof(uint64_t);

	/* read prefix offset */
	data->prefix_offset = malloc(prefix_offset_size_byte);
	fseek(fp, 4, SEEK_SET);
	xfread(data->prefix_offset, 1, prefix_offset_size_byte, fp);
	data->prefix_offset[data->prefix_offset_size] = data->header.total_kmers + 1;

	/* read signature map */
	data->sig_map = malloc(sig_map_size_byte);
	xfread(data->sig_map, data->sig_map_size, sizeof(uint32_t), fp);
	fclose(fp);
}

struct reader_buffer_t {
	int buf_size;
	int cur_pos;
	int rem_size;
	uint8_t *buf;
};

void init_reader_buffer(struct reader_buffer_t *b, int size)
{
	b->buf_size = size;
	b->cur_pos = 0;
	b->rem_size = 0;
	b->buf = malloc(size);
}

void destroy_reader_buffer(struct reader_buffer_t *b)
{
	free(b->buf);
	b->buf = NULL;
	b->buf_size = 0;
}

static inline void KMC_add_prefix_kmer(uint8_t *kmer, int k, int len,
								uint64_t prefix)
{
	int l = (len + 3) / 4, i;
	for (i = 0; i < l; ++i) {
		kmer[k++] = prefix & 0xff;
		prefix >>= 8;
	}
}

static inline void fill_kmer_char(uint8_t k, char *b, int len)
{
	while (len) {
		b[--len] = nt4_char[(uint32_t)(k & 0x3)];
		k >>= 2;
	}
}

void KMC_retrive_kmer(const char *path, struct kmc_info_t *data, void *bundle,
				void (*process)(uint8_t *, uint32_t, void *))
{
	FILE *fp = fopen(path, "rb");
	char sig[4];
	/* read signature to check whether truncated file */
	fseek(fp, 0, SEEK_SET);
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, KMCS, 4))
		__ERROR("Wrong KMC suffix format file");
	fseek(fp, -4, SEEK_END);
	xfread(sig, 4, 1, fp);
	if (strncmp(sig, KMCS, 4))
		__ERROR("[Open] Truncated KMC suffix file");

	/* Calculate record size */
	fseek(fp, 0, SEEK_END);
	uint64_t file_size = ftell(fp) - 8;
	if (file_size % data->header.total_kmers != 0)
		__ERROR("Total kmer does not divide kmer record field total size");
	int record_size = file_size / data->header.total_kmers;
	int suffix_size = (data->header.kmer_length - data->header.lut_prefix_length) / 4;
	if (record_size != suffix_size + (int)data->header.counter_size)
		__ERROR("KMC record size not consistent (%d != %d + %u)",
			record_size, suffix_size, data->header.counter_size);

	/* iterate through buffer */
	uint64_t prefix, prefix_mask, n_kmers, k, i;
	prefix_mask = ((uint64_t)1 << (2 * data->header.lut_prefix_length)) - 1;
	i = 0;
	struct reader_buffer_t b;
	init_reader_buffer(&b, 1 << 24);
	fseek(fp, 4, SEEK_SET);
	b.rem_size = fread(b.buf, 1, b.buf_size, fp);
	b.cur_pos = 0;
	uint8_t *kmer = alloca((data->header.kmer_length + 3) / 4);
	while (data->prefix_offset[i + 1] != data->header.total_kmers + 1) {
		prefix = i & prefix_mask;
		n_kmers = data->prefix_offset[i + 1] - data->prefix_offset[i];
		for (k = 0; k < n_kmers; ++k) {
			/* fill buffer */
			if (b.rem_size < record_size) {
				memmove(b.buf, b.buf + b.cur_pos, b.rem_size);
				b.rem_size += fread(b.buf + b.rem_size, 1,
						b.buf_size - b.rem_size, fp);
				b.cur_pos = 0;
			}
			if (b.rem_size < record_size)
				__ERROR("Truncated KMC suffix file");
			int j;
			for (j = 0; j < suffix_size; ++j)
				kmer[suffix_size - j - 1] = b.buf[b.cur_pos + j];
			KMC_add_prefix_kmer(kmer, suffix_size,
					data->header.lut_prefix_length, prefix);
			uint32_t counter = *((uint32_t *)(b.buf + b.cur_pos + suffix_size));
			process(kmer, counter, bundle);
			b.cur_pos += record_size;
			b.rem_size -= record_size;
		}
		++i;
	}
	destroy_reader_buffer(&b);
	fclose(fp);
}

