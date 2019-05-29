#ifndef __KMC_READER_H__
#define __KMC_READER_H__

#include "attribute.h"
#include "assembly_graph.h"

struct kmc_header_t {
	uint32_t kmer_length;
	uint32_t mode;
	uint32_t counter_size;
	uint32_t lut_prefix_length;
	uint32_t signature_length;
	uint32_t min_count;
	uint32_t max_count;
	uint64_t total_kmers;
	uint8_t  both_strands;
	uint8_t  tmp_char[3];
	uint32_t tmp_uint[6];
	uint32_t KMC_VER;
} __attribute__((packed));

struct kmc_info_t {
	struct kmc_header_t header;
	uint64_t sig_map_size;
	uint64_t prefix_offset_size;
	uint64_t *prefix_offset;
	uint32_t *sig_map;
};

void KMC_retrive_kmer(const char *path, struct kmc_info_t *data, void *bundle,
				void (*process)(uint8_t *, uint32_t, void *));

void KMC_read_prefix(const char *path, struct kmc_info_t *data);

void destroy_kmc_info(struct kmc_info_t *inf);

#endif /* __KMC_GRAPH_H__ */
