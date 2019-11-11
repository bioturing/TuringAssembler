#ifndef __KMER_BUILD_H__
#define __KMER_BUILD_H__
#include"assembly_graph.h"

struct kmbuild_bundle_t {
	struct kmhash_t *h;
	int ksize;
	uint8_t *k1;
	uint8_t *k2;
	uint8_t *k1_rc;
	uint8_t *k2_rc;
};

void save_graph_info(const char *out_dir, struct asm_graph_t *g, const char *suffix);
void build_graph_from_scratch(int ksize, int n_threads, int mmem, int n_files,
				char **files_1, char **files_2, char *work_dir,
						struct asm_graph_t *g);
void build_graph_from_scratch_without_count(int ksize, int n_threads, int mmem, int n_files,
											char **files_1, char **files_2, char *work_dir,
											struct asm_graph_t *g);
void build_0_1(struct asm_graph_t *g0, struct asm_graph_t *g);

void kmbuild_bundle_init(struct kmbuild_bundle_t *b, struct kmhash_t *h,
						 int ksize);
void kmbuild_bundle_destroy(struct kmbuild_bundle_t *b);
#endif /* __KMER_BUILD_H__ */
