#ifndef __KMER_BUILD_H__
#define __KMER_BUILD_H__
#include"assembly_graph.h"
#include "khash.h"
#include "../include/bwa.h"
#include "../include/bwamem.h"

struct kmbuild_bundle_t {
	struct kmhash_t *h;
	int ksize;
	uint8_t *k1;
	uint8_t *k2;
	uint8_t *k1_rc;
	uint8_t *k2_rc;
};

void build_graph_from_scratch(int ksize, int n_threads, int mmem, int n_files,
				char **files_1, char **files_2, char *work_dir,
						struct asm_graph_t *g);
void build_graph_from_scratch_without_count(int ksize, int n_threads, int mmem, int n_files,
											char **files_1, char **files_2, char *work_dir,
											struct asm_graph_t *g);
void kmbuild_bundle_init(struct kmbuild_bundle_t *b, struct kmhash_t *h,
						 int ksize);
void kmbuild_bundle_destroy(struct kmbuild_bundle_t *b);
void build_local_graph_cov(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *lg, int e1, int e2, char *work_dir,
		char *R1_path, char *R2_path);

#endif /* __KMER_BUILD_H__ */
