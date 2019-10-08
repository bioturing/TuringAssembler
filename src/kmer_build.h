#ifndef __KMER_BUILD_H__
#define __KMER_BUILD_H__
#include"assembly_graph.h"

void save_graph_info(const char *out_dir, struct asm_graph_t *g, const char *suffix);
void build_graph_from_scratch(int ksize, int n_threads, int mmem, int n_files,
				char **files_1, char **files_2, char *work_dir,
						struct asm_graph_t *g);
void build_graph_from_scratch_without_count(int ksize, int n_threads, int mmem, int n_files,
											char **files_1, char **files_2, char *work_dir,
											struct asm_graph_t *g);
void build_0_1(struct asm_graph_t *g0, struct asm_graph_t *g);

#endif /* __KMER_BUILD_H__ */
