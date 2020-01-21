//
// Created by che on 08/01/2020.
//

#ifndef SKIPPING_MULTI_KMERS_H
#define SKIPPING_MULTI_KMERS_H

#include "assembly_graph.h"
void build_graph_from_scratch(int ksize, int n_threads, int mmem, int n_files,
			      char **files_1, char **files_2, char *work_dir,
			      struct asm_graph_t *g);
void load_asm_graph(struct asm_graph_t *g, const char *path);
void resolve_multi_kmer(struct opt_proc_t *opt, struct asm_graph_t *g, int lastk);
void save_graph_info(const char *out_dir, struct asm_graph_t *g, const char *suffix);
void resolve_graph_operation(struct asm_graph_t *g0, struct asm_graph_t *g);
void init_logger(int level, const char * file_path);
void set_log_stage(char *stage);
gint_t remove_tips(struct asm_graph_t *g);
void asm_condense(struct asm_graph_t *g0, struct asm_graph_t *g);
void asm_graph_destroy(struct asm_graph_t *g);
gint_t resolve_simple_bubble(struct asm_graph_t *g);
void dirty(struct asm_graph_t *g, struct opt_proc_t *opt);
#endif //SKIPPING_MULTI_KMERS_H
