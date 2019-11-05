#ifndef __BARCODE_RESOLVE2_H__
#define __BARCODE_RESOLVE2_H__
#include "assembly_graph.h"
KHASH_SET_INIT_INT64(gint);

struct opt_local_t {
	char *out_dir;
	struct read_path_t *read_path;
	khash_t(bcpos) *dict;
	int ksize;
	int n_threads;
	int mmem;
};

struct result_local_t {
	uint32_t *seq;
	int len;
	int trim_e1;
	int trim_e2;
};
int get_reads_local_graph(struct read_path_t *reads, struct read_path_t *rpath,
			khash_t(bcpos) *dict, struct asm_graph_t *g,
			gint_t e1, gint_t e2, const char *prefix);
int get_reads_kmer_check(struct opt_proc_t *opt, struct asm_graph_t *g,
		int e1, int e2, struct read_path_t *local_read_path);
void get_shared_reads(struct read_path_t *reads, struct read_path_t *rpath,
		khash_t(bcpos) *dict, struct asm_graph_t *g, gint_t e1,
		gint_t e2, const char *prefix, int contig_level);
void get_union_reads(struct read_path_t *reads, struct read_path_t *rpath,
		khash_t(bcpos) *dict, struct asm_graph_t *g, gint_t e1,
		gint_t e2, const char *prefix, int contig_level);
int check_large_pair_superior(struct asm_graph_t *g, gint_t e1,
							gint_t e2, gint_t e2a);
void construct_read_index(struct read_path_t *rpath, khash_t(bcpos) *h);
khash_t(gint) *get_shared_bc(khash_t(gint) *h1, khash_t(gint) *h2);
khash_t(gint) *get_union_bc(khash_t(gint) *h1, khash_t(gint) *h2);
khash_t(gint) *get_exclude_bc(khash_t(gint) *h_in, khash_t(gint) *h_ex);
khash_t(gint) *barcode_hash_2_khash(struct barcode_hash_t *bc);
void khash_2_arr(khash_t(gint) *h, uint64_t **arr, int *n);
int fill_path_local(struct opt_local_t *opt, struct asm_graph_t *g0, struct asm_graph_t *g,
			gint_t e1, gint_t e2, struct result_local_t *sret);
void resolve_1_2(struct asm_graph_t *g, struct opt_proc_t *opt);
#endif
