#ifndef __BARCODE_BUILDER__
#define __BARCODE_BUILDER__
#include "../include/bwa.h"
#include "../include/bwamem.h"
#include "cluster_molecules.h"
struct asm_align_t {
	int rid;
	int pos;
	int score:30, strand:2;
	int aligned;
};

int ksw_global2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int *n_cigar_, uint32_t **cigar_);
mem_opt_t *asm_memopt_init();
struct asm_align_t asm_reg2aln(const mem_opt_t *opt, const bntseq_t *bns,
	const uint8_t *pac, int l_query, const uint8_t *query, const mem_alnreg_t *ar);
struct read_path_t parse_read_path_from_opt(struct opt_proc_t *opt);

void *rp_count_buffer_iterator(void *data);

void get_all_read_pairs_count(struct opt_proc_t *opt, khash_t(long_int) *rp_count);

struct rp_count_bundle_t{
	struct asm_graph_t *g;
	struct dqueue_t *q;
	bwaidx_t *bwa_idx;
	mem_opt_t *bwa_opt;
	pthread_mutex_t *lock;
	khash_t(int_int) **rp_count;
	pthread_mutex_t *lock_rp;
};
#endif
