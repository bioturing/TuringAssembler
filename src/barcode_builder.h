#ifndef __BARCODE_BUILDER__
#define __BARCODE_BUILDER__
struct asm_align_t {
	int rid;
	int pos;
	int score:30, strand:2;
	int aligned;
};

mem_opt_t *asm_memopt_init();
struct asm_align_t asm_reg2aln(const mem_opt_t *opt, const bntseq_t *bns,
	const uint8_t *pac, int l_query, const uint8_t *query, const mem_alnreg_t *ar);
#endif
