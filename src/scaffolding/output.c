#include <stdio.h>
#include "attribute.h"
#include "alloca.h"
#include "utils.h"
#include <string.h>
#include "assembly_graph.h"
#include <stdlib.h>

gint_t dump_edge_seq_reduce_N(char **seq, uint32_t *m_seq, struct asm_edge_t *e)
{
	for (uint32_t i = 0; i < e->n_holes; i++) {
		if (e->l_holes[i] > 1000)
			e->l_holes[i] = 1000;
	}
	return dump_edge_seq(seq, m_seq, e);
}

void print_seq(FILE *fp, int index, char *seq, int len, int cov)
{
	fprintf(fp, ">SEQ_%lld_length_%lld_count_%llu\n", (long long)index,
		(long long)len, (long long unsigned)cov);
	gint_t k = 0;
	char *buf = alloca(81);
	while (k < len) {
		gint_t l = MIN(80, len - k);
		memcpy(buf, seq + k, l);
		buf[l] = '\0';
		fprintf(fp, "%s\n", buf);
		k += l;
	}
	while (k < len) {
		gint_t l = __min(80, len - k);
		memcpy(buf, seq + k, l);
		buf[l] = '\0';
		fprintf(fp, "%s\n", buf);
		k += l;
	}
}

void print_contig(struct asm_graph_t *g, FILE *out_file, int index, int n_contig, int *list_contig)
{
	char *seq = NULL, *total_seq = NULL, *NNN = NULL;
	const int len_NNN = 300;
	NNN = calloc(len_NNN, sizeof(char));
	for (int i = 0; i < len_NNN; i++) 
		NNN[i] = 'N';
	uint32_t seq_len = 0;
 	int total_len = 0;
	for(int i = 0; i < n_contig; i++) {
		int e = list_contig[i];
		int len_of_contig = dump_edge_seq_reduce_N(&seq, &seq_len, &g->edges[e]);
		total_seq = realloc(total_seq, (total_len + len_of_contig + len_NNN) * sizeof(char));
		memcpy(total_seq + total_len, seq, len_of_contig);
		memcpy(total_seq + total_len + len_of_contig, NNN, len_NNN);
		total_len += len_of_contig + len_NNN;
	}
	print_seq(out_file, index, total_seq, total_len, 1);
	for(int i = 0; i < n_contig; i++) {
		int e = list_contig[i];
	}
	free(seq);
	free(total_seq);
	free(NNN);
}

