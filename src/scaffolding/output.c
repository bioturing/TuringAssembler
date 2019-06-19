#include <stdio.h>
#include "attribute.h"
#include "alloca.h"
#include "utils.h"
#include <string.h>
#include "assembly_graph.h"
#include <stdlib.h>
#include "scaffolding/edge.h"
#include "scaffolding/contig.h"

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
	total_len -= len_NNN;
	print_seq(out_file, index, total_seq, total_len, 1);
	for(int i = 0; i < n_contig; i++) {
		int e = list_contig[i];
	}
	free(seq);
	free(total_seq);
	free(NNN);
}

void print_gfa_from_E(struct asm_graph_t *g, int n_e, struct scaffold_edge *listE, int n_v, int *listV, FILE *out_graph)
{
	struct scaffold_edge *list_one_dir_E = calloc(n_e, sizeof(struct scaffold_edge));
	for (int i = 0; i < n_e; i++) {
		list_one_dir_E[i] = listE[i];
		normalize_min_index(g, list_one_dir_E+i);
	}
	for (int i = 0; i < n_v; i++) {
		struct asm_edge_t *e = &g->edges[listV[i]];
		char *seq = NULL;
		uint32_t seq_len = 0;
		dump_edge_seq_reduce_N(&seq, &seq_len, e);
		fprintf(out_graph,"S\t%d\t%s\tKC:i:%lu\n", listV[i], seq, e->count);
	}

	for (int i = 0; i < n_e; i++) {
		fprintf(out_graph, "L\t%d\t%c\t%d\t%c\t45M\n", 
			list_one_dir_E[i].src, list_one_dir_E[i].rv_src == 0?'+':'-', list_one_dir_E[i].des, list_one_dir_E[i].rv_des == 0?'+':'-');
		fprintf(out_graph, "L\t%d\t%c\t%d\t%c\t45M\n", 
			list_one_dir_E[i].des, list_one_dir_E[i].rv_des == 0?'-':'+', list_one_dir_E[i].src, list_one_dir_E[i].rv_src == 0?'-':'+');
	}
}

