#include "contig_graph.h"
#include <stdlib.h>
#include"verbose.h"
#include"global_params.h"
#include"scaffolding/buck.h"

void reverse_arr(int n_arr, int *arr)
{
	for (int i = 0; i < n_arr/2; i++) {
		int tmp = arr[i];
		arr[i] = arr[n_arr-1-i];
		arr[n_arr-1-i] = tmp;
	}
}

struct contig_graph* build_graph_contig(struct asm_graph_t *g, int *mark_short)
{
	struct contig_graph *a = calloc(1, sizeof(struct contig_graph));
	int n = g->n_v, m = 0;
	a->head = calloc(g->n_v, sizeof(int));
	a->next = calloc(g->n_e * 2, sizeof(int));
	a->edges = calloc(g->n_e * 2, sizeof(int));
	a->edge_index = calloc(g->n_e * 2, sizeof(int));
	for (int i = 0 ; i < n; i++){
		a->head[i] = -1;
	}
	for (int i = 0; i < g->n_e; i++) {
		struct asm_edge_t *e = &(g->edges[i]);
		int x = e->source, y = e->target;
		VERBOSE_FLAG(0, "short edge %d src %d des %d\n", i, x, y);
	}
	for (int i = 0; i < g->n_e; i++) if (mark_short[i]){
		struct asm_edge_t *e = &(g->edges[i]);
		int x = e->source, y = e->target;
		a->edges[m] = y;
		a->next[m] = a->head[x];
		a->head[x] = m;
		a->edge_index[m] = i;
		m++;
	}
	VERBOSE_FLAG(0, "m short %d\n", m);
	a->n_v = g->n_v;
	a->n_e = g->n_e;
	return a;
}

void push_queue(int *n_queue, int *queue, int value)
{
	queue[*n_queue] = value;
	(*n_queue)++;
}

int pop_queue(int *queue, int *bot) {
	(*bot)++;
	return queue[*bot-1];
}

void find_path_short(int src, int des, struct asm_graph_t *g, struct contig_graph *contig_g,
		int *n_insert_path, int **insert_path, int *mark_short, 
		struct barcode_hash_t *hl, struct barcode_hash_t *hr)
{
	int limit_queue = 5000;
	VERBOSE_FLAG(0, "find src %d des %d\n", src, des);
	int *queue = calloc(limit_queue, sizeof(int)), n_queue = 0, bot = 0;
	int *trace = calloc(limit_queue, sizeof(int));
	int *edge_index = calloc(limit_queue, sizeof(int));
	push_queue(&n_queue, queue, src);
	while (bot < n_queue && n_queue < limit_queue){
		int x = pop_queue(queue, &bot); 
		VERBOSE_FLAG(0, "x %d\n", x);
		if (x == des) {
			break;
		}
		for (int i = contig_g->head[x]; i != -1; i = contig_g->next[i]) {
			int y =  contig_g->edges[i];
			int e_index = contig_g->edge_index[i];
			VERBOSE_FLAG(0, "to edge %d\n", e_index);
// can not check barcode because this contig is too short
//			if (!check_share_bc(g, hl, hr, e_index))
//				continue;
			if (mark_short[e_index]){
				VERBOSE_FLAG(0, "e index %d\n" , e_index);
				push_queue(&n_queue, queue, y);
				edge_index[n_queue-1] = e_index;
				trace[n_queue-1] = bot-1;
				mark_short[e_index]--;
				int rc = g->edges[e_index].rc_id;
				mark_short[rc]--;
			}
			if (n_queue == limit_queue)
				break;
		}
	}
	for (int i = 0; i < n_queue; i++){
		int e = edge_index[i];
		mark_short[e]++;
		int rc = g->edges[e].rc_id;
		mark_short[rc]++;
	}
	int pos = bot-1;
	for (pos = n_queue; pos >= 0; pos--) {
		if (queue[pos] == des){
			break;
		}
	}
	if (pos < 0 || queue[pos] != des){
		return;
		free(queue);
		free(trace);
		free(edge_index);
	}
	VERBOSE_FLAG(1, "found short path\n");
	while (pos != 0) {
		(*insert_path) = realloc(*insert_path, (*n_insert_path + 1) * sizeof(int));
		(*insert_path)[(*n_insert_path)++] = edge_index[pos];
		pos = trace[pos];
	}
	reverse_arr(*n_insert_path, *insert_path); 
	VERBOSE_FLAG(0, "n insert %d\n", *n_insert_path);
	for (int i = 0 ; i < *n_insert_path; i++)
		VERBOSE_FLAG(0, "insert path %d ", (*insert_path)[i]);
	VERBOSE_FLAG(0, "\n");
	for (int i = 0; i < *n_insert_path; i++){
		int s = (*insert_path)[i];
		mark_short[s]--;
		int rc = g->edges[s].rc_id;
		mark_short[rc]--;
		VERBOSE_FLAG(0, "insert path %d rc %d", s, rc);
	}
	free(queue);
	free(trace);
	free(edge_index);
}

