#include "cluster_molecules.h"
#include "helper.h"

int get_shortest_path(struct asm_graph_t *g, int source, int target)
{
	struct queue_t *q = calloc(1, sizeof(struct queue_t));
	struct dijkstra_node_t wrapper = {
		.vertex = source,
		.len = g->edges[source].seq_len
	};
	push_queue(q, pointerize(&wrapper, sizeof(struct dijkstra_node_t)));

	khash_t(int_int) *L = kh_init(int_int);
	put_in_map(L, source, 0);
	while (!is_queue_empty(q)){
		int p = q->front;
		for (int i = q->front + 1; i < q->back; ++i){
			if (((struct dijkstra_node_t *)q->data[p])->len >
				((struct dijkstra_node_t *)q->data[i])->len)
				p = i;
		}
		if (p != q->front)
			swap(q->data[q->front], q->data[p], sizeof(struct dijkstra_node_t));
		struct dijkstra_node_t *node = get_queue(q);
		pop_queue(q);
		int v = node->vertex;
		int len = node->len;
		if (get_in_map(L, v) != len)
			goto dijkstra_node_outdated;
		if (v == target)
			break;
		int tg = g->edges[v].target;
		for (int i = 0; i < g->nodes[tg].deg; ++i){
			int u = g->nodes[tg].adj[i];
			if (check_in_map(L, u) == 0)
				put_in_map(L, u, 2e9);
			if (get_in_map(L, u) > len + g->edges[u].seq_len){
				khiter_t it = kh_get(int_int, L, u);
				kh_val(L, it) = len + g->edges[u].seq_len;
				wrapper.vertex = u;
				wrapper.len = len + g->edges[u].seq_len;
				push_queue(q, pointerize(&wrapper,
					sizeof(struct dijkstra_node_t)));
			}
		}
dijkstra_node_outdated:
		free(node);
	}
	kh_destroy(int_int, L);
	free_queue_content(q);
	destroy_queue(q);
}

void get_sub_graph_of_molecules(struct asm_graph_t *g, struct mm_hits_t *hits)
{
	khash_t(set_int) edges = kh_init(set_int);
	//for (khiter_t it = kh_begin(hits); it != 
}

