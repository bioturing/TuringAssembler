#include "cluster_molecules.h"
#include "helper.h"
#define MAX_RADIUS 40000

int get_shortest_path(struct asm_graph_t *g, int source, int target)
{
	struct queue_t *q = calloc(1, sizeof(struct queue_t));
	struct dijkstra_node_t wrapper = {
		.vertex = source,
		.len = g->edges[source].seq_len
	};
	push_queue(q, pointerize(&wrapper, sizeof(struct dijkstra_node_t)));

	khash_t(int_int) *L = kh_init(int_int);
	put_in_map(L, source, wrapper.len);
	while (!is_queue_empty(q)){
		int p = q->front;
		for (int i = q->front + 1; i < q->back; ++i){
			if (((struct dijkstra_node_t *)q->data[p])->len >
				((struct dijkstra_node_t *)q->data[i])->len)
				p = i;
		}
		if (p != q->front){
			struct dijkstra_node_t *tmp = q->data[q->front];
			q->data[q->front] = q->data[p];
			q->data[p] = tmp;
			/*swap(q->data + q->front, q->data + p,
					sizeof(struct dijkstra_node_t *));*/
		}
		struct dijkstra_node_t *node = get_queue(q);
		pop_queue(q);
		int v = node->vertex;
		int len = node->len;
		if (get_in_map(L, v) != len)
			goto dijkstra_node_outdated;
		if (v == target)
			break;
		if (len > MAX_RADIUS)
			break;
		int tg = g->edges[v].target;
		for (int i = 0; i < g->nodes[tg].deg; ++i){
			int u = g->nodes[tg].adj[i];
			if (check_in_map(L, u) == 0)
				put_in_map(L, u, 2e9);
			if ((uint32_t) get_in_map(L, u) > len + g->edges[u].seq_len){
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
	int res = -1;
	if (check_in_map(L, target) != 0)
		res = get_in_map(L, target);
	kh_destroy(int_int, L);
	free_queue_content(q);
	destroy_queue(q);

	return res;
}

void get_sub_graph(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct mm_hits_t *hits)
{
	khash_t(set_int) *edges = kh_init(set_int);
	for (khiter_t it = kh_begin(hits->edges); it != kh_end(hits->edges); ++it){
		if (!kh_exist(hits->edges, it))
			continue;
		int e = kh_key(hits->edges, it);
		put_in_set(edges, e);
		put_in_set(edges, g->edges[e].rc_id);
	}

	int *tmp = calloc(kh_size(edges), sizeof(int));
	int n = 0;
	for (khiter_t it = kh_begin(edges); it != kh_end(edges); ++it){
		if (!kh_exist(edges, it))
			continue;
		tmp[n++] = kh_key(edges, it);
	}
	kh_destroy(set_int, edges);

	FILE *f = fopen(opt->lc, "w");
	fprintf(f, "graph %s{\n", opt->bx_str);
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j){
			if (i == j)
				continue;
			int len = get_shortest_path(g, tmp[i], tmp[j]);
			if (len != -1 && len <= MAX_RADIUS)
				fprintf(f, "%d -> %d;\n", tmp[i], tmp[j]);
		}
	}
	fprintf(f, "}\n");
	fclose(f);
	free(tmp);
}

