//
// Created by BioTuring on 2019-11-24.
//

#include <stdlib.h>
#include <stdio.h>

#include "heap/heap.h"
#include "assembly_graph.h"
#include "log.h"

#define DEFAULT_S_SIZE 1024
#define MINIMUM_DISTANCE 10000

struct dijkstra_t {
    struct kheap_t *h; /* Prioritized queue to store the shortest distance of remain nodes */
    uint32_t *S;       /* Set of edges that has been explored*/
    uint64_t n;        /* Number of node in S */
    uint64_t size;        /* Allocated room for S */

};

/**
 * @brief Init the struct for dijkstra algorithm
 * @param g main assembly graph G
 * @param e the very first element that we want to explore
 * @return the struct
 */
struct dijkstra_t *init_dijkstra(struct asm_graph_t *g)
{
	struct dijkstra_t *d = malloc(sizeof(struct dijkstra_t));
	d->h = kheap_init(g->n_e);
	for (int i = 0; i < g->n_e; ++i)
		d->h->key[i] = UINT32_MAX;
	for (int i = 0; i < g->n_e; ++i) {
		kheap_insert(d->h, i);
		d->h->pos[i] = d->h->n;
	}
	d->S = calloc(DEFAULT_S_SIZE, sizeof(uint32_t));
	d->n = 0;
	d->size = DEFAULT_S_SIZE;
	return d;
}

void add_node(struct dijkstra_t *d, struct asm_graph_t *g, gint_t e)
{
	gint_t dst = g->edges[e].target;
	gint_t u;
	int i;
	if (d->n + 1 > d->size) {
		d->size <<= 1;
		d->S = realloc(d->S, sizeof(uint32_t) * d->size);
	}
	d->S[d->n++] = e;
	for (i = 0; i < g->nodes[dst].deg; ++i) {
		u = g->nodes[dst].adj[i];
		if (d->h->key[u] > d->h->key[e] + g->edges[u].seq_len) {
			d->h->key[u] = d->h->key[e] + g->edges[u].seq_len;
			heapify_up(d->h, d->h->pos[u]);
		}
	}
}

void find_shortest_path(struct asm_graph_t *g, gint_t e)
{
	struct dijkstra_t *d = init_dijkstra(g);
	d->h->key[e] = 0;
	gint_t v = e;
	while (1) {
		add_node(d, g, v);
		log_info("Added v=%llu, d=%d", v, d->h->key[v]);
		v = d->h->H[1]; //1-based
		if (d->h->key[v] > MINIMUM_DISTANCE && d->h->key[v] != UINT32_MAX) {
			break;
		}
		kheap_delete(d->h, 0);
	}
}





