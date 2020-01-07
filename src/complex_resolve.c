#include "complex_resolve.h"
#include "khash.h"
#include "log.h"
#include "resolve.h"
#include "verbose.h"
#include "helper.h"
#include "dqueue.h"
#include "process.h"
#include "utils.h"


static inline void asm_add_node_adj(struct asm_graph_t *g, gint_t u, gint_t e)
{
	g->nodes[u].adj = realloc(g->nodes[u].adj, (g->nodes[u].deg + 1) * sizeof(gint_t));
	g->nodes[u].adj[g->nodes[u].deg++] = e;
}

void init_resolve_bulges(struct asm_graph_t *g, struct resolve_bulges_bundle_t *bundle)
{
	bundle->graph = g;

	bundle->B_vertices = calloc(1, sizeof(struct queue_t));
	init_queue(bundle->B_vertices, bundle->graph->n_v);

	bundle->dom_vertices = calloc(1, sizeof(struct queue_t));
	init_queue(bundle->dom_vertices, bundle->graph->n_v);

	bundle->closest = calloc(1, sizeof(struct queue_t));
	init_queue(bundle->closest, bundle->graph->n_v);

	bundle->source = -1;
	bundle->dom = kh_init(set_int);
	bundle->B = kh_init(set_int);
	bundle->PE = kh_init(int_int);
	bundle->L = kh_init(int_int);
}

void reset_source(struct resolve_bulges_bundle_t *bundle, int s)
{
	bundle->B_vertices->front = bundle->B_vertices->back = 0;
	bundle->dom_vertices->front = bundle->dom_vertices->back = 0;
	bundle->closest->front = bundle->closest->back = 0;
	bundle->source = s;
	kh_destroy(set_int, bundle->dom);
	bundle->dom = kh_init(set_int);

	kh_destroy(set_int, bundle->B);
	bundle->B = kh_init(set_int);

	kh_destroy(int_int, bundle->PE);
	bundle->PE = kh_init(int_int);

	kh_destroy(int_int, bundle->L);
	bundle->L = kh_init(int_int);
}

void bulges_bundle_destroy(struct resolve_bulges_bundle_t *bundle)
{
	destroy_queue(bundle->B_vertices);
	free(bundle->B_vertices);
	destroy_queue(bundle->dom_vertices);
	free(bundle->dom_vertices);
	destroy_queue(bundle->closest);
	free(bundle->closest);

	kh_destroy(set_int, bundle->dom);
	kh_destroy(set_int, bundle->B);
	kh_destroy(int_int, bundle->PE);
	kh_destroy(int_int, bundle->L);
}

void get_dominated_vertices(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	int s = bundle->source;
	khash_t(set_int) *dom = bundle->dom;

	int s_rc = graph->nodes[s].rc_id;
	khash_t(set_int) *s_parents = kh_init(set_int);
	for (int i = 0; i < graph->nodes[s_rc].deg; ++i){
		int e = graph->edges[graph->nodes[s_rc].adj[i]].rc_id;
		int p = graph->edges[e].source;
		put_in_set(s_parents, p);
	}

	struct queue_t *queue = calloc(1, sizeof(struct queue_t));
	init_queue(queue, 1024);
	push_queue(queue, pointerize(&s, sizeof(int)));

	khash_t(int_int) *deg_in = kh_init(int_int);
	while (is_queue_empty(queue) == 0){
		int *v = get_queue(queue);
		pop_queue(queue);
		put_in_set(dom, *v);
		push_queue(bundle->dom_vertices, pointerize(v, sizeof(int)));
		for (int i = 0; i < graph->nodes[*v].deg; ++i){
			int u = get_adj_node(graph, *v, i);
			if (check_in_map(deg_in, u) == 0)
				put_in_map(deg_in, u, 0);
			increase_in_map(deg_in, u, 1);
			int used_deg = get_in_map(deg_in, u);
			int u_rc = graph->nodes[u].rc_id;
			if (used_deg == graph->nodes[u_rc].deg &&
				check_in_set(s_parents, u) == 0)
				push_queue(queue, pointerize(&u, sizeof(int)));
		}
		free(v);
	}
	kh_destroy(int_int, deg_in);
	destroy_queue(queue);
	free(queue);
	kh_destroy(set_int, s_parents);
}

void add_vertex_to_B(struct resolve_bulges_bundle_t *bundle, int v)
{
	//__VERBOSE("Add %d\n", v);
	put_in_set(bundle->B, v);
	push_queue(bundle->B_vertices, pointerize(&v, sizeof(int)));
}

void add_vertex_to_B_dfs(struct resolve_bulges_bundle_t *bundle, int v,
		khash_t(set_int) *in_queue, struct queue_t *q, int depth)
{
	struct asm_graph_t *graph = bundle->graph;
	int int_vertex = 0;
	if (depth == 0){
		for (int i = 0; i < graph->nodes[v].deg; ++i){
			int u = get_adj_node(graph, v, i);
			if (check_in_set(bundle->B, u) != 0)
				int_vertex = 1;
		}
	} else {
		int_vertex = 1;
	}
	if (int_vertex && check_in_set(in_queue, v) == 0){
		put_in_set(in_queue, v);
		push_queue(q, pointerize(&v, sizeof(int)));
		//__VERBOSE("queue in %d\n", v);
	}
	if (check_in_set(bundle->B, v) != 0)
		return;
	add_vertex_to_B(bundle, v);
	//__VERBOSE("add %d\n", v);
	int v_rc = graph->nodes[v].rc_id;
	for (int i = 0; i < graph->nodes[v_rc].deg; ++i){
		int pe = graph->edges[graph->nodes[v_rc].adj[i]].rc_id;
		int p = graph->edges[pe].source;
		add_vertex_to_B_dfs(bundle, p, in_queue, q, depth + 1);
	}
}

int get_closure(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	int s = bundle->source;
	int res = 1;
	struct queue_t q;
	init_queue(&q, 1024);

	khash_t(set_int) *in_queue = kh_init(set_int);
	struct queue_t *B_vertices = bundle->B_vertices;
	/*__VERBOSE("before: ");
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		__VERBOSE("%d ", *v);
	}
	__VERBOSE("\n");*/
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		//__VERBOSE("%d %d\n", i, vg->vertices[i].deg_out);
		for (int i = 0; i < graph->nodes[*v].deg; ++i){
			int u = get_adj_node(graph, *v, i);
			//__VERBOSE("%d %d\n", i, v);
			if (check_in_set(bundle->B, u) != 0){
				put_in_set(in_queue, *v);
				push_queue(&q, pointerize(v, sizeof(int)));
				//__VERBOSE("queue in %d\n", *v);
				break;
			}
		}
	}

	while (res && is_queue_empty(&q) == 0){
		int *v = get_queue(&q);
		pop_queue(&q);
		//__VERBOSE("queue out %d\n", *v);
		for (int i = 0; i < graph->nodes[*v].deg; ++i){
			int u = get_adj_node(graph, *v, i);
			//__VERBOSE("hihi %d %d\n", *v, u);
			if (check_in_set(bundle->dom, u) == 0){
				res = 0;
				break;
			}
			if (check_in_set(bundle->B, u) != 0)
				continue;
			add_vertex_to_B_dfs(bundle, u, in_queue, &q, 0);
		}
		free(v);
	}
	free_queue_content(&q);
	kh_destroy(set_int, in_queue);
	destroy_queue(&q);
	return res;
}

void bfs_to_sinks(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	kh_destroy(int_int, bundle->PE);
	bundle->PE = kh_init(int_int);
	khash_t(int_int) *PE = bundle->PE;

	struct queue_t *q = calloc(1, sizeof(struct queue_t));
	init_queue(q, 1024);
	push_queue(q, pointerize(&bundle->source, sizeof(int)));
	khash_t(set_int) *visited = kh_init(set_int);
	put_in_set(visited, bundle->source);
	put_in_map(PE, bundle->source, -1);

	while (!is_queue_empty(q)){
		int *v = get_queue(q);
		pop_queue(q);
		for (int i = 0; i < graph->nodes[*v].deg; ++i){
			int e = graph->nodes[*v].adj[i];
			int u = get_adj_node(graph, *v, i);
			if (check_in_set(bundle->B, u) == 0)
				continue;
			if (check_in_set(visited, u) == 0){
				put_in_set(visited, u);
				put_in_map(PE, u, e);
				push_queue(q, pointerize(&u, sizeof(int)));
			}
		}
		free(v);
	}
	destroy_queue(q);
	free(q);
	kh_destroy(set_int, visited);
}

void get_distance(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	khash_t(int_int) *L = bundle->L;

	struct queue_t *q = calloc(1, sizeof(struct queue_t));
	init_queue(q, 1024);
	push_queue(q, pointerize(&bundle->source, sizeof(int)));
	put_in_map(L, bundle->source, 0);

	while (!is_queue_empty(q)){
		int *v = get_queue(q);
		pop_queue(q);
		push_queue(bundle->closest, pointerize(v, sizeof(int)));
		for (int i = 0; i < graph->nodes[*v].deg; ++i){
			int u = get_adj_node(graph, *v, i);
			if (check_in_set(bundle->dom, u) == 0)
				continue;
			if (check_in_map(L, u) == 0){
				put_in_map(L, u, get_in_map(L, *v) + 1);
				push_queue(q, pointerize(&u, sizeof(int)));
			}
		}
		free(v);
	}
	destroy_queue(q);
	free(q);
}

int is_complex_closure(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	struct queue_t *B_vertices = bundle->B_vertices;
	int res = 0;
	int s = bundle->source;
	for (int i = 0; i < graph->nodes[s].deg; ++i){
		int v = get_adj_node(graph, s, i);
		if (v == s)
			return 1;
	}

	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		if (check_in_set(bundle->B, graph->nodes[*v].rc_id) != 0)
			return 1;
		for (int j = 0; j < graph->nodes[*v].deg; ++j){
			int e = graph->nodes[*v].adj[j];
			int u = graph->edges[e].target;
			if (check_in_set(bundle->B, u) != 0)
				res = MAX(res, (int) graph->edges[e].seq_len);
		}
	}
	return res >= 1000;
}

int is_closure_tree(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	struct queue_t *B_vertices = bundle->B_vertices;
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		//__VERBOSE("HAHA %d\n", *v);
		int C = 0;
		int v_rc = graph->nodes[*v].rc_id;
		for (int i = 0; i < graph->nodes[v_rc].deg; ++i){
			int pe = graph->edges[graph->nodes[v_rc].adj[i]].rc_id;
			int w = graph->edges[pe].source;
			//__VERBOSE("HAHA %d %d\n", *v, w);
			if (check_in_set(bundle->B, w) != 0)
				++C;
		}
		if (C > 1)
			return 0;
	}
	return 1;
}

int get_next_B_candidate(struct resolve_bulges_bundle_t *bundle)
{
	int res = -1;
	struct queue_t *closest = bundle->closest;
	while (res == -1 && !is_queue_empty(closest)){
		int *v = get_queue(closest);
		pop_queue(closest);
		if (check_in_set(bundle->B, *v) == 0)
			res = *v;
		free(v);
	}
	return res;
}

void supress_bulge(struct resolve_bulges_bundle_t *bundle)
{
	struct asm_graph_t *graph = bundle->graph;
	struct queue_t *B_vertices = bundle->B_vertices;
	khash_t(set_int) *mark = kh_init(set_int);
	put_in_set(mark, bundle->source);
	//__VERBOSE("sink: ");
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		int is_sink = 1;
		for (int j = 0; j < graph->nodes[*v].deg; ++j){
			int u = get_adj_node(graph, *v, j);
			if (check_in_set(bundle->B, u) != 0){
				is_sink = 0;
				break;
			}
		}
		if (is_sink){
			//__VERBOSE("%d ", *v);
			int w = *v;
			while (check_in_set(mark, w) == 0){
				put_in_set(mark, w);
				int e = get_in_map(bundle->PE, w);
				w = graph->edges[e].source;
			}
		}
	}
	//__VERBOSE("\n");
	khash_t(set_int) *rm_edges = kh_init(set_int);
	for (int i = B_vertices->front; i < B_vertices->back; ++i){
		int *v = B_vertices->data[i];
		for (int j = 0; j < graph->nodes[*v].deg; ++j){
			int u = get_adj_node(graph, *v, j);
			int e = graph->nodes[*v].adj[j];
			int rc = graph->edges[e].rc_id;
			//__VERBOSE("HEHE %d %d %d %d\n", *v, u, mark[*v], mark[u]);
			if (check_in_set(bundle->B, u) == 0)
				continue;
			//__VERBOSE("HAHA %d %d %d %d\n", *v, u, mark[*v], mark[u]);
			if (check_in_set(mark, *v) == 0
				|| check_in_set(mark, u) == 0
				|| (get_in_map(bundle->PE, u) != e
					&& get_in_map(bundle->PE, u) != rc)){
				put_in_set(rm_edges, e);
				put_in_set(rm_edges, rc);
				//__VERBOSE("remove %d %d\n", e, rc);
			}
		}
	}
	for (khiter_t it = kh_begin(rm_edges); it != kh_end(rm_edges); ++it){
		if (!kh_exist(rm_edges, it))
			continue;
		int e = kh_key(rm_edges, it);
		//__VERBOSE("remove %d\n", e);
		asm_remove_edge(graph, e);
	}
	kh_destroy(set_int, rm_edges);
	kh_destroy(set_int, mark);
}

int resolve_bulges(struct asm_graph_t *g)
{
	int res = 0;
	struct resolve_bulges_bundle_t bundle;
	init_resolve_bulges(g, &bundle);
	struct asm_graph_t *graph = bundle.graph;
	for (int i = 0; i < graph->n_v; ++i){
		reset_source(&bundle, i);
		get_dominated_vertices(&bundle);
		/*__VERBOSE("dom: ");
		for (int i = 0; i < bundle.dom_vertices->back; ++i)
			__VERBOSE("%d ", *(int *)bundle.dom_vertices->data[i]);
		__VERBOSE("\n");*/
		get_distance(&bundle);

		add_vertex_to_B(&bundle, i);
		free(get_queue(bundle.closest));
		pop_queue(bundle.closest);
		while(1){
			int next_cand = get_next_B_candidate(&bundle);
			if (next_cand == -1)
				break;

			add_vertex_to_B(&bundle, next_cand);
			int ret = get_closure(&bundle);
			if (!ret)
				break;
			if (is_complex_closure(&bundle))
				break;
			if (is_closure_tree(&bundle))
				continue;
			bfs_to_sinks(&bundle);
			/*__VERBOSE("B\n");
			for (int j = bundle.B_vertices->front; j < bundle.B_vertices->back; ++j)
				__VERBOSE("%d ", *(int *)bundle.B_vertices->data[j]);
			__VERBOSE("\n");*/
			/*for (int i = 0; i < g->n_v; ++i){
				for (int j = 0; j < g->nodes[i].deg; ++j){
					int e = g->nodes[i].adj[j];
					int sr = g->edges[e].source;
					int tg = g->edges[e].target;
					if (bundle.B[sr] == 0 || bundle.B[tg] == 0){
						int rc = g->edges[e].rc_id;
						asm_remove_edge(g, rc);
						asm_remove_edge(g, e);
					}
				}
			}*/
			supress_bulge(&bundle);
			log_debug("Bulge detected at %d", i);
			++res;
			break;
		}
		free_queue_content(bundle.B_vertices);
		free_queue_content(bundle.dom_vertices);
		free_queue_content(bundle.closest);
	}
	bulges_bundle_destroy(&bundle);
	test_asm_graph(g);
	return res;
}

int asm_resolve_complex_bulges_ite(struct asm_graph_t *g)
{
	int ite = 0;
	int res = 0;
	do{
		int resolved = resolve_bulges(g);
		if (!resolved)
			break;
		res += resolved;
		++ite;
		log_debug("%d-th iteration: %d complex bulge(s) resolved", ite, resolved);
		struct asm_graph_t g1;
		asm_condense(g, &g1);
		asm_graph_destroy(g);
		*g = g1;

	} while(1);
	log_info("%d complex bulge(s) resolved after %d iterations", res, ite);
	return res;
}

int get_adj_node(struct asm_graph_t *g, int v, int id)
{
	int e = g->nodes[v].adj[id];
	return g->edges[e].target;
}

void create_super_nodes(struct asm_graph_t *g, int e, struct asm_graph_t *supg,
		khash_t(long_int) *node_map)
{
	int pu = g->edges[e].source;
	int pv = g->edges[e].target;
	if (g->edges[e].seq_len > g->ksize + 1){
		supg->nodes = realloc(supg->nodes, (supg->n_v + 2)
				* sizeof(struct asm_node_t));
		memset(supg->nodes + supg->n_v, 0, sizeof(struct asm_node_t) * 2);
		supg->edges = realloc(supg->edges, (supg->n_e + 1)
				* sizeof(struct asm_edge_t));
		int u = supg->n_v;
		int v = supg->n_v + 1;
		asm_add_node_adj(supg, u, supg->n_e);
		asm_clone_seq(supg->edges + supg->n_e, g->edges + e);
		supg->edges[supg->n_e].source = u;
		supg->edges[supg->n_e].target = v;
		kh_long_int_set(node_map, GET_CODE(pu, e), u);
		kh_long_int_set(node_map, GET_CODE(e, pv), v);
		supg->n_v += 2;
		++supg->n_e;
	} else {
		supg->nodes = realloc(supg->nodes, (supg->n_v + 1)
				* sizeof(struct asm_node_t));
		memset(supg->nodes + supg->n_v, 0, sizeof(struct asm_node_t));
		int u = supg->n_v;
		kh_long_int_set(node_map, GET_CODE(pu, e), u);
		kh_long_int_set(node_map, GET_CODE(e, pv), u);
		++supg->n_v;
	}
}

void get_big_kmer(int e1, int e2, struct asm_graph_t *g, char **big_kmer)
{
	char *seq1, *seq2;
	decode_seq(&seq1, g->edges[e1].seq, g->edges[e1].seq_len);
	decode_seq(&seq2, g->edges[e2].seq, g->edges[e2].seq_len);
	*big_kmer = calloc(g->ksize + 3, sizeof(char));
	int len1 = strlen(seq1);
	strcpy(*big_kmer, seq1 + len1 - (g->ksize + 1));
	(*big_kmer)[g->ksize + 1] = seq2[g->ksize];
	free(seq1);
	free(seq2);
}

int get_big_kmer_count(char *big_kmer, int (*kmer_count)(char *))
{
	return 1000;
	return kmer_count(big_kmer) >= 2;
}

void add_super_edge(int mid, int e1, int e2, struct asm_graph_t *supg,
		char *big_kmer, int count, khash_t(long_int) *node_map)
{
	int u = kh_long_int_get(node_map, GET_CODE(e1, mid));
	int v = kh_long_int_get(node_map, GET_CODE(mid, e2));
	uint32_t *bin_seq;
	encode_seq(&bin_seq, big_kmer);
	struct asm_edge_t *e = calloc(1, sizeof(struct asm_edge_t));
	e->count = count;
	e->seq = bin_seq;
	e->seq_len = strlen(big_kmer);
	e->n_holes = 0;
	e->p_holes = NULL;
	e->l_holes = NULL;
	e->source = u;
	e->target = v;
	e->rc_id = 0;
	e->barcodes = NULL;
	supg->edges = realloc(supg->edges, (supg->n_e + 1)
			* sizeof(struct asm_edge_t));
	supg->edges[supg->n_e] = *e;
	asm_add_node_adj(supg, u, supg->n_e);
	++supg->n_e;
}

void create_super_edges(struct asm_graph_t *g, struct asm_graph_t *supg,
		khash_t(long_int) *node_map, int (*kmer_count)(char *))
{
	for (int u = 0; u < g->n_v; ++u){
		if (100 * (u + 1) / g->n_v > 100 * u / g->n_v)
			log_debug("Processed %d/%d vertices (%d\%)", u + 1,
					g->n_v, 100 * (u + 1) / g->n_v);
		for (int i = 0; i < g->nodes[u].deg; ++i){
			int e1 = g->nodes[u].adj[i];
			int v = g->edges[e1].target;
			for (int j = 0; j < g->nodes[v].deg; ++j){
				int e2 = g->nodes[v].adj[j];
				int w = g->edges[e2].target;
				char *big_kmer;
				get_big_kmer(e1, e2, g, &big_kmer);
				int count = get_big_kmer_count(big_kmer, kmer_count);
				if (count >= 2){
					add_super_edge(v, e1, e2, supg, big_kmer,
							count, node_map);
				}
				free(big_kmer);
			}
		}
	}
}

void assign_reverse_complement(struct asm_graph_t *g, struct asm_graph_t *supg,
		khash_t(long_int) *node_map)
{
	for (int u = 0; u < g->n_v; ++u){
		for (int i = 0; i < g->nodes[u].deg; ++i){
			int e = g->nodes[u].adj[i];
			int u_rc = g->nodes[u].rc_id;
			int e_rc = g->edges[e].rc_id;

			int supu = kh_long_int_get(node_map, GET_CODE(u, e));
			int supu_rc = kh_long_int_get(node_map,
					GET_CODE(e_rc, u_rc));
			supg->nodes[supu].rc_id = supu_rc;
		}
	}

	for (int u = 0; u < supg->n_v; ++u){
		for (int i = 0; i < supg->nodes[u].deg; ++i){
			int e = supg->nodes[u].adj[i];
			int v = supg->edges[e].target;
			int v_rc = supg->nodes[v].rc_id;
			int u_rc = supg->nodes[u].rc_id;

			int e_rc = -1;
			for (int j = 0; j < supg->nodes[v_rc].deg; ++j){
				e_rc = supg->nodes[v_rc].adj[j];
				int tmp = supg->edges[e_rc].target;
				if (tmp == u_rc)
					break;
			}
			if (e_rc == -1)
				log_error("Something went wrong, e_rc not found at node %d, edge %d",
						u, e);
			supg->edges[e].rc_id = e_rc;
		}
	}
}

void upsize_graph(struct asm_graph_t *g, struct asm_graph_t *supg,
		int (*kmer_count)(char *))
{
	khash_t(long_int) *node_map = kh_init(long_int);
	log_info("Creating super nodes");
	for (int e = 0; e < g->n_e; ++e){
		if (100 * (e + 1) / g->n_e > 100 * e / g->n_e)
			log_info("Processed %d/%d edges (%d\%)", e + 1,
					g->n_e, 100 * (e + 1) / g->n_e);
		create_super_nodes(g, e, supg, node_map);
	}
	log_info("Creating super edges");
	create_super_edges(g, supg, node_map, kmer_count);
	log_info("Assigning reverse complement id for nodes and edges");
	assign_reverse_complement(g, supg, node_map);
	kh_destroy(long_int, node_map);

	struct asm_graph_t g1;
	asm_condense(supg, &g1);
	asm_graph_destroy(supg);
	*supg = g1;
}

void resolve_multi_kmer(char *out_dir, struct asm_graph_t *g, int lastk, int (*kmer_count)(char *))
{
	for (int k = g->ksize; k <= lastk; ++k){
		log_info("Resolving using kmer of size %d", k);
		struct asm_graph_t *supg = calloc(1, sizeof(struct asm_graph_t));
		upsize_graph(g, supg, kmer_count);
	}
	save_graph_info(out_dir, g, "kmer_resolve");
}

