#include "minimizers/count_barcodes.h"
#include "complex_resolve.h"
#include "khash.h"
#include "log.h"
#include "resolve.h"
#include "verbose.h"
#include "atomic.h"
#include "helper.h"
#include "dqueue.h"
#include "process.h"
#include "utils.h"
#include "kmhash.h"
#include "kmer.h"
#include "fastq_producer.h"


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
		khash_t(long_int) *node_map_fw, khash_t(long_int) *node_map_bw,
		struct mini_hash_t *kmer_table)
{
	int e_rc = g->edges[e].rc_id;
	int pu = g->edges[e].source;
	if (pu == -1)
		return;
	int pv = g->edges[e].target;
	int ok = 1;
	//if (g->edges[e].seq_len > g->ksize + 1 && g->edges[e].seq_len <= g->ksize + 10){
	//	int unmatch = 0;
	//	char *seq;
	//	char *kmer = calloc(g->ksize + 3, sizeof(char));
	//	char *kmer_rc = calloc(g->ksize + 3, sizeof(char));
	//	decode_seq(&seq, g->edges[e].seq, g->edges[e].seq_len);
	//	for (int i = 0; i + g->ksize + 2 <= g->edges[e].seq_len; ++i){
	//		strncpy(kmer, seq + i, g->ksize + 2);
	//		strcpy(kmer_rc, kmer);
	//		flip_reverse(kmer_rc);
	//		if (!get_big_kmer_count(kmer, kmer_table) &&
	//			!get_big_kmer_count(kmer_rc, kmer_table)){
	//			++unmatch;
	//			int total_kmer = g->edges[e].seq_len - (g->ksize + 2) + 1;
	//			//if (1.0 * unmatch / total_kmer > 0.05){
	//			if (unmatch){
	//				ok = 0;
	//				log_info("Delete edge %d len %d match %d",
	//						e, g->edges[e].seq_len,
	//						total_kmer - unmatch);
	//				break;
	//			}
	//		}
	//	}
	//	free(seq);
	//	free(kmer);
	//	free(kmer_rc);
	//}
	if ((int) g->edges[e].seq_len > g->ksize + 1){
		int u = supg->n_v;
		int v = supg->n_v + 1;
		kh_long_int_set(node_map_fw, GET_CODE(pu, e), u);
		kh_long_int_set(node_map_bw, GET_CODE(e, pv), v);
		supg->n_v += 2;

		if (ok){
			asm_add_node_adj(supg, u, supg->n_e);
			asm_clone_seq(supg->edges + supg->n_e, g->edges + e);
			supg->edges[supg->n_e].source = u;
			supg->edges[supg->n_e].target = v;
			//log_debug("(pu, e) = (%d, %d) -> (u) = (%d)", pu, e, u);
			//log_debug("(e, pv) = (%d, %d) -> (v) = (%d)", e, pv, v);
			++supg->n_e;
		}
	} else {
		int u = supg->n_v;
		kh_long_int_set(node_map_fw, GET_CODE(pu, e), u);
		kh_long_int_set(node_map_bw, GET_CODE(e, pv), u);
		++supg->n_v;
	}
}

void get_big_kmer(int e1, int e2, struct asm_graph_t *g, char **big_kmer)
{
	//char *seq1, *seq2;
	//decode_seq(&seq1, g->edges[e1].seq, g->edges[e1].seq_len);
	//decode_seq(&seq2, g->edges[e2].seq, g->edges[e2].seq_len);
	//*big_kmer = calloc(g->ksize + 3, sizeof(char));
	*big_kmer = malloc((g->ksize + 3) * sizeof(char));

	for (int i = 0; i < g->ksize + 1; ++i)
		(*big_kmer)[i] = int_to_base(__binseq_get(g->edges[e1].seq,
				g->edges[e1].seq_len - (g->ksize + 1) + i));
	(*big_kmer)[g->ksize + 1] = int_to_base(__binseq_get(g->edges[e2].seq,
						g->ksize));
	//int len1 = strlen(seq1);
	//strcpy(*big_kmer, seq1 + len1 - (g->ksize + 1));
	//(*big_kmer)[g->ksize + 1] = seq2[g->ksize];
	//seq2[g->ksize] = 0;
	//free(seq1);
	//free(seq2);
	(*big_kmer)[g->ksize + 2] = '\0';
}

int get_big_kmer_count(char *big_kmer, struct mini_hash_t *kmer_table)
{
	int ksize = strlen(big_kmer);
	uint8_t *seq = calloc((ksize + 3) >> 2, 1);
	for (int i = 0; i < ksize; i++) {
		km_shift_append(seq, ksize, (ksize+3)>>2, nt4_table[big_kmer[i]]);
	}
	uint64_t key_new = MurmurHash3_x64_64(seq, (ksize +3 ) >> 2);
	uint64_t *slot = mini_get(kmer_table, key_new);
	int res = 0;
	if (slot != (uint64_t *)EMPTY_SLOT)
		res = *slot;
	free(seq);
	return res;
}

void add_super_edge(int mid, int e1, int e2, struct asm_graph_t *supg,
		char *big_kmer, int count, khash_t(long_int) *node_map_fw,
		khash_t(long_int) *node_map_bw)
{
	int u = kh_long_int_get(node_map_bw, GET_CODE(e1, mid));
	int v = kh_long_int_get(node_map_fw, GET_CODE(mid, e2));
	uint32_t *bin_seq;
	encode_seq(&bin_seq, big_kmer);
	struct asm_edge_t *e = supg->edges + supg->n_e;
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
	asm_add_node_adj(supg, u, supg->n_e);
	++supg->n_e;
}

void add_super_edge_multi(int mid, int e1, int e2, struct asm_graph_t *supg,
		char *big_kmer, int count, khash_t(long_int) *node_map_fw,
		khash_t(long_int) *node_map_bw, struct super_edges_bundle_t *bundle)
{
	int u = kh_long_int_get(node_map_bw, GET_CODE(e1, mid));
	int v = kh_long_int_get(node_map_fw, GET_CODE(mid, e2));
	uint32_t *bin_seq;
	encode_seq(&bin_seq, big_kmer);
	pthread_mutex_lock(bundle->edge_id_lock);
	int edge_id = supg->n_e;
	++(supg->n_e);
	pthread_mutex_unlock(bundle->edge_id_lock);

	struct asm_edge_t *e = supg->edges + edge_id;
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
	pthread_mutex_lock(bundle->supg_node_locks + u);
	asm_add_node_adj(supg, u, edge_id);
	pthread_mutex_unlock(bundle->supg_node_locks + u);
	//pthread_mutex_lock(bundle->supg_node_locks);
	//asm_add_node_adj(supg, u, edge_id);
	//pthread_mutex_unlock(bundle->supg_node_locks);
}

void create_super_edges(struct asm_graph_t *g, struct asm_graph_t *supg,
		khash_t(long_int) *node_map_fw, khash_t(long_int) *node_map_bw,
		struct mini_hash_t *kmer_table)
{
	int *accept = calloc(g->n_v, sizeof(int));
	int count_resolve = 0;
	for (int e1 = 0; e1 < g->n_e; ++e1){
		//log_debug("-----");
		if (g->edges[e1].source == -1)
			continue;
		int u = g->edges[e1].target;
		int u_rc = g->nodes[u].rc_id;
		for (int i = 0; i < g->nodes[u].deg; ++i){
			int e2 = g->nodes[u].adj[i];

			char *big_kmer;
			get_big_kmer(e1, e2, g, &big_kmer);

			char *big_kmer_rc = calloc(g->ksize + 3, sizeof(char));
			strcpy(big_kmer_rc, big_kmer);
			flip_reverse(big_kmer_rc);
			int count1 = get_big_kmer_count(big_kmer, kmer_table) ;
			int count2 = get_big_kmer_count(big_kmer_rc, kmer_table);
			int count = count1 + count2;
			//if (g->nodes[u].deg > 1) {
			//	if (count1 + count2 > 5000) {
			//		log_debug("this kmer exist>5000: %s", big_kmer);
			//	}
			//	log_debug("    To %d: %d %d", e2, count1, count2);
			//}
			if ((g->nodes[u].deg == 1 && g->nodes[u_rc].deg == 1)
				|| count >= 3){
				add_super_edge(u, e1, e2, supg, big_kmer,
						count * (g->ksize + 2 - g->ksize_count),
						node_map_fw, node_map_bw);
				++accept[u];
			} else {
				++count_resolve;
			}

			free(big_kmer);
			free(big_kmer_rc);
		}
	}
	for (int u = 0; u < g->n_v; ++u){
		if (accept[u] > 0)
			continue;
		int u_rc = g->nodes[u].rc_id;
		for (int i = 0; i < g->nodes[u_rc].deg; ++i){
			int e1 = g->edges[g->nodes[u_rc].adj[i]].rc_id;
			for (int j = 0; j < g->nodes[u].deg; ++j){
				int e2 = g->nodes[u].adj[j];
				char *big_kmer;
				get_big_kmer(e1, e2, g, &big_kmer);
				add_super_edge(u, e1, e2, supg, big_kmer, 0,
						node_map_fw, node_map_bw);
				free(big_kmer);
			}
		}
	}
	free(accept);
	log_info("Total resolve using big kmer: %d edges", count_resolve);
}

void assign_reverse_complement(struct asm_graph_t *g, struct asm_graph_t *supg,
		khash_t(long_int) *node_map_fw, khash_t(long_int) *node_map_bw)
{
	for (int u = 0; u < g->n_v; ++u){
		for (int i = 0; i < g->nodes[u].deg; ++i){
			int e = g->nodes[u].adj[i];
			int u_rc = g->nodes[u].rc_id;
			int e_rc = g->edges[e].rc_id;

			int supu = kh_long_int_get(node_map_fw, GET_CODE(u, e));
			int supu_rc = kh_long_int_get(node_map_bw,
					GET_CODE(e_rc, u_rc));
			supg->nodes[supu].rc_id = supu_rc;
		}
	}

	for (int v = 0; v < g->n_v; ++v){
		int v_rc = g->nodes[v].rc_id;
		for (int i = 0; i < g->nodes[v_rc].deg; ++i){
			int e_rc = g->nodes[v_rc].adj[i];
			int e = g->edges[e_rc].rc_id;

			int supv = kh_long_int_get(node_map_bw, GET_CODE(e, v));
			int supv_rc = kh_long_int_get(node_map_fw,
					GET_CODE(v_rc, e_rc));
			supg->nodes[supv].rc_id = supv_rc;
		}
	}

	for (int u = 0; u < supg->n_v; ++u){
		for (int i = 0; i < supg->nodes[u].deg; ++i){
			int e1 = supg->nodes[u].adj[i];
			int v = supg->edges[e1].target;
			int v_rc = supg->nodes[v].rc_id;
			int u_rc = supg->nodes[u].rc_id;

			int e_rc = -1;
			for (int j = 0; j < supg->nodes[v_rc].deg; ++j){
				int e2 = supg->nodes[v_rc].adj[j];
				if (u_rc != supg->edges[e2].target)
					continue;
				if (is_reverse_complement(supg, e1, e2)){
					e_rc = e2;
					break;
				}
			}
			if (e_rc == -1)
				log_error("Something went wrong, e_rc not found at node %d, edge %d, u_rc %d v_rc %d",
						u, e1, u_rc, v_rc);
			supg->edges[e1].rc_id = e_rc;
		}
	}
}

void estimate_something(struct asm_graph_t *g, int *count_edge, int *count_node)
{
	int c_e = 0, c_n = 0;
	for (int i_e = 0; i_e < g->n_e; i_e++) {
		if (g->edges[i_e].source == -1 )
			continue;
		if ((int) g->edges[i_e].seq_len > g->ksize + 1){
			c_n += 2;
			c_e+=1;

		} else {
			c_n++;
		}
	}

	for (int e = 0; e < g->n_e; ++e){
		if (g->edges[e].source == -1)
			continue;
		c_e += g->nodes[g->edges[e].target].deg;
	}

	*count_edge = c_e;
	*count_node = c_n;
}

void upsize_graph(struct opt_proc_t *opt, int super_k, struct asm_graph_t *g,
		struct asm_graph_t *supg)
{
	khash_t(long_int) *node_map_fw = kh_init(long_int);
	khash_t(long_int) *node_map_bw = kh_init(long_int);
	struct mini_hash_t *kmer_table = get_kmer_count_from_kmc(super_k,
					opt->n_files, opt->files_1, opt->files_2,
					opt->n_threads, opt->mmem, opt->out_dir);
	//struct mini_hash_t *kmer_table = count_kmer_simple(opt, super_k,
	//					opt->files_1[0], opt->files_2[0]);
	log_info("Creating super nodes");
	int estimate_node = 0, estimate_edge = 0;
	estimate_something(g, &estimate_edge, &estimate_node);
	supg->nodes = calloc(estimate_node, sizeof(struct asm_node_t));
	supg->edges = calloc(estimate_edge, sizeof(struct asm_edge_t));
	for (int e = 0; e < g->n_e; ++e){
		//if (100 * (e + 1) / g->n_e > 100 * e / g->n_e)
		//	log_info("Processed %d/%d edges (%d\%)", e + 1,
		//			g->n_e, 100 * (e + 1) / g->n_e);
		create_super_nodes(g, e, supg, node_map_fw, node_map_bw,
				kmer_table);
	}

	log_info("Creating super edges");

	create_super_edges_multi(opt, g, supg, node_map_fw, node_map_bw,
			kmer_table);
	destroy_mini_hash(kmer_table);

	supg->nodes = realloc(supg->nodes, supg->n_v * sizeof(struct asm_node_t));
	supg->edges = realloc(supg->edges, supg->n_e * sizeof(struct asm_edge_t));
	log_info("node %d estimated node %d, edge %d estimated edge %d",
			supg->n_v, estimate_node, supg->n_e, estimate_edge);

	log_info("Assigning reverse complement id for nodes and edges");

	assign_reverse_complement_multi(opt, g, supg, node_map_fw, node_map_bw);

	kh_destroy(long_int, node_map_fw);
	kh_destroy(long_int, node_map_bw);
	supg->ksize = super_k;
	if (super_k % 2){
		log_info("Condesing graph");
		struct asm_graph_t g1;
		asm_condense(supg, &g1);
		asm_graph_destroy(supg);
		*supg = g1;

		//log_info("Resolving graph");
		//struct asm_graph_t g2;
		//resolve_graph_operation(supg, &g2);
	}
	test_asm_graph(supg);
	supg->ksize_count = g->ksize_count;
}

void resolve_multi_kmer(struct opt_proc_t *opt, struct asm_graph_t *g, int lastk)
{
	for (int k = g->ksize + 1; k <= lastk; ++k){
		log_info("Resolving using kmer of size %d", k);
		struct asm_graph_t supg;
		memset(&supg, 0, sizeof(struct asm_graph_t));
		upsize_graph(opt, k, g, &supg);

		log_info("Resolving done kmer of size %d", k);
		asm_graph_destroy(g);
		*g = supg;
		save_graph_info(opt->out_dir, g, "kmer_resolve");
	}
}

int is_reverse_complement(struct asm_graph_t *g, int e1, int e2)
{
	int u = g->edges[e1].source;
	int v = g->edges[e1].target;
	int u_rc = g->edges[e2].target;
	int v_rc = g->edges[e2].source;
	if (g->nodes[u].rc_id != u_rc)
		return 0;
	if (g->nodes[v].rc_id != v_rc)
		return 0;
	int base = __binseq_get(g->edges[e1].seq, g->ksize);
	int base_rc = __binseq_get(g->edges[e2].seq,
					g->edges[e2].seq_len - g->ksize - 1);
	return (base ^ base_rc) == 3;
}

void *assign_kmer_pos_eid(void *data)
{
	struct kmer_eid_bundle_t *bundle = (struct kmer_eid_bundle_t *) data;

	struct asm_graph_t *g = bundle->g;
	char *kmer = calloc(g->ksize + 1, sizeof(char));
	do{
		pthread_mutex_lock(bundle->e_id_lock);
		int e = *(bundle->e_id);
		++(*bundle->e_id);
		pthread_mutex_unlock(bundle->e_id_lock);
		if (e >= g->n_e)
			break;
		char *seq;
		decode_seq(&seq, g->edges[e].seq, g->edges[e].seq_len);
		for (int i = 0; i + g->ksize <= (int) g->edges[e].seq_len; ++i){
			strncpy(kmer, seq + i, g->ksize);
			uint64_t hash = MurmurHash3_x64_64((uint8_t *) kmer,
								g->ksize);
			kh_long_int_set(bundle->kmer_pos, hash, e);
			//log_debug("Kmer %s hash code is %ld", kmer, hash);
		}
		free(seq);
	} while(1);
	free(kmer);
	return NULL;
}

void get_kmer_edge_id(struct opt_proc_t *opt, struct asm_graph_t *g, khash_t(long_int) *kmer_eid)
{
	struct kmer_eid_bundle_t *worker_bundle = calloc(opt->n_threads,
			sizeof(struct kmer_eid_bundle_t));
	khash_t(long_int) **each_kmer_eid = calloc(opt->n_threads,
			sizeof(khash_t(long_int) *));
	for (int i = 0; i < opt->n_threads; ++i)
		each_kmer_eid[i] = kh_init(long_int);

	int e_id = 0;

	pthread_mutex_t e_id_lock;
	pthread_mutex_init(&e_id_lock, NULL);

	for (int i = 0; i < opt->n_threads; ++i){
		worker_bundle[i].g = g;
		worker_bundle[i].kmer_pos = each_kmer_eid[i];
		worker_bundle[i].e_id = &e_id;
		worker_bundle[i].e_id_lock = &e_id_lock;
	}

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_t *threads = calloc(opt->n_threads, sizeof(pthread_t));
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(threads + i, &attr, assign_kmer_pos_eid,
				worker_bundle + i);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(threads[i], NULL);

	for (int i = 0; i < opt->n_threads; ++i){
		for (khiter_t it = kh_begin(each_kmer_eid[i]);
			it != kh_end(each_kmer_eid[i]); ++it){
			if (!kh_exist(each_kmer_eid[i], it))
				continue;
			uint64_t key = kh_key(each_kmer_eid[i], it);
			int val = kh_val(each_kmer_eid[i], it);
			kh_long_int_set(kmer_eid, key, val);
		}
		kh_destroy(long_int, each_kmer_eid[i]);
	}

	//for (khiter_t it = kh_begin(kmer_eid); it != kh_end(kmer_eid); ++it){
	//	if (!kh_exist(kmer_eid, it))
	//		continue;
	//	log_debug("Kmer %ld appears in edge %d", kh_key(kmer_eid, it),
	//			kh_val(kmer_eid, it));
	//}
	free(each_kmer_eid);
	free(worker_bundle);

}

void create_super_edges_multi(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *supg, khash_t(long_int) *node_map_fw,
		khash_t(long_int) *node_map_bw, struct mini_hash_t *kmer_table)
{
	struct super_edges_bundle_t *worker_bundle = calloc(opt->n_threads,
			sizeof(struct super_edges_bundle_t));

	int node_id = 0;
	pthread_mutex_t node_id_lock;
	pthread_mutex_init(&node_id_lock, NULL);

	pthread_mutex_t edge_id_lock;
	pthread_mutex_init(&edge_id_lock, NULL);

	pthread_mutex_t *supg_node_locks = calloc(supg->n_v,
						sizeof(pthread_mutex_t));
	for (int i = 0; i < supg->n_v; ++i)
		pthread_mutex_init(supg_node_locks + i, NULL);

	int *count_resolve = calloc(30, sizeof(int));
	for (int i = 0; i < opt->n_threads; ++i){
		worker_bundle[i].g = g;
		worker_bundle[i].supg = supg;
		worker_bundle[i].node_id = &node_id;
		worker_bundle[i].node_id_lock = &node_id_lock;
		worker_bundle[i].node_map_fw = node_map_fw;
		worker_bundle[i].node_map_bw = node_map_bw;
		worker_bundle[i].kmer_table = kmer_table;
		worker_bundle[i].edge_id_lock = &edge_id_lock;
		worker_bundle[i].supg_node_locks = supg_node_locks;
		worker_bundle[i].count_resolve = &count_resolve;
	}

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_t *threads = calloc(opt->n_threads, sizeof(pthread_t));
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(threads + i, &attr, create_super_edges_ite,
				worker_bundle + i);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(threads[i], NULL);
	log_info("count resolve at kmer%d is: %d", g->ksize, count_resolve[0]);
	log_info("count node disconnect %d is: %d", g->ksize, count_resolve[1]);
	log_info("count node keep many many %d is: %d", g->ksize, count_resolve[2]);
	for (int i = 3; i < 30; i++) {
		log_info("count deg in out %d: %d", i-3, count_resolve[i]);
	}

	free(threads);
	free(supg_node_locks);
	free(worker_bundle);
}

void *create_super_edges_ite(void *data)
{
	struct super_edges_bundle_t *bundle = (struct super_edges_bundle_t *) data;
	struct asm_graph_t *g = bundle->g;
	struct asm_graph_t *supg = bundle->supg;
	khash_t(long_int) *node_map_fw = bundle->node_map_fw;
	khash_t(long_int) *node_map_bw = bundle->node_map_bw;
	struct mini_hash_t *kmer_table = bundle->kmer_table;
	do{
		pthread_mutex_lock(bundle->node_id_lock);
		int u = *bundle->node_id;
		++(*bundle->node_id);
		pthread_mutex_unlock(bundle->node_id_lock);
		if (u >= g->n_v)
			break;

		int accept = 0;
		int u_rc = g->nodes[u].rc_id;
		for (int i = 0; i < g->nodes[u_rc].deg; ++i){
			int e1 = g->edges[g->nodes[u_rc].adj[i]].rc_id;
			for (int j = 0; j < g->nodes[u].deg; ++j){
				int e2 = g->nodes[u].adj[j];
				char *big_kmer;
				get_big_kmer(e1, e2, g, &big_kmer);

				char *big_kmer_rc = calloc(g->ksize + 3, sizeof(char));
				strcpy(big_kmer_rc, big_kmer);
				flip_reverse(big_kmer_rc);
				int count1 = get_big_kmer_count(big_kmer, kmer_table) ;
				int count2 = get_big_kmer_count(big_kmer_rc, kmer_table);
				int count = count1 + count2;
				if ((g->nodes[u].deg == 1 && g->nodes[u_rc].deg == 1)
					|| count >= 3){
					add_super_edge_multi(u, e1, e2, supg, big_kmer,
							count * (g->ksize + 2 - g->ksize_count),
							node_map_fw, node_map_bw,
							bundle);
					accept++;
				} else {
					atomic_add_and_fetch32(*bundle->count_resolve, 1);
				}
				free(big_kmer);
				free(big_kmer_rc);
			}
		}
		int t = g->nodes[u_rc].deg * 5 + g->nodes[u].deg;
		atomic_add_and_fetch32(*bundle->count_resolve + t + 3, 1);
		if (accept != 0) {
			atomic_add_and_fetch32(*bundle->count_resolve + 1, 1);
		} else if (accept > 1) {
			atomic_add_and_fetch32(*bundle->count_resolve + 2, 1);
		}

		//if (accept)
		//	continue;


		//for (int i = 0; i < g->nodes[u_rc].deg; ++i){
		//	int e1 = g->edges[g->nodes[u_rc].adj[i]].rc_id;
		//	for (int j = 0; j < g->nodes[u].deg; ++j){
		//		int e2 = g->nodes[u].adj[j];
		//		char *big_kmer;
		//		get_big_kmer(e1, e2, g, &big_kmer);

		//		add_super_edge_multi(u, e1, e2, supg, big_kmer, 0,
		//				node_map_fw, node_map_bw,
		//				bundle);
		//		free(big_kmer);
		//	}
		//}
	}while (1);
	return NULL;
}

void assign_reverse_complement_multi(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct asm_graph_t *supg, khash_t(long_int) *node_map_fw,
		khash_t(long_int) *node_map_bw)
{
	struct assign_rc_bundle_t *worker_bundle = calloc(opt->n_threads,
			sizeof(struct assign_rc_bundle_t));

	int node_id = 0;
	pthread_mutex_t node_id_lock;
	pthread_mutex_init(&node_id_lock, NULL);

	for (int i = 0; i < opt->n_threads; ++i){
		worker_bundle[i].g = g;
		worker_bundle[i].supg = supg;
		worker_bundle[i].node_map_fw = node_map_fw;
		worker_bundle[i].node_map_bw = node_map_bw;
		worker_bundle[i].node_id = &node_id;
		worker_bundle[i].node_id_lock = &node_id_lock;
	}

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_t *threads = calloc(opt->n_threads, sizeof(pthread_t));
	log_info("Assigning reverse complement for nodes");

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(threads + i, &attr, assign_node_rc_ite,
				worker_bundle + i);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(threads[i], NULL);

	log_info("Assigning reverse complement for edges");
	node_id = 0;
	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(threads + i, &attr, assign_edge_rc_ite,
				worker_bundle + i);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(threads[i], NULL);

	free(threads);
	free(worker_bundle);
}

void *assign_node_rc_ite(void *data)
{
	struct assign_rc_bundle_t *bundle = (struct assign_rc_bundle_t *) data;
	struct asm_graph_t *g = bundle->g;
	struct asm_graph_t *supg = bundle->supg;
	khash_t(long_int) *node_map_fw = bundle->node_map_fw;
	khash_t(long_int) *node_map_bw = bundle->node_map_bw;
	
	do{
		pthread_mutex_lock(bundle->node_id_lock);
		int u = *bundle->node_id;
		++(*bundle->node_id);
		pthread_mutex_unlock(bundle->node_id_lock);

		if (u >= g->n_v)
			break;

		for (int i = 0; i < g->nodes[u].deg; ++i){
			int e = g->nodes[u].adj[i];
			int u_rc = g->nodes[u].rc_id;
			int e_rc = g->edges[e].rc_id;

			int supu = kh_long_int_get(node_map_fw, GET_CODE(u, e));
			int supu_rc = kh_long_int_get(node_map_bw,
					GET_CODE(e_rc, u_rc));
			supg->nodes[supu].rc_id = supu_rc;
		}

		int u_rc = g->nodes[u].rc_id;
		for (int i = 0; i < g->nodes[u_rc].deg; ++i){
			int e_rc = g->nodes[u_rc].adj[i];
			int e = g->edges[e_rc].rc_id;

			int supu = kh_long_int_get(node_map_bw, GET_CODE(e, u));
			int supu_rc = kh_long_int_get(node_map_fw,
					GET_CODE(u_rc, e_rc));
			supg->nodes[supu].rc_id = supu_rc;
		}
	} while(1);
	return NULL;
}

void *assign_edge_rc_ite(void *data)
{
	struct assign_rc_bundle_t *bundle = (struct assign_rc_bundle_t *) data;
	struct asm_graph_t *supg = bundle->supg;
	khash_t(long_int) *node_map_fw = bundle->node_map_fw;
	khash_t(long_int) *node_map_bw = bundle->node_map_bw;
	
	do{
		pthread_mutex_lock(bundle->node_id_lock);
		int u = *bundle->node_id;
		++(*bundle->node_id);
		pthread_mutex_unlock(bundle->node_id_lock);

		if (u >= supg->n_v)
			break;

		for (int i = 0; i < supg->nodes[u].deg; ++i){
			int e1 = supg->nodes[u].adj[i];
			int v = supg->edges[e1].target;
			int v_rc = supg->nodes[v].rc_id;
			int u_rc = supg->nodes[u].rc_id;

			int e_rc = -1;
			for (int j = 0; j < supg->nodes[v_rc].deg; ++j){
				int e2 = supg->nodes[v_rc].adj[j];
				if (u_rc != supg->edges[e2].target)
					continue;
				if (is_reverse_complement(supg, e1, e2)){
					e_rc = e2;
					break;
				}
			}
			if (e_rc == -1)
				log_error("Something went wrong, e_rc not found at node %d, edge %d, u_rc %d v_rc %d",
						u, e1, u_rc, v_rc);
			supg->edges[e1].rc_id = e_rc;
		}
	} while(1);
	return NULL;
}

struct mini_hash_t *count_kmer_simple(struct opt_proc_t *opt, int ksize,
			char *R1_path, char *R2_path)
{
	struct producer_bundle_t *producer_bundles;
	producer_bundles = init_fastq_pair(opt->n_threads, 1, &R1_path,
						&R2_path);

	struct count_kmer_bundle_t *worker_bundles = calloc(opt->n_threads,
			sizeof(struct count_kmer_bundle_t));

	struct mini_hash_t *h;
	init_mini_hash(&h, 20);

	pthread_mutex_t lock;
	pthread_mutex_init(&lock, NULL);
	for (int i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = producer_bundles->q;
		worker_bundles[i].ksize = ksize;
		worker_bundles[i].lock = &lock;
		worker_bundles[i].h = &h;
	}

	pthread_t *producer_threads, *worker_threads;
	/* FIXME: not actually opt->n_files, noob */

	opt->n_files = 1;
	producer_threads = calloc(opt->n_files, sizeof(pthread_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	for (int i = 0; i < opt->n_files; ++i)
		pthread_create(producer_threads + i, &attr, fastq_producer,
		               producer_bundles + i);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_create(worker_threads + i, &attr, count_kmer_ite,
		               worker_bundles + i);

	for (int i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (int i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	free_fastq_pair(producer_bundles, 1);
	free(worker_bundles);

	free(producer_threads);
	free(worker_threads);
	return h;
}

void *count_kmer_ite(void *data)
{
	struct count_kmer_bundle_t *bundle = (struct count_kmer_bundle_t *) data;
	struct dqueue_t *q = bundle->q;
	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2, input_format, mapper_algo;

	int64_t n_reads;
	int64_t *gcnt_reads;
	uint64_t barcode;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf) {
			break;
		}
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		buf1 = ext_buf->R1_buf;
		buf2 = ext_buf->R2_buf;
		input_format = ext_buf->input_format;

		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
			      get_read_from_fq(&read1, buf1, &pos1) :
			      get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
			      get_read_from_fq(&read2, buf2, &pos2) :
			      get_read_from_fa(&read2, buf2, &pos2);


			if (rc1 == READ_FAIL || rc2 == READ_FAIL) {
				log_warn("buf1 %s", buf1);
				log_warn("buf2 %s", buf2);
				log_error("\nWrong format file when build barcode scaffold\n");
			}

			add_count(&read1, bundle->ksize, bundle->h);
			add_count(&read2, bundle->ksize, bundle->h);
			if (rc1 == READ_END)
				break;
		}
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

void add_count(struct read_t *r, int ksize, struct mini_hash_t **h)
{
	for (int i = 0; i + ksize <= r->len; ++i){
		uint64_t hash = MurmurHash3_x64_64(r->seq + i, ksize);
		mini_put(h, hash);
	}
}

