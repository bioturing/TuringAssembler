#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define RESET "\x1B[0m"

#include "io_utils.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"
#include "khash.h"
#include "queue.h"

#include "assembly_graph.h"

#define MIN_COMPONENT 250
#define MIN_BRIDGE_LEG 5000
#define MIN_SHARED_BARCODE 50
#define MIN_LAYER 100
#define MIN_VISITED_NODES 1
#define MIN_COMPLEX_REGION 10000
#define MIN_ITER_RESOLVE 100

#define is_lg_leg(e, v, i, cnt) (((e[v[i].adj[0]].seq_len > MIN_BRIDGE_LEG) + \
		(e[v[i].adj[1]].seq_len > MIN_BRIDGE_LEG)) == cnt)
#define ISOLATE_NODE(v, u) \
	v[u].deg = 0; \
	v[v[u].rc_id].deg = 0;
#define SAME_EDGE(e, i, j) (e[i].target == e[j].target && e[i].source == e[j].source)

#define REFINE_EDGE(e, v, i) \
	v[e[i].target].deg = 1; \
	v[v[e[i].target].rc_id].deg = 1; \
	v[e[i].source].deg = 1; \
	v[v[e[i].source].rc_id].deg = 1; \

KHASH_SET_INIT_INT(khInt);

static int is_bridge(struct asm_graph_t *g, gint_t i);
static void asm_edge_cc(struct asm_graph_t *g, gint_t *id_edge, gint_t **ret_size);


static void condense_3_seq(gint_t e1, gint_t e2, gint_t e3, struct asm_edge_t *e, 
			struct asm_node_t *v, gint_t ksize)
{
	v[e[e1].target].deg = 0;
	v[e[e3].source].deg = 0;
	uint32_t *seq;
	gint_t len;

	concat_3_seq(e[e1], e[e2], e[e3], ksize, &len, &seq);
	e[e1].count = e[e1].count + e[e2].count + e[e3].count;
	e[e1].seq = seq;
	e[e1].seq_len = len;
	clone_edge_rev(e + e[e1].rc_id, e + e1);
	__VERBOSE("Condensed length of %d: %d\n", e1, len);

	e[e1].target = e[e3].target;

	v[e[e1].target].deg = 0;
	v[e[e3].source].deg = 0;
	v[v[e[e1].target].rc_id].deg = 0;
	v[v[e[e3].source].rc_id].deg = 0;
	
	e[e3].seq_len = 0;
	e[e[e3].rc_id].seq_len = 0;
}

static int is_bridge(struct asm_graph_t *g, gint_t i)
{
	struct asm_edge_t *e = g->edges;
	struct asm_node_t *v = g->nodes; 
	gint_t src = e[i].source;
	gint_t dest = e[i].target;
	gint_t r_src = v[src].rc_id;
	int k, j;
	if (!(v[r_src].deg == 2 && v[dest].deg == 2))
		return 0;
	if (is_lg_leg(e, v, r_src, 2) && is_lg_leg(e, v, dest, 2)){
		for (k = 0; k < v[r_src].deg; ++k)
			for (j = 0; j < v[dest].deg; ++j){
				if (v[r_src].adj[k] == v[dest].adj[j] || 
					v[r_src].adj[k] == e[v[dest].adj[j]].rc_id)
					return 0;
			}
		return 1;
	}
	return 0;
}

static void tandem_helper(struct asm_edge_t *e, struct asm_node_t *v,
			  gint_t u, khash_t(khInt) *set_v, uint32_t *comp_sz,
			  aqueue_t *q, khash_t(khInt) *lg, int is_first)
{
	gint_t next_e;
	int i, missing;
	gint_t rc_src, dest;
	rc_src = v[e[u].source].rc_id;
	dest = e[u].target;

	for (i = 0 ; i < v[dest].deg; ++i) {
		next_e = v[dest].adj[i]; //iterate outgoing edge of the target of u
		if (kh_get(khInt, set_v, next_e) == kh_end(set_v)){ // visited or not
			if (e[next_e].seq_len > MIN_BRIDGE_LEG){
				kh_put(khInt, lg, next_e, &missing);
				 continue;
			}
			aqueue_add(q, next_e); //add neighbor edge
			kh_put(khInt, set_v, next_e, &missing);
			*comp_sz += e[next_e].seq_len; //add size of the edge to total size
		}
	}
	if (!is_first)
		for (i = 0 ; i < v[rc_src].deg; ++i) {
			next_e = e[v[rc_src].adj[i]].rc_id; //iterate outgoing edge of the target of u
			if (kh_get(khInt, set_v, next_e) == kh_end(set_v)){ // visited or not
				if (e[next_e].seq_len > MIN_BRIDGE_LEG){
					kh_put(khInt, lg, next_e, &missing);
					 continue;
				}
				aqueue_add(q, next_e); //add neighbor edge
				kh_put(khInt, set_v, next_e, &missing);
				*comp_sz += e[next_e].seq_len; //add size of the edge to total size
		}
	}
}

static void print_hash_lg(khash_t(khInt) *h)
{
	khiter_t k;
	uint32_t key;
	for (k = kh_begin(h); k != kh_end(h); ++k){
		if (kh_exist(h, k)){
			key = kh_key(h, k);
			__VERBOSE("Key = %d\n", key);
		}
	}
}

static int is_source_edge(struct asm_graph_t *g, gint_t e_i, khash_t(khInt) *h)
{
	int i;
	struct asm_edge_t *e = g->edges;
	struct asm_node_t *v = g->nodes;
	gint_t dest = e[e_i].target;
	for (i = 0; i < v[dest].deg; ++i){
		if (kh_get(khInt, h, v[dest].adj[i]) == kh_end(h))
			return 0;
	}
	return 1;
}

static int is_sink_edge(struct asm_graph_t *g, gint_t e_i, khash_t(khInt) *h)
{
	int i;
	struct asm_edge_t *e = g->edges;
	struct asm_node_t *v = g->nodes;
	gint_t rc_src = v[e[e_i].source].rc_id;
	for (i = 0; i < v[rc_src].deg; ++i){
		if (kh_get(khInt, h, e[v[rc_src].adj[i]].rc_id) == kh_end(h))
			return 0;
	}
	return 1;
}

static int is_large_loop(struct asm_graph_t *g, khash_t(khInt) *h, khash_t(khInt) *comp_set)
{
	khiter_t k;
	uint32_t key;
	for (k = kh_begin(h); k != kh_end(h); ++k){
		if (kh_exist(h, k)){
			key = kh_key(h, k);
			if (is_source_edge(g, key, comp_set) && is_sink_edge(g, key, comp_set)){
				__VERBOSE(KRED "Self loop detect e.i is both source and sink %d\n" RESET, key);
				return 1;
			}
		}
	}
	return 0;
}

static void isolate_source_region(struct asm_edge_t *e, struct asm_node_t *v, gint_t *ei,
		int nei)
{
	int i, j;
	gint_t dest, rc_e;
	for (i = 0; i < nei; ++i){
		dest = e[ei[i]].target;
		for (j = 0; j < v[dest].deg; ++j){
			rc_e = e[v[dest].adj[j]].rc_id;
			v[e[rc_e].target].deg = 0;
		}
		v[dest].deg = 0;
	}
}

static void isolate_sink_region(struct asm_edge_t *e, struct asm_node_t *v, gint_t *eo,
		int neo)
{
	int i, j;
	gint_t rc_src, rc_e;
	for (i = 0; i < neo; ++i){
		rc_src = v[e[eo[i]].target].rc_id;
		for (j = 0; j < v[rc_src].deg; ++j){
			rc_e = e[v[rc_src].adj[j]].rc_id;
			v[e[rc_e].target].deg = 0;
		}
		v[rc_src].deg = 0;
	}
}

static int simple_tandem(struct asm_graph_t *g, gint_t e_i, uint32_t *comp_sz, khash_t(khInt) *lg,
		khash_t(khInt) *set_v) // for keeping visited nodes
{
	struct asm_edge_t *e = g->edges;
	struct asm_node_t *v = g->nodes;
	int j, i, n_items, cnt = 0;
	int missing = 0;
	aqueue_t *q = init_aqueue();
	uint32_t _comp_sz = 0;
	gint_t next_e, u, src, dest;
	*comp_sz = 0;

	if (e[e_i].seq_len < MIN_BRIDGE_LEG)
		return 0;

	aqueue_add(q, e_i); //add the big edge to queue
	if (e[e_i].seq_len > MIN_BRIDGE_LEG)
		kh_put(khInt, lg, e_i, &missing);
	while (1){
		n_items = q->n;
		if (n_items == 0)
			break;
		for (j = 0; j < n_items; j++) {
			u = aqueue_pop(q); //u is an edge
			tandem_helper(e, v, u, set_v, &_comp_sz, q, lg, (u == e_i));
		}
		if (++cnt > MIN_LAYER) //very complex region, never end
			break;
	}
	
	//must be visited at least 3 node and completed tour
	if (kh_size(set_v) < MIN_VISITED_NODES || cnt > MIN_LAYER)
		return 0;
	for (i = q->p; i < q->n; i++) {
		u = q->e[i]; //u is an edge
		dest = e[u].target;
		for (j = 0; j < v[dest].deg; j++) {
			if (kh_get(khInt, set_v, v[dest].adj[j]) == kh_end(set_v)){
				__VERBOSE("%d - One remain node* have outgoing edge* or region is too complex\n", e_i);
				return 0;
			}
		}
	}

	*comp_sz = _comp_sz;
	return 1;
}

static void collapse_bulge(struct asm_edge_t *e, struct asm_node_t *v, gint_t e_i, gint_t e_j,
		gint_t ei, gint_t eo)
{
	gint_t u;
	v[e[ei].target].deg = 1;
	v[e[ei].target].adj[0] = e_i;
	v[e[eo].target].deg = 1;
	v[e[eo].target].adj[0] = e[e_i].rc_id;

	REFINE_EDGE(e, v, e_i)
	if (SAME_EDGE(e, e_i, e_j))
		e[e_j].seq_len = 0;
	else{
		ISOLATE_NODE(v, e[e_j].target)
		ISOLATE_NODE(v, v[e[e_j].source].rc_id)
	}
}

static int detect_bugle(struct asm_edge_t *e, struct asm_node_t *v, gint_t e_i)
{
	gint_t dest = e[e_i].target;
	gint_t rc_src  = v[e[e_i].source].rc_id;

	if (v[dest].deg != 1 || v[rc_src].deg != 1)
		return 0;

	gint_t eo = v[dest].adj[0];
	gint_t *eo_out = v[v[e[eo].source].rc_id].adj;
	if (v[v[e[eo].source].rc_id].deg != 2)
		return 0;
	gint_t ex = e[eo_out[0]].rc_id == e_i ? e[eo_out[1]].rc_id : e[eo_out[0]].rc_id;

	gint_t ei = e[v[rc_src].adj[0]].rc_id;
	gint_t *ei_out = v[e[ei].target].adj;
	if (v[e[ei].target].deg != 2)
		return 0;
	gint_t ey = ei_out[0] == e_i ? ei_out[1] : ei_out[0];

	if (ex == ey){
		collapse_bulge(e, v, e_i, ex, ei, e[eo].rc_id);
		return 1;
	}
	return 0;
}

static int resolve_loop_wrapper(struct asm_graph_t *g, gint_t *id_edge, 
			gint_t *cc_size, khash_t(khInt) *set_v, struct asm_graph_t **g_ret)
{
	int cnt;
	struct asm_edge_t *e;
	struct asm_node_t *v;
	uint32_t ksize;
	struct asm_graph_t *g1 = calloc(1, sizeof(struct asm_graph_t));
	struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
	int i, ret;
	gint_t cc_id;
	g0 = g1 = g;
	int iter = 0;
	int flag = 1;

	do {
		asm_edge_cc(g0, id_edge, &cc_size);
		kh_clear(khInt, set_v);
		cnt = 0;
		e = g0->edges;
		v = g0->nodes;
		for(i = 0; i < g0->n_e; ++i){
			cc_id = id_edge[i];
			if (cc_size[cc_id] < MIN_COMPONENT || 
				kh_get(khInt, set_v, i) != kh_end(set_v)) //must came from large component and not be found
				continue;
			if (is_trivial_loop(e, v, i)){
				flag = 0;
				cnt++;
				__VERBOSE("Edge %d - %d\n - loop\n", i, e[i].rc_id);
				resolve_loop(g, i, g0->ksize);
				kh_put(khInt, set_v, i, &ret);
			}
		}
		if (cnt > 0)
			asm_condense(g0, g1);
		g0 = g1;
		free(cc_size);
		__VERBOSE("Number of resolved loop in " KGRN "%d iteration" KGRN " %d\n" RESET, ++iter, cnt);
	} while(cnt > 0);
	*g_ret = g1;
	return flag;
}

static int resolve_bugle_wrapper(struct asm_graph_t *g, gint_t *id_edge, 
			gint_t *cc_size, khash_t(khInt) *set_v, struct asm_graph_t **g_ret)
{
	int cnt;
	struct asm_edge_t *e;
	struct asm_node_t *v;
	uint32_t ksize;
	struct asm_graph_t *g1 = calloc(1, sizeof(struct asm_graph_t));
	struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
	int i, ret;
	gint_t cc_id;
	g0 = g1 = g;
	int iter = 0;
	int flag = 1;

	do {
		asm_edge_cc(g0, id_edge, &cc_size);
		kh_clear(khInt, set_v);
		cnt = 0;
		e = g0->edges;
		v = g0->nodes;
		for(i = 0; i < g0->n_e; ++i){
			cc_id = id_edge[i];
			if (cc_size[cc_id] < MIN_COMPONENT || 
				kh_get(khInt, set_v, i) != kh_end(set_v)) //must came from large component and not be found
				continue;
			if (detect_bugle(e, v, i)){
				cnt++;
				flag = 0;
				__VERBOSE(KGRN "Bulge %d\n" RESET, i);
				kh_put(khInt, set_v, i, &ret);
			}
		}
		if (cnt > 0)
			asm_condense(g0, g1);
		g0 = g1;
		free(cc_size);
		__VERBOSE("Number of bulge in " KGRN "%d iteration" KGRN " %d\n" RESET, ++iter, cnt);
	} while(cnt > 0);
	*g_ret = g1;
	return flag;
}

void resolve_loop_bugle(struct asm_graph_t *g0, struct asm_graph_t **g_ret)
{
	khash_t(khInt) *set_v;     // for keeping visited nodes
	khash_t(khInt) *lg;        // for keeping large contigs around the complex region 
	khash_t(khInt) *comp_set;  // for keeping edges in the complex region 
	set_v = kh_init(khInt);
	lg = kh_init(khInt);
	comp_set = kh_init(khInt);

	gint_t *id_edge, *cc_size; // for the connected component
	gint_t cc_id;
	gint_t comp_sz;
	id_edge = malloc(g0->n_e * sizeof(gint_t));
	asm_edge_cc(g0, id_edge, &cc_size);
	struct asm_graph_t *g1, *g_tmp;
	int iter = 0;
	int clean_loop, clean_bugle;

	g1 = g0;
	do {
		++iter;
		clean_loop = resolve_loop_wrapper(g1, id_edge, cc_size, set_v, &g1);
		__VERBOSE("%sClean loop!\n" RESET, clean_loop == 0?KRED:KGRN);
		clean_bugle = resolve_bugle_wrapper(g1, id_edge, cc_size, set_v, &g1);
		__VERBOSE("%sClean bugle!\n" RESET, clean_bugle == 0?KRED:KGRN);
	} while((!clean_loop || !clean_bugle) && iter < MIN_ITER_RESOLVE);

	kh_destroy(khInt, set_v);
	free(id_edge);
	free(cc_size);
	//write_fasta(g1, "step_4.fasta");
	//save_asm_graph_simple(g1, "Results/Staph_k55/graph_step4.bin");
	//write_gfa(g1, "step_4.gfa");
	*g_ret = g1;
}

static void asm_edge_cc(struct asm_graph_t *g, gint_t *id_edge, gint_t **ret_size)
{
	memset(id_edge, 255, g->n_e * sizeof(gint_t));
	gint_t m_cc = 0x10000;
	gint_t *size = malloc(m_cc * sizeof(gint_t));
	gint_t n_cc = 0;
	gint_t k, l, r;
	gint_t *q = malloc(g->n_e * sizeof(gint_t));

	for (k = 0; k < g->n_e; ++k) {
		if (id_edge[k] != -1)
			continue;
		id_edge[k] = id_edge[g->edges[k].rc_id] = n_cc;
		gint_t cur_size = 0;
		l = r = 0;
		q[0] = k;
		while (l <= r) {
			gint_t e = q[l++];
			gint_t e_rc = g->edges[e].rc_id;
			cur_size += 2 * (g->edges[e].seq_len - g->ksize);

			gint_t u, c;
			u = g->edges[e].target;
			if (g->nodes[u].deg == 0)
				cur_size += g->ksize;
			for (c = 0; c < g->nodes[u].deg; ++c) {
				gint_t next_e = g->nodes[u].adj[c];
				gint_t next_e_rc = g->edges[next_e].rc_id;
				if (id_edge[next_e] == -1) {
					id_edge[next_e] = id_edge[next_e_rc] = n_cc;
					q[++r] = next_e;
				}
			}

			u = g->edges[e_rc].target;
			if (g->nodes[u].deg == 0)
				cur_size += g->ksize;
			for (c = 0; c < g->nodes[u].deg; ++c) {
				gint_t next_e = g->nodes[u].adj[c];
				gint_t next_e_rc = g->edges[next_e].rc_id;
				if (id_edge[next_e] == -1) {
					id_edge[next_e] = id_edge[next_e_rc] = n_cc;
					q[++r] = next_e;
				}
			}
		}
		if (n_cc + 1 == m_cc) {
			m_cc <<= 1;
			size = realloc(size, m_cc * sizeof(gint_t));
		}
		size[n_cc++] = cur_size / 2;
	}
	size = realloc(size, n_cc * sizeof(gint_t));
	*ret_size = size;
	free(q);
}

//int main(int argc, char *argv[])
//{
//	struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
//	char *path = argv[1];
//	load_asm_graph(g0, path);
//	resolve_loop_bugle(g0);
//	return 0;
//}

