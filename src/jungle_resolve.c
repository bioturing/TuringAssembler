#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define RESET "\x1B[0m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"


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
#define MAX_NUMBER_LEGS 9
#define UNIT_GENOME_COV_LOW 0.51
#define UNIT_GENOME_COV_HIGH 1.51


KHASH_SET_INIT_INT(khInt)

// https://github.com/attractivechaos/klib/issues/49
// shorthand way to get the key from hashtable or defVal if not found
#define kh_get_val(kname, hash, key, defVal) ({k=kh_get(kname, hash, key);(k!=kh_end(hash)?kh_val(hash,k):defVal);})

// shorthand way to set value in hash with single line command.  Returns value
// returns 0=replaced existing item, 1=bucket empty (new key), 2-adding element previously deleted
#define kh_set(kname, hash, key, val) ({int ret; k = kh_put(kname, hash,key,&ret); kh_value(hash,k) = val; ret;})

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
				log_info(KRED "Self loop detect e.i is both source and sink %d\n" RESET, key);
				return 1;
			}
		}
	}
	return 0;
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

void tandem_helper(struct asm_edge_t *e, struct asm_node_t *v,
			  gint_t u, khash_t(khInt) *set_v, uint32_t *comp_sz,
			  aqueue_t *q, khash_t(khInt) *lg,
			  int is_first)
{
	gint_t next_e;
	int i, missing;
	gint_t rc_src, dest;
	rc_src = v[e[u].source].rc_id;
	dest = e[u].target;
	khiter_t k;

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
			next_e = v[rc_src].adj[i]; //iterate outgoing edge of the target of u
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

int simple_tandem(struct asm_graph_t *g, gint_t e_i, uint32_t *comp_sz, khash_t(khInt) *lg,
		khash_t(khInt) *set_v, gint_t *is_visited) // for keeping visited nodes
{
	assert(kh_size(set_v) == 0);
	struct asm_edge_t *e = g->edges;
	struct asm_node_t *v = g->nodes;
	int j, i, n_items, cnt = 0;
	int missing = 0;
	aqueue_t *q = init_aqueue();
	uint32_t _comp_sz = 0;
	gint_t next_e, u, src, dest;
	khiter_t k;
	*comp_sz = 0;

	if (e[e_i].seq_len < MIN_BRIDGE_LEG)
		return 0;

	aqueue_add(q, e_i); //add the big edge to queue
	if (e[e_i].seq_len > MIN_BRIDGE_LEG){
		k = kh_put(khInt, lg, e[e_i].rc_id, &missing);
	}

	while (1){
		n_items = q->n;
		if (n_items == 0)
			break;
		for (j = 0; j < n_items; j++) {
			u = aqueue_pop(q); //u is an edge
			if (u != e_i && is_visited[u] == 1)
				return 0;
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
				log_info("%d - One remain node* have outgoing edge* or region is too complex\n", e_i);
				return 0;
			}
		}
	}

	*comp_sz = _comp_sz;
	return 1;
}

#define glue_seq_process(g, e1, e2) \
	g->edges[e1].l_holes = realloc(g->edges[e1].l_holes,	/* one more for gap */ \
				(g->edges[e1].n_holes + 1) * sizeof(uint32_t)); \
	g->edges[e1].p_holes = realloc(g->edges[e1].p_holes,	/* one more for gap */ \
				(g->edges[e1].n_holes + 1) * sizeof(uint32_t)); \
	g->edges[e1].l_holes[g->edges[e1].n_holes] = gap_size; \
	g->edges[e1].p_holes[g->edges[e1].n_holes] = g->edges[e1].seq_len; \
	++g->edges[e1].n_holes; \
\
	asm_append_edge_seq(g->edges + e1, g->edges + e2, g->ksize); \
	g->edges[e1].count += g->edges[e2].count; \
	g->edges[e1].target = g->edges[e2].target;

#define glue_2seq_procedure(g, e1, e2, e_rc1, e_rc2) \
	glue_seq_process(g, e1, e2); \
	glue_seq_process(g, e_rc2, e_rc1); \
	\
	g->edges[e1].rc_id = e_rc2; \
	g->edges[e_rc2].rc_id = e1; \
	asm_remove_edge(g, e_rc1); \
	asm_remove_edge(g, e2);

int resolve_jungle(struct asm_graph_t *g, khash_t(khInt) *h,
			khash_t(khInt) *comp_set, float gcov)
{
	uint32_t gap_size = 0;
	khiter_t k;
	for (k = kh_begin(comp_set); k != kh_end(comp_set); ++k) {
		if (!kh_exist(comp_set, k))
			continue;
		gint_t e = kh_key(comp_set, k);
		if (g->edges[e].source == -1)
			continue;
		gint_t len = get_edge_len(g->edges + e);
		int cov = (int)(__get_edge_cov(g->edges + e, g->ksize) / gcov + 0.499999);
		gap_size += cov * (len - g->ksize);
		gint_t e_rc = g->edges[e].rc_id;
		/* Remove edges , e.i isolate the nodes */
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		fprintf(stderr, "Remove edge %ld-%ld\n", e, e_rc);
	}
	gint_t *tmp = alloca(10 * sizeof(gint_t));
	int j = 0;
	for (k = kh_begin(h); k != kh_end(h); ++k) {
		if (kh_exist(h, k))
			tmp[j++] = kh_key(h, k);
	}
	assert(j == 2);
	gint_t e1, e2, e_rc1, e_rc2;
	e2 = tmp[1];
	e_rc2 = g->edges[e2].rc_id;
	e_rc1 = tmp[0];
	e1 = g->edges[e_rc1].rc_id;

	if (g->edges[e1].source == -1 || e1 == e2 || e1 == e_rc2)
		return 0;
	
	glue_2seq_procedure(g, e1, e2, e_rc1, e_rc2);
	
	fprintf(stderr, "Merge edges %ld - %ld\n", e1, e2);
	return 1;
}

int resolve_jungle4(struct asm_graph_t *g, khash_t(khInt) *h,
			khash_t(khInt) *comp_set, float gcov,
			gint_t *bx_bin)
{
	uint32_t gap_size = 0;
	khiter_t k;

	/* Get the big legs */
	gint_t *tmp = alloca(10 * sizeof(gint_t));
	int j = 0, i;
	for (k = kh_begin(h); k != kh_end(h); ++k) {
		if (kh_exist(h, k))
			tmp[j++] = kh_key(h, k);
	}
	assert(j == 4); // Must be 4 legs
	/* Testing for the significant of barcode */
	int cnt = 0;
	int x, y;
	uint32_t max = 0, score;

	for (j = 0; j < 4; ++j) {
		for (i = j + 1; i < 4; ++i) {
			int ret = test_edge_barcode2(g, tmp[j], tmp[i], bx_bin, &score);
			log_info("Test connection %ld <-> %ld: %s%d\n" RESET,
				tmp[j], tmp[i], ret == 1?KGRN:KRED, ret);
			//print_test_barcode_edge(g, tmp[j], tmp[i]);
			if (ret == 1) {
				if (score > max) {
					x = i;
					y = j;
					max = score;
				}
				++cnt;
			}
		}
	}

	/* Must be two ins and two outs */
	if (cnt == 0)
		return 0;

	/* Clean the complex */
	for (k = kh_begin(comp_set); k != kh_end(comp_set); ++k) {
		if (!kh_exist(comp_set, k))
			continue;
		gint_t e = kh_key(comp_set, k);
		if (g->edges[e].source == -1)
			continue;
		gint_t len = get_edge_len(g->edges + e);
		int cov = (int)(__get_edge_cov(g->edges + e, g->ksize) / gcov + 0.499999);
		gap_size += cov * (len - g->ksize);
		gint_t e_rc = g->edges[e].rc_id;
		asm_remove_edge(g, e);
		asm_remove_edge(g, e_rc);
		fprintf(stderr, "Remove edge %ld-%ld\n", e, e_rc);
	}

		
	/* Glue process */
	gint_t e1, e2, e_rc1, e_rc2;
	int n_iter = 0;

	do {
		e2 = tmp[x];
		e_rc2 = g->edges[e2].rc_id;
		e_rc1 = tmp[y];
		e1 = g->edges[e_rc1].rc_id;
		
		log_info(KCYN "One pair %d - %d\n" RESET, e1, e2);
		if (g->edges[e1].source == -1 || e1 == e2 || e1 == e_rc2)
			return 0;
		
		glue_2seq_procedure(g, e1, e2, e_rc1, e_rc2);

		if ((x + y) == 3) {
			x ^= 1;
			y ^= 1;
		} else {
			x = 3 - x;
			y = 3 - y;
		}
		++n_iter;
	} while (n_iter < 2);
	return 1;
}

static inline gint_t find_adj_idx(gint_t *adj, gint_t deg, gint_t id)
{
	gint_t i, ret;
	ret = -1;
	for (i = 0; i < deg; ++i) {
		if (adj[i] == id)
			ret = i;
	}
	return ret;
}

static void resolve_baby_flow(struct asm_graph_t *g, gint_t e, float gcov) 
{
	gint_t src = g->edges[e].source;
	gint_t ei_rc = g->nodes[g->nodes[src].rc_id].adj[0];
	int cov_i = (int)(__get_edge_cov(g->edges + ei_rc, g->ksize) / gcov + 0.499999);
	int j, i = find_adj_idx(g->nodes[src].adj, g->nodes[src].deg, e);
	assert(i >= 0);
	for (j = 0; j < g->nodes[src].deg; ++j){
		int cov_o = (int)(__get_edge_cov(g->edges + g->nodes[src].adj[j], g->ksize) / gcov + 0.499999);
		if (cov_o == 1 && cov_i == 1){
			log_info(KCYN "Remove the baby %d\n " RESET, e);
			asm_remove_edge(g, e);
			asm_remove_edge(g, g->edges[e].rc_id);
			break;
		}
	}
}

int jungle_resolve_flow(struct asm_graph_t *g, khash_t(khInt) *h,
			khash_t(khInt) *comp_set, float gcov)
{
	khiter_t k;
	/* Find the baby */
	for (k = kh_begin(comp_set); k != kh_end(comp_set); ++k) {
		if (!kh_exist(comp_set, k))
			continue;
		gint_t e = kh_key(comp_set, k);
		if (g->edges[e].source == -1) //FIXME: I think we don't need this in this case
			continue;
		gint_t len = get_edge_len(g->edges + e);
		int cov = (int)(__get_edge_cov(g->edges + e, g->ksize) / gcov + 0.499999);
		if (!cov){
			log_info(KCYN "Found the baby %d\n " RESET, e);
			resolve_baby_flow(g, e, gcov);
		}
	}
}


#define set_max(best_e, max, e, score) if (score > max){ best_e = e; max = score; } 
static gint_t share_barcode_profile(struct asm_graph_t *g, gint_t e_i, khash_t(khInt) *comp_set,
		gint_t *bx_bin)
{
	khiter_t k;
	int ret[2];
	khash_t(khInt) *set;
	set = kh_init(khInt);
	int missing;
	gint_t best_e;
	uint32_t max = 0;
	uint32_t score;

	for (k = kh_begin(comp_set); k != kh_end(comp_set); ++k) {
		if (!kh_exist(comp_set, k))
			continue;
		gint_t key = kh_key(comp_set, k);
		if (kh_get(khInt, set, key) != kh_end(set))
			continue;
		if (g->edges[key].source == -1) 
			continue;
		if (key == e_i || key == g->edges[e_i].rc_id)
			continue;

		ret[0] = test_edge_barcode2(g, e_i, key, bx_bin, &score);
		if (ret[0])
			set_max(best_e, max, key, score);
		ret[1] = test_edge_barcode2(g, e_i, g->edges[key].rc_id, bx_bin, &score);
		if (ret[1])
			set_max(best_e, max, g->edges[key].rc_id, score);
		kh_put(khInt, set, key, &missing);
		kh_put(khInt, set, g->edges[key].rc_id, &missing);
	}
	if (max >= 7){
		return best_e;
	} else{
		return -1;
	}
}

static int asm_append_pair(struct asm_graph_t *g0, gint_t e1, gint_t e2)
{
	asm_append_edge2(g0, e1, g0->edges[e2].rc_id);
	asm_append_edge2(g0, e2, g0->edges[e1].rc_id);
	asm_remove_edge(g0, g0->edges[e2].rc_id);
	asm_remove_edge(g0, g0->edges[e1].rc_id);
	g0->edges[e1].rc_id = e2;
	g0->edges[e2].rc_id = e1;
	assert(g0->edges[e1].seq != NULL);
	assert(g0->edges[e2].seq != NULL);
}

static void asm_append_triple(struct asm_graph_t *g0, gint_t e1, gint_t e2, gint_t e3)
{
	asm_append_edge2(g0, e2, e3);
	asm_append_edge2(g0, g0->edges[e3].rc_id, g0->edges[e2].rc_id);
	g0->edges[e2].rc_id = g0->edges[e3].rc_id;
	g0->edges[g0->edges[e3].rc_id].rc_id = e2;
	asm_remove_edge(g0, e3);
	asm_remove_edge(g0, g0->edges[e2].rc_id);

	asm_append_edge2(g0, e1, e2);
	asm_append_edge2(g0, g0->edges[e2].rc_id, g0->edges[e1].rc_id);
	g0->edges[e1].rc_id = g0->edges[e2].rc_id;
	g0->edges[g0->edges[e2].rc_id].rc_id = e1;
	asm_remove_edge(g0, e2);
	asm_remove_edge(g0, g0->edges[e1].rc_id);
}

static void set_visited_edge(struct asm_graph_t *g, gint_t *is_visited, khash_t(khInt) *h)
{
	khiter_t k;
	gint_t key;
	for (k = kh_begin(h); k != kh_end(h); ++k){
		if (kh_exist(h, k)){
			key = kh_key(h, k);
			assert(key < g->n_e);
			is_visited[key] = 1;
			is_visited[g->edges[key].rc_id] = 1;
		}
	}
}

static void polish_link(struct asm_graph_t *g0, gint_t *fw_link, gint_t *rv_link)
{
	int i;
	gint_t tmp;
	for (i = 0; i < g0->n_e; ++i){
		if (fw_link[i] == fw_link[g0->edges[i].rc_id] && fw_link[i] != 0){
			log_info(KGRN "Duplicated %d - %d - %d\n" RESET, i, g0->edges[i].rc_id,
					fw_link[i]);
			if (fw_link[fw_link[i]] == i){
				// i is more confident than rev(i)
				// so find which node is linked to rev(i), then assign for the foward link of rev(i)
				fw_link[g0->edges[i].rc_id] = rv_link[g0->edges[i].rc_id];
			}
			if (fw_link[fw_link[i]] == g0->edges[i].rc_id){
				// the same as above
				fw_link[i] = rv_link[i];
			}
		}
	}
}

static void resolve_one_complex(struct asm_graph_t *g0, khash_t(khInt) *lg, gint_t *bx_bin,
		gint_t *fw_link, gint_t *rv_link)
{
	khiter_t k;
	gint_t key;
	int ret, i;
	int flag;

	for (k = kh_begin(lg); k != kh_end(lg); ++k){
		if (kh_exist(lg, k)){
			key = kh_key(lg, k);
			if (key == 111310 || key == 111311)
				flag = 1;
			ret = share_barcode_profile(g0, key, lg, bx_bin);
			if (ret != -1){
				fw_link[key] = ret; // WARN: cannot to be link here
				rv_link[ret] = key; // WARN: cannot to be link here
			}
		}
	}
}

static void scaffold_one_step(struct asm_graph_t *g0, gint_t *fw_link, 
		gint_t ei, khash_t(khInt) *set)
{
	int i, ret;
	gint_t e_it = fw_link[ei];
	gint_t tmp;
	printf("%d-", ei);
	kh_put(khInt, set, ei, &ret);
	while (e_it != 0) {
		printf("%d-", e_it);
		log_info(KCYN "%d - %d\n" RESET, ei, e_it);
		tmp = fw_link[g0->edges[e_it].rc_id];
		//asm_append_pair(g0, ei, g0->edges[e_it].rc_id);
		e_it = tmp;
		kh_put(khInt, set, e_it, &ret);
	}
	printf("\n");
}

static void iterative_scaffolding(struct asm_graph_t *g0, gint_t *fw_link)
{
	khash_t(khInt) *set;
	set = kh_init(khInt);
	int i;
	for (i = 0; i < g0->n_e; ++i){
		if (fw_link[i] == 0 && fw_link[g0->edges[i].rc_id] != 0){
			log_info(KCYN "Found %d\n" RESET, i);
			if (kh_get(khInt, set, g0->edges[i].rc_id) == kh_end(set)){
				scaffold_one_step(g0, fw_link, g0->edges[i].rc_id, set);
			}
		}
	}
}

void detect_simple_tandem(struct asm_graph_t *g0)
{
	float gcov = get_genome_coverage(g0);
	khash_t(khInt) *set_v;     // for keeping visited nodes
	khash_t(khInt) *lg;        // for keeping large contigs around the complex region 
	khash_t(khInt) *comp_set;  // for keeping edges in the complex region 
	int i, ret, missing;
	set_v = kh_init(khInt);
	lg = kh_init(khInt);
	comp_set = kh_init(khInt);

	gint_t *bx_bin;            // To store number of bin of each edge before merging
	gint_t *id_edge, *cc_size; // for the connected component
	gint_t cc_id;
	uint32_t comp_sz;
	id_edge = malloc(g0->n_e * sizeof(gint_t));
	gint_t *is_visited = calloc(g0->n_e, sizeof(gint_t));
	bx_bin = calloc(g0->n_e, sizeof(gint_t));
	memset(is_visited, 0, g0->n_e * sizeof(gint_t));
	memset(bx_bin, 0, g0->n_e * sizeof(gint_t));
	barcode_bin_profiling(g0, bx_bin);
	asm_edge_cc(g0, id_edge, &cc_size);

	gint_t *fw_link = calloc(g0->n_e, sizeof(gint_t));
	memset(fw_link, 0, sizeof(gint_t) * g0->n_e);
	gint_t *rv_link = calloc(g0->n_e, sizeof(gint_t));
	memset(rv_link, 0, sizeof(gint_t) * g0->n_e);

	int cnt = 0;
	
	for (i = 0; i < g0->n_e; ++i){
		kh_clear(khInt, lg);
		kh_clear(khInt, comp_set);
		cc_id = id_edge[i];
		if (cc_size[cc_id] < MIN_COMPONENT || 
			kh_get(khInt, set_v, i) != kh_end(set_v)) //must came from large component and not be found
			continue;

		if (simple_tandem(g0, i, &comp_sz, lg, comp_set, is_visited)){
			set_visited_edge(g0, is_visited, comp_set);
			log_info(KGRN "Complex Tandem %d\n" RESET, i);
			log_info(KBLU "Numer of keys %d\n" RESET, kh_size(lg));
			log_info(KMAG "Numer of edges in the complex %d\n" RESET, kh_size(comp_set));
			log_info(KWHT "Size of the complex %d\n" RESET, comp_sz);
			if (kh_size(lg) <= MAX_NUMBER_LEGS){
				if (is_large_loop(g0, lg, comp_set)){
					log_info(KRED "Is a self loop complex\n" RESET);
					jungle_resolve_flow(g0, lg, comp_set, gcov);
				}
				/* resolve 1-1 complex */
				if (kh_size(lg) == 2)
					cnt += resolve_jungle(g0, lg, comp_set, gcov);
				else if (kh_size(lg) == 4) {
					cnt += resolve_jungle4(g0, lg, comp_set, gcov, bx_bin);
				}
			}
			kh_put(khInt, set_v, i, &missing);
		}
	}

	kh_clear(khInt, set_v);
	memset(is_visited, 0, g0->n_e * sizeof(gint_t));
	for (i = 0; i < g0->n_e; ++i){
		kh_clear(khInt, lg);
		kh_clear(khInt, comp_set);
		cc_id = id_edge[i];
		if (cc_size[cc_id] < MIN_COMPONENT || 
			kh_get(khInt, set_v, i) != kh_end(set_v)) //must came from large component and not be found
			continue;

		if (simple_tandem(g0, i, &comp_sz, lg, comp_set, is_visited)){
			set_visited_edge(g0, is_visited, comp_set);
			log_info(KGRN "Complex Tandem %d\n" RESET, i);
			log_info(KBLU "Numer of keys %d\n" RESET, kh_size(lg));
			log_info(KMAG "Numer of edges in the complex %d\n" RESET, kh_size(comp_set));
			log_info(KWHT "Size of the complex %d\n" RESET, comp_sz);
			resolve_one_complex(g0, lg, bx_bin, fw_link, rv_link);
			kh_put(khInt, set_v, i, &missing);
		}
	}
	polish_link(g0, fw_link, rv_link);
	iterative_scaffolding(g0, fw_link);
	log_info("Number of resolved jungle: %d\n", cnt);
}

