#include <stdlib.h>
#include <string.h>

#include "io_utils.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"
#include "queue.h"

#include "assembly_graph.h"
#include "khash.h"

#define MIN_CON_LEN 7000
#define MIN_RATIO_COV 1.7
#define MIN_BRIDGE_LEG 10000
#define MIN_COMPONENT 250
#define MAX_COMPONENT_REGION 3000
#define MIN_LAYER 100
#define MIN_VISITED_NODES 3
#define MIN_RATIO_LOOP_COV 1.5
#define MAX_MARGIN_RATIO 1.75

#define __shift_seq_i(i) ((i & 15) << 1)
#define __get_c_bin(s, i) (s[i >> 4] >> __shift_seq_i(i)) & (uint32_t)0x3
#define is_lg_leg(e, v, i, cnt) (((e[v[i].adj[0]].seq_len > MIN_BRIDGE_LEG) + \
		(e[v[i].adj[1]].seq_len > MIN_BRIDGE_LEG)) == cnt)

#define get_seq_cov(e, ksize) ((e).count / (e).seq_len - (ksize))
#define get_seq_cov(e, ksize) (e.count / (e.seq_len - ksize))

#define HAS_LOOP(e, u, t) e[u].target == e[t].source && e[u].source == e[t].target

static uint64_t g_cov;
static uint32_t n_edges;
static uint32_t n_nodes;

KHASH_SET_INIT_INT(khInt);

static void dump_b_seq(uint32_t *s, size_t l);
void isolate_edge(struct asm_edge_t *e, struct asm_node_t *v, gint_t e_i);
static uint32_t *concat_b_seq(uint32_t *dst, gint_t dlen, uint32_t *src,
					gint_t slen, int skip);

void clone_edge_rev(struct asm_edge_t *d, struct asm_edge_t *s);

static void __set_c_bin(uint32_t *s, int i, uint32_t c)
{
	s[i >> 4] = s[i >> 4] & ~((uint32_t)3 << __shift_seq_i(i)) | c << __shift_seq_i(i);
}

uint32_t *rev_bin(uint32_t *seq, size_t l)
{
	int i, n;
	uint32_t t;
	n = (l + 15) >> 4;
	uint32_t u, v;
	uint32_t *s = (uint32_t *)calloc(n, sizeof(uint32_t));
	memcpy(s, seq, n * sizeof(uint32_t));

	for (i = 0; i < l >> 1; i++) {
		u = __get_c_bin(s, i);
		v = __get_c_bin(s, l - i - 1);
		__set_c_bin(s, l - i - 1, u ^ 3);
		__set_c_bin(s, i, v ^ 3);
	}
	return s;
}

static void rev_holes(struct asm_edge_t *e)
{
	int i = 0;
	uint32_t t;
	uint32_t len = 0;
	uint32_t n = e->n_holes;
	while (i < e->n_holes)
		len += e->l_holes[i++];
	len += e->seq_len;
	for (i = 0; i < n >> 1; ++i) {
		t = e->p_holes[i];
		e->p_holes[i] = len - (e->p_holes[n - i - 1] + 1 + e->l_holes[n - i - 1]);
		e->p_holes[n - i - 1] = len - (t + 1 +  e->l_holes[i]);
		t = e->l_holes[i];
		e->l_holes[i] = e->l_holes[n - i - 1];
		e->l_holes[n - i - 1] = e->l_holes[i];
	}
}

uint32_t get_genome_cov(struct asm_graph_t *g0)
{
	struct asm_edge_t *e = g0->edges;
	int i, cnt = 0;
	uint64_t cov = 0;

	for (i = 0 ; i < g0->n_e; i++) {
		if (e[i].seq_len > MIN_CON_LEN){
			cov += get_seq_cov(e[i], g0->ksize);
			cnt++;
		}
	}
	assert(cnt > 0);
	return cov / cnt;
}


int is_trivial_loop(struct asm_edge_t *e, struct asm_node_t *v, gint_t i)
{
	if (e[i].seq == NULL)
		return 0;

	gint_t src = e[i].source;
	gint_t dest = e[i].target;
	gint_t r_src = v[src].rc_id;

	if (!(v[r_src].deg == 2 && v[dest].deg == 2))
		return 0;
	gint_t u;
	u = e[i].target;
	if (HAS_LOOP(e, i, v[u].adj[0]) || HAS_LOOP(e, i, v[u].adj[1]))
		return 1;
	return 0;
}

static int is_bridge(struct asm_edge_t *e, struct asm_node_t *v, gint_t i)
{
	gint_t src = e[i].source;
	gint_t dest = e[i].target;
	gint_t r_src = v[src].rc_id;
	if (!(v[r_src].deg == 2 && v[dest].deg == 2))
		return 0;
	if (is_lg_leg(e, v, r_src, 2) && is_lg_leg(e, v, dest, 2))
		return 1;
	return 0;
}

void glue_2_seq(struct asm_edge_t *e, struct asm_node_t *v, gint_t e1, gint_t e2,
		size_t comp_sz)
{
	gint_t u = e[e1].target;
	gint_t t = e[e2].source;
	v[v[u].rc_id].deg = 0;
	v[v[t].rc_id].deg = 0;
	v[u].deg = 0;
	v[t].deg = 0;

	e[e1].n_holes++;
	e[e1].p_holes = realloc(e[e1].p_holes, e[e1].n_holes * sizeof(uint32_t));
	e[e1].l_holes = realloc(e[e1].l_holes, e[e1].n_holes * sizeof(uint32_t));
	e[e1].p_holes[e[e1].n_holes - 1] = e[e1].seq_len - 1;
	e[e1].l_holes[e[e1].n_holes - 1] = comp_sz;

	e[e1].seq = concat_b_seq(e[e1].seq, e[e1].seq_len, e[e2].seq, e[e2].seq_len, 0);
	e[e1].seq_len = e[e1].seq_len + e[e2].seq_len;
	e[e1].count = e[e1].count + e[e2].count;

	e[e1].target = e[e2].target;
	e[e2].target = e[e2].source;
	e[e2].seq_len = 0;
	e[e2].count = 0;
	
	clone_edge_rev(e + e[e1].rc_id, e + e1);
}

void clone_edge_rev(struct asm_edge_t *d, struct asm_edge_t *s)
{
	d->seq = rev_bin(s->seq, s->seq_len);
	d->seq_len = s->seq_len;
	d->count = s->count;

	if (d->n_holes != s->n_holes){
		d->n_holes = s->n_holes;
		if (!d->p_holes)
			d->p_holes = malloc(d->n_holes * sizeof(uint32_t));
		if (!d->l_holes)
			d->l_holes = malloc(d->n_holes * sizeof(uint32_t));
	}
	memcpy(d->p_holes, s->p_holes, d->n_holes * sizeof(uint32_t));
	memcpy(d->l_holes, s->l_holes, d->n_holes * sizeof(uint32_t));
	rev_holes(d);
}

void isolate_edge(struct asm_edge_t *e, struct asm_node_t *v, gint_t e_i)
{
	gint_t dst = e[e_i].target;
	gint_t src = e[e_i].source;

	v[dst].deg = 0;
	v[src].deg = 0;
	v[v[dst].rc_id].deg = 0;
	v[v[src].rc_id].deg = 0;
	e[e_i].seq_len = 0;
	e[e[e_i].rc_id].seq_len = 0;
}

static stretch_edge(struct asm_edge_t *e, struct asm_node_t *v, gint_t e_i,
		gint_t ei, gint_t eo)
{
	int k;
	if (e_i == 673518)
		k = 1;
	gint_t dst = e[e_i].target;
	gint_t rc_src = v[e[e_i].source].rc_id;
	gint_t rc_dest = v[e[e_i].target].rc_id;
	gint_t src = e[e_i].source;

	v[dst].deg = 1;
	v[dst].adj[0] = eo;
	v[v[dst].rc_id].deg = 1;
	v[v[dst].rc_id].adj[0] = e[e_i].rc_id;

	v[src].deg = 1;
	v[src].adj[0] = e_i;
	v[rc_src].deg = 1;
	v[rc_src].adj[0] = e[ei].rc_id;
}

uint32_t *concat_b_seq(uint32_t *dst, gint_t dlen, uint32_t *src,
			gint_t slen, int skip)
{
	{
		gint_t i;
		for (i = 0; i < skip; ++i) {
			uint32_t c1, c2;
			c1 = __get_c_bin(src, i);
			c2 = __get_c_bin(dst, dlen - (skip - i));
			if (c1 != c2) {
				dump_b_seq(dst, dlen);
				dump_b_seq(src, slen + skip);
				assert(0 && "Fail when append");
			}
		}
	}
	uint32_t *new_ptr;
	gint_t cur_m, m, i, k;
	cur_m = (dlen + 15) >> 4;
	m = (dlen + slen + 15) >> 4;
	if (m > cur_m) {
		new_ptr = (uint32_t *)malloc(m * sizeof(uint32_t));
		memcpy(new_ptr, dst, ((dlen + 15) >> 4) * sizeof(uint32_t));
		memset(new_ptr + cur_m, 0, (m - cur_m) * sizeof(uint32_t));
	} else {
		new_ptr = dst;
	}
		uint32_t u, v;
	for (i = skip; i < skip + slen; ++i) {
		k = dlen + (i - skip);
		u = __get_c_bin(src, i);
		__set_c_bin(new_ptr, k , u);
		v = __get_c_bin(new_ptr, k);
		assert(v == u);	//check assign successfully
		}
	return new_ptr;
}

void concat_3_seq(struct asm_edge_t e1, struct asm_edge_t e2,
				struct asm_edge_t e3, gint_t ksize, gint_t *len, uint32_t **seq)
{
	*seq = concat_b_seq(e1.seq, e1.seq_len, e2.seq, e2.seq_len - ksize, ksize);
	*seq = concat_b_seq(*seq , e1.seq_len + e2.seq_len - ksize, e3.seq, e3.seq_len - ksize, ksize); 
	*len = e1.seq_len + e2.seq_len + e3.seq_len - 2*ksize;
}

static void condense_3_seq(gint_t e1, gint_t e2, gint_t e3, struct asm_edge_t *e, 
			struct asm_node_t *v, gint_t ksize)
{
	v[e[e1].target].deg = 0;
	v[e[e3].source].deg = 0;
	uint32_t *seq;
	gint_t len;

	concat_3_seq(e[e1], e[e2], e[e3], ksize, &len, &seq);
	e[e2].count = e[e1].count + e[e2].count + e[e3].count;
	e[e2].seq = seq;
	e[e2].seq_len = len;
	clone_edge_rev(e + e[e2].rc_id, e + e2);
	__VERBOSE("Condensed length of %d: %d\n", e2, len);

	e[e2].source = e[e1].source;
	e[e2].target = e[e3].target;

	e[e1].source = v[e[e1].target].rc_id;
	e[e3].target = v[e[e3].source].rc_id;

	e[e1].seq_len = 0;
	e[e[e1].rc_id].seq_len = 0;
	e[e3].seq_len = 0;
	e[e[e3].rc_id].seq_len = 0;

}

int resolve_loop(struct asm_graph_t *g, gint_t u,
		gint_t ksize)
{
	struct asm_edge_t *e = g->edges;
	struct asm_node_t *v = g->nodes;
	uint32_t *seq;
	gint_t t, dest;
	gint_t rc_src;
	gint_t len;
	gint_t e_in, e_out;

	dest = e[u].target;
	rc_src = v[e[u].source].rc_id;

	t = HAS_LOOP(e, u, v[dest].adj[0]) ? v[dest].adj[0] : v[dest].adj[1];

	__VERBOSE("Concat the loop and the shared\n");

	concat_3_seq(e[u], e[t], e[u], ksize, &len, &seq);
	e[u].count += e[t].count + e[u].count;
	e[u].seq_len = len;
	free(e[u].seq);
	e[u].seq = seq;
 
	clone_edge_rev(e + e[u].rc_id, e + u);
	asm_remove_edge(g, t);
	asm_remove_edge(g, e[t].rc_id);
}
/*
static void find_forest(struct asm_graph_t *g0)
{
	init_clock();
	g_cov = get_genome_cov(g0);

	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips\n");
	__VERBOSE_LOG("kmer_%d_graph_#1", "Number of nodes: %lld\n", g0->ksize,
							(long long)g0->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#1", "Number of edges: %lld\n", g0->ksize,
							(long long)g0->n_e);
	__VERBOSE("Genome walk coverage %d\n", (int)g_cov);

	n_edges = g0->n_e;
	n_nodes = g0->n_v;

	int i;
	int ret, flag;
	gint_t r_src;
	gint_t src, dest;
	gint_t lg;
	struct asm_edge_t *e = g0->edges;
	struct asm_node_t *v = g0->nodes;

	khash_t(khInt) *set_v; // for keeping visited nodes
	set_v = kh_init(khInt);

	gint_t *id_edge, *cc_size; // for the connected component
	gint_t cc_id;
	gint_t comp_sz;
	id_edge = malloc(g0->n_e * sizeof(gint_t));
	cc_size = NULL;
	asm_edge_cc(g0, id_edge, &cc_size);

	//MAIN STUFF HERE
	for (i = 0 ; i < g0->n_e; i++) {
		cc_id = id_edge[i];
		if (cc_size[cc_id] < MIN_COMPONENT || kh_get(khInt, set_v, i) != kh_end(set_v)) //must came from large component and not be found
			continue;
		src = e[i].source;
		dest = e[i].target;
		r_src = v[src].rc_id;
		flag = 0;
		if (is_bridge(e, v, i)){
			__VERBOSE("Edge %d - %d\n - bridge\n", i, e[i].rc_id);
			flag = 1;
		}
		if (is_trivial_loop(e, v, i)){
			__VERBOSE("Edge %d - %d\n - loop\n", i, e[i].rc_id);
			resolve_loop(e, v, i, g0->ksize);
			flag = 1;
		}
		if (is_simple_tandem(e, v, i, &lg, &comp_sz)){
			__VERBOSE("Edge %d - %d\n - simple complex\n", i, e[i].rc_id);
			glue_2_seq(e, v, i, lg, comp_sz);
			flag = 1;
		}
		if (flag){
			kh_put(khInt, set_v, i, &ret);
			kh_put(khInt, set_v, e[i].rc_id, &ret);
		}
	}

	write_gfa(g0, "resolved_loop.gfa");
	dump_fasta(g0, "resolved_loop.fasta");
	kh_destroy(khInt, set_v);
	free(id_edge);

} */

static void dump_b_seq(uint32_t *s, size_t l)
{
	extern char *nt4_char;
	size_t i;
	for(i = 0; i < l; i++) {
	 putc(nt4_char[__get_c_bin(s, i)], stderr);
	}
	putc('\n', stderr);
}

//int main(int argc, char *argv[])
//{
//	struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
//	char *path = argv[1];
//	load_asm_graph(g0, path);
//	find_forest(g0);
//
//	return 0;
//}
