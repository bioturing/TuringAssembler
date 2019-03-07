#include <stdlib.h>
#include <string.h>

#include "simple.h"

#include "io_utils.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"
#include "queue.h"

#define MIN_CON_LEN 7000
#define MIN_RATIO_COV 1.7
#define MIN_BRIDGE_LEG 3000
#define MIN_COMPONENT 250
#define MAX_COMPONENT_REGION 3000
#define MIN_LAYER 100
#define MIN_VISITED_NODES 3
#define MIN_RATIO_LOOP_COV 1.5


static uint64_t g_cov;
static uint32_t n_edges;
static uint32_t n_nodes;


uint32_t get_seq_cov(struct asm_edge_t *e, int ksize)
{
  uint32_t cov = e->count / (e->seq_len + 2 * ksize);
  return cov;
}

uint32_t get_genome_cov(struct asm_graph_t *g0)
{
  struct asm_edge_t *e = g0->edges;
  int i, cnt = 0;
  uint64_t cov = 0;

  for (i = 0 ; i < g0->n_e; i++){
    if (e[i].seq_len > MIN_CON_LEN){
      cov += get_seq_cov(e + i, g0->ksize);
      cnt++;
    }
  }
  assert(cnt > 0);
  return cov / cnt;
}

int is_lg_leg(struct asm_edge_t *e, struct asm_node_t *v , gint_t i, int cnt)
{
  assert(v[i].deg == 2);
  return (((e[v[i].adj[0]].seq_len > MIN_BRIDGE_LEG) + 
          (e[v[i].adj[1]].seq_len > MIN_BRIDGE_LEG)) == cnt);
}

/*            ______
 *           /      \
 *          / self   \
 *          \ loop   /
 *           \______/
 *___________/      \_____________
 *
 */
int is_trivial_loop(struct asm_edge_t *e, struct asm_node_t *v, gint_t i)
{
  gint_t src = e[i].source;
  gint_t dest = e[i].target;
  gint_t r_src = v[src].rc_id;

  if (!(v[r_src].deg == 2 && v[dest].deg == 2))
    return 0;
  if (!is_lg_leg(e, v, r_src, 1) || !is_lg_leg(e, v, dest, 1))// needed condition
    return 0;

  gint_t u;
  u = e[i].target;
  gint_t j = e[v[u].adj[0]].seq_len > MIN_BRIDGE_LEG ? v[u].adj[1] : v[u].adj[0];
  if (e[j].target == e[i].source && e[i].source == e[j].target)
    return 1;
  return 0;
}

/*     \                         /
 *      \                       /  
 *       \                     / 
 *        \_______bridge______/
 *        /                   \
 *       /                     \
 *      /                       \
 *     /                         \
 */                           
int is_bridge(struct asm_edge_t *e, struct asm_node_t *v, gint_t i)
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


void simple_tandem_helper(struct asm_edge_t *e, struct asm_node_t *v,
                          gint_t u, khash_t(khInt) *set_v, uint32_t *comp_sz,
                          aqueue_t *q, int *n_larges, gint_t *lg)
{
  gint_t next_e;
  int i, missing;
  gint_t src, dest;
  src = e[u].source;
  dest = e[u].target;

  assert(src < n_nodes);
  assert(dest < n_nodes);

  for (i = 0 ; i < v[dest].deg; i++){
    next_e = v[dest].adj[i]; //iterate outgoing edge of the target of u
    if (kh_get(khInt, set_v, next_e) == kh_end(set_v)){ // visited or not
      if (e[next_e].seq_len > MIN_BRIDGE_LEG){
         (*n_larges)++;
         *lg = next_e;
         continue;
      }
      assert(next_e <= n_edges);
      aqueue_add(q, next_e); //add neighbor edge
      kh_put(khInt, set_v, next_e, &missing);
      *comp_sz += e[next_e].seq_len; //add size of the edge to total size
    }
  }
}

/*              ___________
 *             /           \
 *            /    /\       \
 *___________/____/  \______ \_____________
 *           \    \____  /   /
 *            \        \/   /
 *             \___________/
 *
 */

int is_simple_tandem(struct asm_edge_t *e, struct asm_node_t *v, gint_t e_i,
                     gint_t *k)
{
  int j, i, n_items, cnt = 0;
  int n_larges = 0;
  int missing = 0;
  khash_t(khInt) *set_v; // for keeping visited nodes
  set_v = kh_init(khInt);
  aqueue_t *q = init_aqueue();
  uint32_t comp_sz = 0;
  gint_t next_e, u, src, dest, lg;

  *k = NULL;

  if (e[e_i].seq_len < MIN_BRIDGE_LEG)
    return 0;

  aqueue_add(q, e_i); //add the big edge to queue
  while (1){
    n_items = q->n;
    if (n_items == 0)
      break;
    for (j = 0; j < n_items; j++){
      u = aqueue_pop(q); //u is an edge
      if (u >= n_edges){
        printf("Here\n");
      }
      assert(u < n_edges);
      simple_tandem_helper(e, v, u, set_v, &comp_sz, q, &n_larges, &lg);
    }
    if (++cnt > MIN_LAYER) //very complex region, never end
      break;
  }
  
  //must be visited at least 3 node and completed tour
  if (kh_size(set_v) < MIN_VISITED_NODES || cnt > MIN_LAYER)
    return 0;
  for (i = q->p; i < q->n; i++){
    u = q->e[i]; //u is an edge
    dest = e[u].target;
    for (j = 0; j < v[dest].deg; j++){
      if (kh_get(khInt, set_v, v[dest].adj[j]) == kh_end(set_v)){
        __VERBOSE("%d - One remain node* have outgoing edge* or region is too complex\n", e_i);
        return 0;
      }
    }
  }

  if (n_larges > 1 || n_larges == 0)
    return 0;
  else if (n_larges == 1){
    *k = lg;
    return 1;
  }
}

void resolve_loop(struct asm_edge_t *e, struct asm_node_t *v, gint_t u,
                  gint_t ksize)
{
  uint32_t cov, t_cov;
  uint32_t *seq;
  gint_t t, dest;
  gint_t rc_src;

  cov = get_seq_cov(e + u, ksize);
  dest = e[u].target;
  rc_src = v[e[u].source].rc_id;

  t = e[v[dest].adj[0]].seq_len >= MIN_BRIDGE_LEG ? v[dest].adj[1]:v[dest].adj[0];
  t_cov = get_seq_cov(e + t, ksize);
  
  if (!(cov/t_cov < MIN_RATIO_COV && (uint32_t)(cov/g_cov) == 2))
    __VERBOSE("Loop is not double of the genome walk\n");
  
  e[u].seq_len = 2*e[u].seq_len + e[t].seq_len - 2*ksize;
  e[u].count += e[t].count;
  seq = append_bin_seq(e[t].seq, e[t].seq_len, e[u], e[u].seq_len, ksize);
  seq = append_bin_seq(e[u].seq, e[u].seq_len, seq , e[u].seq_len + e[t].seq_len - ksize, ksize); 

  e[u].seq = seq;
  e[e[u].rc_id].seq = seq;

  t = t == v[dest].adj[1] ? v[dest].adj[0] : v[dest].adj[1]; //t is now the big leg;
  v[dest].deg = 1;
  v[dest].adj[0] = t;

  t = e[v[rc_src].adj[0]].seq_len >= MIN_BRIDGE_LEG ? v[dest].adj[1]:v[dest].adj[0];

}

void find_forest(struct asm_graph_t *g0)
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
  int ret, flag, loop;
  gint_t r_src;
  gint_t src, dest;
  gint_t lg;
  struct asm_edge_t *e = g0->edges;
  struct asm_node_t *v = g0->nodes;

  khash_t(khInt) *set_v; // for keeping visited nodes
  set_v = kh_init(khInt);

  gint_t *id_edge, *cc_size; // for the connected component
  gint_t cc_id;
  id_edge = malloc(g0->n_e * sizeof(gint_t));
  cc_size = NULL;
  asm_edge_cc(g0, id_edge, &cc_size);

  /* MAIN STUFF HERE */
  for (i = 0 ; i < g0->n_e; i++){
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
      flag = 1;
    }
    if (is_simple_tandem(e, v, i, &lg)){
      __VERBOSE("Edge %d - %d\n - simple complex\n", i, e[i].rc_id);
      flag = 1;
    }
    if (flag){
      kh_put(khInt, set_v, i, &ret);
      kh_put(khInt, set_v, e[i].rc_id, &ret);
    }
  }
  kh_destroy(khInt, set_v);
  free(id_edge);

}


int main(int argc, char *argv[])
{
  struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
  char *path = argv[1];
  load_asm_graph(g0, path);
  gint_t k;
  find_forest(g0);
  return 0;
}
