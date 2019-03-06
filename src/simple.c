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
#define MIN_BRIDGE_LEG 7000
#define MIN_COMPONENT 250
#define MAX_COMPONENT_REGION 3000
#define MIN_LAYER 10

static uint64_t g_cov;

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
  gint_t next_e, u, dest, lg;

  *k = NULL;

  if (e[e_i].seq_len < MIN_BRIDGE_LEG)
    return 0;

  aqueue_add(q, e[e_i].target); //add the big edge to queue
  while (1){
    n_items = q->n;
    if (n_items == 0)
      break;
    for (j = 0; j < n_items; j++){
      u = aqueue_pop(q);
      for (i = 0 ; i < v[u].deg; i++){
        next_e = v[u].adj[i];
        if (kh_get(khInt, set_v, e[next_e].target) == kh_end(set_v)){ // visited or not
          if (e[next_e].seq_len > MIN_BRIDGE_LEG){
             n_larges++;
             lg = next_e;
          }
          aqueue_add(q, e[next_e].target); //add neighbor edge
          kh_put(khInt, set_v, e[next_e].target, &missing);
          comp_sz += e[next_e].seq_len; //add size of the edge to total size
        }
      }
    }
    if (++cnt > MIN_LAYER)
      break;
  }

  printf("Visisted nodes: %d, n_larges: %d, edge %d\n", kh_size(set_v), n_larges, e_i);
  for (i = q->p; i < q->n; i++){
    u = q->e[i];
    for (j = 0; j < v[u].deg; j++){
      if (kh_get(khInt, set_v, e[v[u].adj[j]].target) == kh_end(set_v)){
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

  int i;
  int ret, flag, loop;
  gint_t r_src;
  gint_t src, dest;
  gint_t lg;
  struct asm_edge_t *e = g0->edges;
  struct asm_node_t *v = g0->nodes;

  khash_t(khInt) *set_v; // for keeping visited nodes
  set_v = kh_init(khInt);

  gint_t *id_node, *id_edge, *cc_size; // for the connected component
  gint_t cc_id;
  id_edge = malloc(g0->n_e * sizeof(gint_t));
  cc_size = NULL;
  asm_edge_cc(g0, id_edge, &cc_size);

  for (i = 0 ; i < g0->n_e; i++){
    cc_id = id_edge[i];
    if (cc_size[cc_id] < MIN_COMPONENT || kh_get(khInt, set_v, i) != kh_end(set_v)) //must came from large component and not be found
      continue;
    src = e[i].source;
    dest = e[i].target;
    r_src = v[src].rc_id;
    flag = 0;
    if (is_bridge(e, v, i)){
      __VERBOSE("Edge %d\n - bridge\n", i);
      flag = 1;
    }
    if (is_trivial_loop(e, v, i)){
      __VERBOSE("Edge %d\n - loop\n", i);
      flag = 1;
    }
    if (is_simple_tandem(e, v, i, &lg)){
      __VERBOSE("Edge %d\n - simple complex\n", i);
      flag = 1;
    }
    if (flag){
      kh_put(khInt, set_v, i, &ret);
      kh_put(khInt, set_v, e[i].rc_id, &ret);
    }
  }
  kh_destroy(khInt, set_v);
  free(id_node);
  free(id_edge);

}

int main(int argc, char *argv[])
{
  struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
  char *path = argv[1];
  load_asm_graph(g0, path);
  gint_t k;
  is_simple_tandem(g0->edges, g0->nodes, 661461, &k);
  find_forest(g0);
  return 0;
}
