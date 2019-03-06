#include <stdlib.h>
#include <string.h>

#include "simple.h"

#include "io_utils.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

#define MIN_CON_LEN 7000
#define MIN_RATIO_COV 1.7
#define MIN_BRIDGE_LEG 7000
#define MIN_COMPONENT 250

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
  gint_t u;
  u = e[i].target;
  gint_t j = e[v[u].adj[0]].seq_len > MIN_BRIDGE_LEG ? v[u].adj[1] : v[u].adj[0];
  if (e[j].target == e[i].source && e[i].source == e[j].target)
    return 1;
  return 0;
}

/* \                         /
 *  \                       /  
 *   \                     / 
 *    \_______bridge______/
 *    /                   \
 *   /                     \
 *  /                       \
 * /                         \
 */                           
void find_bridge(struct asm_graph_t *g0)
{
  int i;
  int ret, flag, loop;
  gint_t r_src;
  gint_t src, dest;
  struct asm_edge_t *e = g0->edges;
  struct asm_node_t *v = g0->nodes;

  khash_t(khInt) *set_v; // for keeping visited nodes
  set_v = kh_init(khInt);

  gint_t *id_node, *id_edge, *cc_size; // for the connected component
  gint_t cc_id;
  id_node = malloc(g0->n_v * sizeof(gint_t));
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
    if (v[r_src].deg == 2 && v[dest].deg == 2){
      if (is_lg_leg(e, v, r_src, 2) && is_lg_leg(e, v, dest, 2)){ //neighbor edges must be large
          __VERBOSE("Edge %d\n - bridge\n", i);
          flag = 1;
      } else if (is_lg_leg(e, v, r_src, 1) && is_lg_leg(e, v, dest, 1)){ //1 neighbor edge must be large
          if (is_trivial_loop(e, v, i)){
            __VERBOSE("Edge %d\n - loop\n", i);
            flag = 1;
          }
      }
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

void find_forest(struct asm_graph_t *g0)
{
	init_clock();
    g_cov = get_genome_cov(g0);
    find_bridge(g0);

	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips\n");
	__VERBOSE_LOG("kmer_%d_graph_#1", "Number of nodes: %lld\n", g0->ksize,
							(long long)g0->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#1", "Number of edges: %lld\n", g0->ksize,
							(long long)g0->n_e);
    __VERBOSE("Genome walk coverage %d\n", (int)g_cov);
}

int main(int argc, char *argv[])
{
  struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
  char *path = argv[1];
  load_asm_graph(g0, path);
  find_forest(g0);
  return 0;
}
