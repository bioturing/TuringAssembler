#include <stdlib.h>
#include <string.h>

#include "simple.h"

#include "io_utils.h"
#include "k31hash.h"
#include "k63hash.h"
#include "k31_count.h"
#include "k63_count.h"
#include "utils.h"
#include "time_utils.h"
#include "verbose.h"

#define MIN_CON_LEN 7000

uint32_t get_genome_cov(struct asm_graph_t *g0)
{
  struct asm_edge_t *e = g0->edges;
  int i, cnt = 0;
  uint64_t cov = 0;

  for (i = 0 ; i < g0->n_e; i++){
    if (e[i].seq_len > MIN_CON_LEN){
      cov += e[i].count / (e[i].seq_len + 2 * g0->ksize);
      cnt++;
    }
  }
  assert(cnt > 0);
  return cov / cnt;
}

void find_forest(struct asm_graph_t *g0)
{
	init_clock();
    uint64_t g_cov = get_genome_cov(g0);
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
