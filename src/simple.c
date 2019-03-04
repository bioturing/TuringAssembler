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

void find_forest(struct asm_graph_t *g0)
{
	init_clock();
	__VERBOSE_LOG("INFO", "kmer size: %d\n", g0->ksize);
	__VERBOSE("\n+------------------------------------------------------------------------------+\n");
	__VERBOSE("Removing tips\n");
	__VERBOSE_LOG("kmer_%d_graph_#1", "Number of nodes: %lld\n", g0->ksize,
							(long long)g0->n_v);
	__VERBOSE_LOG("kmer_%d_graph_#1", "Number of edges: %lld\n", g0->ksize,
							(long long)g0->n_e);
}

int main(int argc, char *argv[])
{
  struct asm_graph_t *g0 = calloc(1, sizeof(struct asm_graph_t));
  char *path = argv[1];
  load_asm_graph(g0, path);
  find_forest(g0);
  return 0;
}
