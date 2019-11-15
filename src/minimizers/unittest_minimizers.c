//
// Created by BioTuring on 2019-11-11.
//

#include "minimizers.h"
#include "utils.h"
#include "khash.h"
#include "assembly_graph.h"

int mm_find_linear(struct mm_db_t *db, uint64_t mm)
{
	uint32_t i;
	for (i = 0; i < db->n; ++i) {
		if (db->mm[i] == mm) {
			return 1;
		}
	}
	return 0;
}

void test_bin_minimizers()
{
	char *data[3] = {"ACTAGTT", "AGTTCGT", "GGCAATC"};
	char ref[10] = "ACAGTTCCGT";
	struct mm_db_t *mm_db[4];
	uint32_t *s;
	int i;

	for (i = 0; i < 3; ++i) {
		s = seq2uint32t(data[i], 7);
		mm_db[i] = mm_index_bin_str(s, 3, 3, 7);
	}
	s = seq2uint32t(ref, 10);
	mm_db[3] = mm_index_bin_str(s, 3, 3, 10);

	for (i = 0; i < mm_db[3]->n; ++i) {
		printf("minimiers %d %lu, pos_ref: %d,  %d %d %d\n", i,
			mm_db[3]->mm[i] , mm_db[3]->p[i],
			mm_find_linear(mm_db[0], mm_db[3]->mm[i]),
			mm_find_linear(mm_db[1], mm_db[3]->mm[i]),
			mm_find_linear(mm_db[2], mm_db[3]->mm[i]));
	}
	printf("EXPECT\n");
	printf("minimiers 0 5188146770730811392, pos_ref: 1,  0 0 0\n");
	printf("minimiers 1 3170534137668829184, pos_ref: 2,  1 0 0\n");
	printf("minimiers 2 17582052945254416384, pos_ref: 4,  0 1 0\n");
	printf("minimiers 3 15276209936040722432, pos_ref: 5,  0 0 0\n");
	printf("minimiers 4 6341068275337658368, pos_ref: 6,  0 0 0\n");
}

void test_str_minimizers()
{
	char *data[3] = {"ACTAGTT", "AGTTCGT", "GGCAATC"};
	char ref[10] = "ACAGTTCCGT";
	struct mm_db_t *mm_db[4];
	uint32_t *s;
	int i;

	for (i = 0; i < 3; ++i) {
		mm_db[i] = mm_index_char_str(data[i], 3, 3, 7);
	}
	mm_db[3] = mm_index_char_str(ref, 3, 3, 10);

	for (i = 0; i < mm_db[3]->n; ++i) {
		printf("minimiers %d %lu, pos_ref: %d,  %d %d %d\n", i,
		       mm_db[3]->mm[i] , mm_db[3]->p[i],
		       mm_find_linear(mm_db[0], mm_db[3]->mm[i]),
		       mm_find_linear(mm_db[1], mm_db[3]->mm[i]),
		       mm_find_linear(mm_db[2], mm_db[3]->mm[i]));
	}
	printf("EXPECT\n");
	printf("minimiers 0 5188146770730811392, pos_ref: 1,  0 0 0\n");
	printf("minimiers 1 3170534137668829184, pos_ref: 2,  1 0 0\n");
	printf("minimiers 2 17582052945254416384, pos_ref: 4,  0 1 0\n");
	printf("minimiers 3 15276209936040722432, pos_ref: 5,  0 0 0\n");
	printf("minimiers 4 6341068275337658368, pos_ref: 6,  0 0 0\n");
}

void test_index_mm_graph()
{
	char *data[3] = {"ACTAGTTGTGTGCA", "AGTAGCGTAGTTCA", "GGCAATAGTTTCA"};
	int i;
	struct asm_graph_t g;
	g.n_e = 3;
	g.edges = calloc(3, sizeof(struct asm_edge_t));

	for (i = 0; i < g.n_e; ++i) {
		g.edges[i].seq = seq2uint32t(data[i], 14);
		g.edges[i].seq_len = 14;
	}
	struct mm_db_edge_t *db = mm_index_edges(&g, 3, 4);
	printf("EXPECT\n");
	printf("Redundant minimizers AGT: contig: 0, pos: 3\n");
	printf("Singleton TTG: contig: 0, pos: 5\n");
	printf("Singleton GTA: contig: 1, pos: 1\n");
	printf("Singleton TTT: contig: 2, pos: 8\n");
	printf("Singleton GCG: contig: 1, pos: 4\n");
	printf("Singleton GTG: contig: 0, pos: 7\n");
	printf("Singleton TGC: contig: 0, pos: 10\n");
	printf("Singleton GCA: contig: 2, pos: 1\n");
	printf("Singleton AGC: contig: 1, pos: 3\n");
}

