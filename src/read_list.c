#include "read_list.h"

void read_list_init(struct read_list_t *rlist)
{
	rlist->n = 0;
	rlist->m = 0;
	rlist->reads = NULL;
}

void load_reads_from_fastq(char *r_path, struct read_list_t *rlist)
{
	__VERBOSE("Loading reads\n");
	read_list_init(rlist);
	FILE *f = fopen(r_path, "r");

	while (!feof(f)){
		char *name = (char *) calloc(MAX_READ_BUF, sizeof(char));
		fscanf(f, "@%[^\n]\n", name);
		name = (char *) realloc(name, (strlen(name) + 1) * sizeof(char));

		char *seq = (char *) calloc(MAX_READ_BUF, sizeof(char));
		fscanf(f, "%s\n", seq);
		seq = (char *) realloc(seq, (strlen(seq) + 1) * sizeof(char));

		char *info = (char *) calloc(MAX_READ_BUF, sizeof(char));
		fscanf(f, "+%s\n", info);
		info = (char *) realloc(info, (strlen(info) + 1) * sizeof(char));

		struct read_t r = {.name = name,
				.seq = seq,
				.info = info,
				.len = strlen(seq)};
		push_read_to_list(r, rlist);
	}
	fclose(f);
}

void push_read_to_list(struct read_t read, struct read_list_t *rlist)
{
	if (rlist->n == rlist->m){
		rlist->m = rlist->m == 0? 1 : (rlist->m << 1);
		rlist->reads = realloc(rlist->reads,
				rlist->m * sizeof(struct read_t));
	}
	rlist->reads[rlist->n] = read;
	++(rlist->n);
}

void load_read_pairs_from_fastq(char *r1_path, char *r2_path,
		struct read_list_t *list1, struct read_list_t *list2)
{
	load_reads_from_fastq(r1_path, list1);
	load_reads_from_fastq(r2_path, list2);
}

void read_list_destroy(struct read_list_t *rlst)
{
	for (int i = 0; i < rlst->n; ++i){
		free(rlst->reads[i].name);
		free(rlst->reads[i].seq);
		free(rlst->reads[i].info);
	}
	free(rlst->reads);
}
