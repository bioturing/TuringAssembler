//
// Created by che on 02/12/2019.
//

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "yeast_analyze_utils.h"

void load_ground_truth_file(int ((**a)[2]), int *n) {
	FILE *in = fopen("/home/che/bioturing/data/yeast/metadata/neo4j/truth_neo4j.txt", "r");
	char *s = calloc(20, 1);
	fscanf(in, "%s\n", s);
	int x, y;
	int (*res)[2] = NULL, n_res = 0;
	while (fscanf(in, "%d,%d\n", &x, &y) != EOF) {
		res = realloc(res, (n_res + 1) * sizeof(int[2]));
		res[n_res][0] = x;
		res[n_res][1] = y;
		n_res++;
	}
	*a = res;
	*n = n_res;
}

void load_share_barcode_file(int (**share_bx_list)[3], int *n_share)
{
	FILE *in = fopen("/home/che/bioturing/data/yeast/metadata/all_barcodes_count.txt", "r");
	int x, y, z;
	int (*res_share)[3] = NULL, n_res = 0;
	while (fscanf(in, "%d %d %d\n", &x, &y, &z) != EOF) {
		res_share = realloc(res_share, (n_res+1) * sizeof(int[3]));
		res_share[n_res][0] = x;
		res_share[n_res][1] = y;
		res_share[n_res][2] = z;
		n_res++;
	}
	*share_bx_list = res_share;
	*n_share = n_res;
}

int lessa2(int *a, int *b)
{
	return (a[0] < b[0]) || (a[0] == b[0] && a[1] < b[1]);
}

int get_share_barcode_fromfile(int a[2], int (*share_bx_list)[3], int n_share)
{
	int l = 0, r = n_share, mid;
	while (l < r-1) {
		mid = (l+r+1)/2;
		if (lessa2(a, share_bx_list[mid])) {
			r = mid;
		} else {
			l = mid;
		}
	}
	if (share_bx_list[l][0] == a[0] && share_bx_list[l][1] == a[1]) {
		return share_bx_list[l][2];
	}
	return -1;
//	assert(0 && "not found");
}

void write_pair_info(FILE *out, struct info_pair *t)
{
	fprintf(out, "x:%d y:%d {\n", t->x, t->y);
	fprintf(out, "\tlenx: %d\n", t->len_x);
	fprintf(out, "\tleny: %d\n", t->len_y);
	fprintf(out, "\tn_shortest: %d\n", t->n_shortest);
	fprintf(out, "\tlen_shortest: %d\n", t->len_shortest);
	fprintf(out, "\tshare_bx: %d\n", t->share_bx);
	fprintf(out, "}\n");
}
