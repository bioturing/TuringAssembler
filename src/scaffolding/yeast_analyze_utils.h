//
// Created by che on 02/12/2019.
//

#ifndef SKIPPING_YEAST_ANALYZE_UTILS_H
#define SKIPPING_YEAST_ANALYZE_UTILS_H

struct info_pair{
    int x,y;
    int len_x, len_y;
    int n_shortest, len_shortest;
    int share_bx;
};

void load_ground_truth_file(int ((**a)[2]), int *n);
void write_pair_info(FILE *out, struct info_pair *t);
int get_share_barcode_fromfile(int a[2], int (*share_bx_list)[3], int n_share);
void load_share_barcode_file(int (**share_bx_list)[3], int *n_share);

#endif //SKIPPING_YEAST_ANALYZE_UTILS_H
