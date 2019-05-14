#include "scaffolding/score.h"
#include <stdlib.h>
void destroy_matrix_score(struct matrix_score *x)
{
	free(x->A);
	free(x);
}

int detect_anomal_diagonal(struct matrix_score *score, float threshold)
{
	int n_bucks = score->n_bucks;
	float sum_dia, sum_all, avg_dia, avg_all;
	sum_dia = 0;
	sum_all = 0;
	avg_dia = 0;
	avg_all = 0;
	int count_dia = 0, count_all = 0;
	for (int i = 0; i < n_bucks; i++) {
		for (int j = 0; j < n_bucks; j++) {
			if (i == j) {
				float sc = score->A[i*n_bucks+j];
				if (sc > -0.000001){
					sum_dia += sc;
					count_dia++;
				}
			} else {
				float sc = score->A[i*n_bucks+j];
				if (sc > -0.000001){
					sum_all += sc;
					count_all++;
				}
			}
		}
	}
	avg_all = sum_all / count_all;
	avg_dia = sum_dia / count_dia;
	return (sum_all + sum_dia > threshold && avg_dia > avg_all * 2);
}

