#ifndef __SCAFFOLDING_SCORE_H__
#define __SCAFFOLDING_SCORE_H__
struct bucks_score {
	float score;
};

struct matrix_score{
	int n_bucks;
	float *A;
};

void destroy_matrix_score(struct matrix_score *x);
int detect_anomal_diagonal(struct matrix_score *score, float threshold);
#endif /* __SCAFFOLDING_SCORE_H__ */
