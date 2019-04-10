#include <stdlib.h>
#include "assembly_graph.h"
#include "k31hash.h"
#include "pthread.h"
#include "verbose.h"
#include <string.h>

struct bucks_score {
	float score;
	float bin_distance;
};

uint32_t min(uint32_t a, uint32_t b)
{
	if (a < b) return a;
	else return b;
}

uint32_t roundint(float x)
{ 
	return (x)>=0?(int)((x)+0.5):(int)((x)-0.5);
}

float getScoreBucks(struct barcode_hash_t *buck0,struct barcode_hash_t *buck1) 
{
	const uint32_t thres_cnt = 30;
	uint32_t cnt0 = 0, cnt1 = 0, res2 = 0;

	for (uint32_t i = 0; i < buck1->size; ++i) {
		if (buck1->cnts[i] != (uint32_t)(-1) && buck1->cnts[i] >= thres_cnt) {
			cnt1++;
		}
	}

	for (uint32_t i = 0; i < buck0->size; ++i) {
		if ((buck0->keys[i]) != (uint64_t)(-1) && buck0->cnts[i] >= thres_cnt) {
			cnt0++;
			uint32_t tmp = barcode_hash_get(buck1, buck0->keys[i]);
			if (tmp != BARCODE_HASH_END(buck1) && buck1->cnts[tmp] >= thres_cnt) {
				res2++;
			}
		}
	}
	return 1.0 * res2 / (cnt0 + cnt1);
}

uint32_t abssub(uint32_t a, uint32_t b) {
	if (a>b) 
		return a-b;
	else 
		return b-a;
}

struct bucks_score getScore(uint32_t i0, uint32_t i1, struct asm_graph_t *g) {
	uint32_t n_bucks = 5;
	struct asm_edge_t *e0 = &g->edges[i0], *e1 = &g->edges[i1];
	uint32_t n0_bucks = (get_edge_len(e0) + g->bin_size-1) / g->bin_size;
	uint32_t n1_bucks = (get_edge_len(e1) + g->bin_size-1) / g->bin_size;
	float score2=0;
	struct bucks_score score;

	if (n0_bucks < n_bucks || n1_bucks < n_bucks) {
		score.bin_distance = -1;
		return score;
	}
	float res = 0;
	float maxtmp =0 ;
	float sum = 0;

	uint32_t rev_i0 = g->edges[i0].rc_id;
	struct asm_edge_t *rev_e0 = &g->edges[rev_i0];
	for (uint32_t i = 0; i < n_bucks; ++i) {
		for (uint32_t j = 0; j < n_bucks; ++j) {
			float tmp = 0;
			if (i!=0 || j!=0) {
				tmp = getScoreBucks(&rev_e0->bucks[i], &e1->bucks[j]) / (n_bucks * n_bucks-1);
			} else{
			}
			sum += tmp;
			score2 += tmp*(i+j); 
			res += tmp;
		}
//		__VERBOSE("\n");
	}
//	__VERBOSE("\n");
	score.score = res;
	score.bin_distance = score2/res;
	return score;
}

void listContig(struct asm_graph_t *g, FILE *out_file) {
	const uint32_t thres_len_e = 20000; 
	float thres_score = 0.025;
	uint32_t *listE = NULL;
	uint32_t n_e=0;
	for (uint32_t e = 0; e < g->n_e; ++e) {
		uint32_t len = get_edge_len(&g->edges[e]);
		if (len > thres_len_e) {
			++n_e;
			listE = realloc(listE, n_e*sizeof(uint32_t));
			listE[n_e-1] = e; 
		}
	}
	__VERBOSE("n_e: %d\n", n_e);
	for (uint32_t i = 0; i < n_e; i++) {
		__VERBOSE("%d\n",i);
		uint32_t e0 = listE[i];
		for (uint32_t i1 = 0; i1 < n_e; i1++) {
			uint32_t e1 = listE[i1];
			struct bucks_score score = getScore(e0, e1, g);
			uint32_t check = 0;
			if (e0 == e1){
				float cvr = get_genome_coverage(g);
				if (__get_edge_cov(&g->edges[e0], g->ksize)/cvr > 1.8) {
					check = 1;
				}
			} else if (e0 == g->edges[e1].rc_id) {
			} else {
				if (score.bin_distance > 0 ) {
					check = 1;
				}
			}
			if (check) {
				if (score.score > thres_score) {
					fprintf(out_file, "score: %f center:%f edge: %d %d\n", score.score, score.bin_distance, e0, e1); 
				}
			}
		}
	}
	fclose(out_file);
}

struct contig_edge {
	uint32_t src, des;
	uint32_t rv_src, rv_des;
	float score0, score1;
};

void swap(uint32_t *a, uint32_t *b)
{ 
	uint32_t tmp;
	tmp = *a; 
	*a = *b; 
	*b = tmp;
}

void normalize_min_index(struct asm_graph_t *g, struct contig_edge *e)
{
	uint32_t rc_id_src = g->edges[e->src].rc_id;
	if (rc_id_src < e->src) {
		e->src = rc_id_src;
		e->rv_src ^= 1;
	}
	uint32_t rc_id_des = g->edges[e->des].rc_id;
	if (rc_id_des < e->des) {
		e->des = rc_id_des ;
		e->rv_des ^= 1;
	}
	if (e->src > e->des) {
		swap(&e->src, &e->des);
		swap(&e->rv_src, &e->rv_des);
		e->rv_src ^= 1 ;
		e->rv_des ^= 1;
	}
}

void normalize_one_dir(struct asm_graph_t *g, struct contig_edge *e)
{
	if (e->rv_src == 1) {
		e->src = g->edges[e->src].rc_id;
		e->rv_src = 0;
	}
	if (e->rv_des == 1) {
		e->des = g->edges[e->des].rc_id;
		e->rv_des = 0;
	}
}

void add_contig_edge(struct asm_graph_t *g,struct contig_edge *listE, uint32_t pos, uint32_t src, uint32_t des, float score0, float score1)
{
	struct contig_edge e;
	e.src = src;
	e.des = des;
	e.score0 = score0;
	e.score1 = score1;
	e.rv_src = 0;
	e.rv_des = 0;
	normalize_min_index(g, &e);
	listE[pos] = e;
}

int less_contig_edge(const void *e0,const void *e1)
{
	struct contig_edge * v0 = (struct contig_edge *) e0;
	struct contig_edge * v1 = (struct contig_edge *) e1;
	return (v0->src > v1->src || (v0->src == v1->src && v0->des > v1->des));
}

int less_uint32(const void *e0,const void *e1)
{
	return *(uint32_t*)(e0) > *(uint32_t*)(e1);
}

uint32_t equal_contig_edge(struct contig_edge *e0,struct contig_edge *e1)
{
	return (e0->src == e1->src && e0->des == e1->des);
}

uint32_t better_contig_edge(struct contig_edge *e0, struct contig_edge *e1)
{
	return (e0->score0 > e1->score0);
}

void unique_edge(struct contig_edge *listE, uint32_t *n_e)
{
	uint32_t new_n_e = 0;
	for (uint32_t i = 0; i < *n_e; ) {
		uint32_t j = i;
		struct contig_edge best;
		best.score0 = 0;
		while (equal_contig_edge(&listE[j], &listE[i])){
			if (better_contig_edge(&listE[j], &best)){
				best = listE[j];
			}
			++j;
		}
		listE[new_n_e++] = best;
		i = j;
	}
	*n_e = new_n_e;
}

uint32_t *binary_search(uint32_t *list, uint32_t n, uint32_t value)
{
	uint32_t l = 0, r = n - 1;
	while (l != r) {
		uint32_t m = (l + r) / 2;
		if (value > list[m]) 
			l = m + 1;
		else
			r = m;
	}
	if (list[l] != value)
		return list-1;
	return &list[l];
}

void unique_vertex(uint32_t *listV, uint32_t *n_v)
{
	uint32_t new_n_v = 0;
	for (uint32_t i = 0; i < *n_v; ) {
		uint32_t j = i;
		while (listV[i] == listV[j]) j++;
		listV[new_n_v++] = listV[i];
		i = j;
	}
	*n_v = new_n_v;
}

// **listV
void build_V_from_E(struct contig_edge *listE, uint32_t n_e, uint32_t **listV, uint32_t *n_v)
{
	*listV = NULL; 
	*n_v = 0;
	for (uint32_t i = 0; i < n_e; i++) {
		__VERBOSE("%d %d BBBB\n", listE[i].src, listE[i].des);

		++(*n_v);
		uint32_t t = (*n_v)*sizeof(uint32_t);
		*listV = realloc(*listV, t);
		(*listV)[(*n_v)-1] = listE[i].src;

		++(*n_v);
		*listV = realloc(*listV, (*n_v)*(sizeof(uint32_t)));
		(*listV)[(*n_v)-1] = listE[i].des;
	}
	for (uint32_t i = 0; i < *n_v; i++) {
		__VERBOSE("%d ", (*listV)[i]);
	}
	__VERBOSE("\n");
	qsort(*listV, *n_v, sizeof(uint32_t), less_uint32);
	for (uint32_t i = 0; i < *n_v; i++) {
		__VERBOSE("%d ", (*listV)[i]);
	}
	__VERBOSE("\n");
	qsort(*listV, *n_v, sizeof(uint32_t), less_uint32);
	unique_vertex(*listV, n_v);
//	__VERBOSE("listV %d\n", listV);
}

void dfs_eulerian(uint32_t x, uint32_t *E, int *head, uint32_t *next, uint32_t n_v, uint32_t *res, uint32_t *n_res)
{
	// not verify coressness yet
	uint32_t y = E[head[x]];
	while (head[x] != -1) {
		head[x] = next[head[x]];
		dfs_eulerian(y, E, head, next, n_v, res, n_res);
	}
	res[*n_res++] = x;
}

void algo_find_eulerian(uint32_t *E, uint32_t n_e, uint32_t *head, uint32_t *next, uint32_t n_v, uint32_t *res, uint32_t *n_res)
{
	int *tmp_head = calloc(n_v, sizeof(uint32_t));
	for (uint32_t i = 0; i < n_v; ++i) {
		tmp_head[i] = head[i];
	}
	uint32_t *deg = calloc(n_v, sizeof(uint32_t));
	for (uint32_t i = 0; i < n_v; i++) {
		deg[i] = 0;
	}
	for (uint32_t i = 0; i < n_e; i++) {
		deg[E[i]]++;
	}
	uint32_t count = 0, start = 0;
	for (uint32_t i = 0; i < n_v; i++) if (deg[i]&1)
		++count;
	if (count & 1) {
		__VERBOSE("BUG");
	} else if (count == 0) {
		start = 0;
	} else if (count == 2) {
		for (uint32_t i = 0; i < n_v; i++) if (deg[i]&1) {
			start = i;
			break;
		}
	}
	dfs_eulerian(start, E, tmp_head, next, n_v, res, n_res);
	for (uint32_t i = 0; i < *n_res; i++) {
		__VERBOSE("%d ",res[i]);
	}
}

void dfs_hamiltonian(int x, uint32_t depth, uint32_t *E, int *head, uint32_t *next, uint32_t *mark, uint32_t n_v, 
		uint32_t *res, uint32_t *best_res, uint32_t *best_n_res)
{
	mark[x]--;
	res[depth-1] = x;
	if (depth > *best_n_res) {
		__VERBOSE("Find a path");
		*best_n_res = depth;
		for (uint32_t i = 0; i < depth; i++){
			best_res[i] = res[i];
		}
	} 
	for (int i = head[x]; i != -1; i = next[i]) if (mark[E[i]] != 0) {
		dfs_hamiltonian(E[i], depth + 1, E, head, next, mark, n_v, res, best_res, best_n_res);
	}
	mark[x]++;
}

void print_seq(FILE *fp, uint32_t index, char *seq, uint32_t len, uint32_t cov)
{
	fprintf(fp, ">SEQ_%lld_length_%lld_count_%llu\n", (long long)index,
		(long long)len, (long long unsigned)cov);
	gint_t k = 0;
	char *buf = alloca(81);
//	printf("%s\n", seq);
	while (k < len) {
		gint_t l = __min(80, len - k);
		__VERBOSE("%ld %ld %d\n",l, k, len);
		memcpy(buf, seq + k, l);
		buf[l] = '\0';
		fprintf(fp, "%s\n", buf);
		k += l;
	}
	while (k < len) {
		gint_t l = __min(80, len - k);
		__VERBOSE("%ld %ld %d\n",l, k, len);
		memcpy(buf, seq + k, l);
		buf[l] = '\0';
		fprintf(fp, "%s\n", buf);
		k += l;
	}
}

void algo_find_hamiltonian(FILE *out_file, struct asm_graph_t *g, uint32_t *E, uint32_t n_e, int *head, uint32_t *next, uint32_t n_v, uint32_t *res, uint32_t *n_res, uint32_t *listV)
{
	uint32_t *mark = calloc(n_v, sizeof(uint32_t));
	for (uint32_t i = 0; i < n_v; i++) {
		float cvr = get_genome_coverage(g);
		float cov_times = (__get_edge_cov(&g->edges[listV[i]], g->ksize)/cvr) ;
		__VERBOSE ("%f xxx\n", cov_times);
		mark[i] = roundint(cov_times);
	}
	uint32_t *best_res = calloc(n_v, sizeof(uint32_t)), best_n_res = 0 ;
	for (uint32_t i = 0; i < n_v; i++){
		__VERBOSE("dfs from %d\n", i);
		dfs_hamiltonian(i,  1, E, head, next, mark, n_v, res, best_res, &best_n_res);
	}
	__VERBOSE("best length %d\n", best_n_res);
	char *seq = NULL, *total_seq = NULL, *NNN = NULL;
	NNN = calloc(1000, sizeof(char));
	for (uint32_t i = 0; i < 1000; i++) 
		NNN[i] = 'N';
	uint32_t seq_len = 0, total_len = 0;
	for(uint32_t i = 0; i < best_n_res; i++) {
		uint32_t e = listV[best_res[i]];
		uint32_t len = dump_edge_seq(&seq, &seq_len, &g->edges[e]);
		total_seq = realloc(total_seq, (total_len + len) * sizeof(char));
		memcpy(total_seq+total_len, seq, len);
		total_len += len;
	}
	print_seq(out_file, 0, total_seq, total_len, 1);
	uint32_t count = 0;
//	for (uint32_t e = 0; e < g->n_e; ++e) {
//		uint32_t len = get_edge_len(&g->edges[e]);
//		if (len < 20000) {
//			count++;
//			seq_len = 0;
//			uint32_t len = dump_edge_seq(&seq, &seq_len, &g->edges[e]);
//			print_seq(out_file, count, seq, seq_len, 1); 
//		}
//	}
	for(uint32_t i = 0; i < best_n_res; i++) {
		uint32_t e = listV[best_res[i]];
		__VERBOSE("%d ", e);
	}
	fclose(out_file);
}

void find_hamiltonian_contig_edge(FILE *out_file, struct asm_graph_t *g, struct contig_edge *listE_ori, uint32_t n_e){
	struct contig_edge *list_one_dir_E = calloc(2*n_e, sizeof(struct contig_edge));
	for (uint32_t i = 0; i < n_e; i++) {
		list_one_dir_E[i] = listE_ori[i];
		normalize_one_dir(g, list_one_dir_E+i);
	}
	for (uint32_t i = 0; i < n_e; i++) {
		list_one_dir_E[i + n_e] = listE_ori[i];
		list_one_dir_E[i + n_e].src = g->edges[list_one_dir_E[i].des].rc_id;
		list_one_dir_E[i + n_e].des = g->edges[list_one_dir_E[i].src].rc_id;
	}
	n_e *=2;
	uint32_t *listV = NULL , n_v; 
	build_V_from_E(list_one_dir_E, n_e, &listV, &n_v);
	__VERBOSE("n_v: %d \n", n_v);
	for (uint32_t i = 0; i < n_v; i++) {
		__VERBOSE("%d ", listV[i]);
	}

	uint32_t m = 0;
	uint32_t *E = calloc(n_e * 2, sizeof(uint32_t)), *next = calloc(n_e * 2, sizeof(uint32_t));
	int *head = calloc(n_v, sizeof(uint32_t));
	for (uint32_t i = 0; i < n_v; i++) {
		head[i] = -1;
	}
	for (uint32_t i = 0; i < n_e; i++) {
		int u = list_one_dir_E[i].src;
		u = binary_search(listV, n_v, u) - listV;
		int v = list_one_dir_E[i].des;
		v = binary_search(listV, n_v, v) - listV;
		__VERBOSE("edge: %d %d \n", u, v);
		if (u == -1 || v == -1){
			__VERBOSE("%d %d ERRRRRR\n", list_one_dir_E[i].src, list_one_dir_E[i].des);
		}
//		E[m] = u;
//		next[m] = head[v];
//		head[v] = m;
//		m++;
		E[m] = v;
		next[m] = head[u];
		head[u] = m;
		m++;
	}
	__VERBOSE("\nE:\n"); 
	__VERBOSE("\n"); 
	uint32_t *res = calloc(n_v, sizeof(uint32_t)), n_res=0;
	algo_find_hamiltonian(out_file, g,E, m, head, next, n_v, res, &n_res, listV);
	// todo free all pointer
}

void print_gfa_from_E(struct asm_graph_t *g, struct contig_edge *listE, uint32_t n_e)
{
	struct contig_edge *list_one_dir_E = calloc(2*n_e, sizeof(struct contig_edge));
	for (uint32_t i = 0; i < n_e; i++) {
		list_one_dir_E[i] = listE[i];
		normalize_min_index(g, list_one_dir_E+i);
	}
	uint32_t *listV = NULL , n_v; 
	build_V_from_E(list_one_dir_E, n_e, &listV, &n_v);
	FILE *fp = fopen("gr" , "w");
	for (uint32_t i = 0; i < n_v; i++) {
		fprintf(fp,"S\t%d\tAAA\tKC:i:1\n", listV[i]);
//		__VERBOSE("S\t%d\tAAA\tKC:i:1\n", g->edges[listV[i]].rc_id);
	}

	for (uint32_t i = 0; i < n_e; i++) {
		fprintf(fp, "L\t%d\t%c\t%d\t%c\t45M\n", 
			list_one_dir_E[i].src, list_one_dir_E[i].rv_src == 0?'+':'-', list_one_dir_E[i].des, list_one_dir_E[i].rv_des == 0?'+':'-');
		fprintf(fp, "L\t%d\t%c\t%d\t%c\t45M\n", 
			list_one_dir_E[i].des, list_one_dir_E[i].rv_des == 0?'-':'+', list_one_dir_E[i].src, list_one_dir_E[i].rv_src == 0?'-':'+');
	}
}

void build_graph_2(FILE *fp, FILE *out_file, struct asm_graph_t *g)
{
	for (uint32_t i = 0 ; i < 5; ++i) {
		uint32_t t = fscanf(fp, "%*[^\n]\n");
	}
	float score, center;
	uint32_t src, des;
	struct contig_edge *listE = NULL;
	uint32_t n_e=0;
	while (fscanf(fp,"score: %f center:%f edge: %d %d\n", &score, &center, &src, &des) !=EOF){
//		__VERBOSE("%f", score);
		n_e += 1;
		listE = realloc(listE, n_e*sizeof(struct contig_edge));
		add_contig_edge(g, listE, n_e-1, src, des, score, center);
	}
//	__VERBOSE("n_e: %d\n", n_e);
	qsort(listE, n_e, sizeof(struct contig_edge), less_contig_edge);
	unique_edge(listE, &n_e);

	// finding hamiltonian path
	find_hamiltonian_contig_edge(out_file, g, listE, n_e);

//	print_gfa_from_E(g, listE, n_e);
	// print as gfa graph
//	__VERBOSE("n_v: %d \n",n_v);
//	__VERBOSE("new_n_v: %d \n",n_v);
}

