//
// Created by che on 04/12/2019.
//

#include <stdlib.h>
#include "minimizers/count_barcodes.h"
#include "cluster_molecules.h"
#include "log.h"
#include "utils.h"
#include "io_utils.h"

#define MIN_SHARE_BARCODE_COUNT 100
#define MIN_READ_PAIR_COUNT 1
#define VERY_SHORT_EDGE_LEN 250
#define LONG_PATH 10
#define SHORT_PATH 2
#define MIN_PAIR_SUPPORT_PAIR_END 1
#define MIN_PAIR_SUPPORT_PAIR_END_SOFT 0
#define MIN_SHARED_BARCODE_RATIO 0.005
#define MOLECULE_DENSITY 5000

struct barcode_graph {
    int *first_index, *next_index, *edges, *is_del;
    int n_nodes, n_edges;
};

struct barcode_graph *init_for_barcode_graph(int n)
{
	struct barcode_graph *bg = calloc(1, sizeof(struct barcode_graph));
	bg->n_nodes = n;
	bg->first_index = calloc(n, sizeof(int));
	for (int i = 0; i < bg->n_nodes; i++)
		bg->first_index[i] = -1;
	return bg;
}

void destroy_barcode_graph(struct barcode_graph *bg)
{
	free(bg->first_index);
	free(bg->next_index);
	free(bg->edges);
	free(bg->is_del);
	free(bg);
}

void append_list_edge(struct barcode_graph *bg, int n_pairs, int *list_edges)
{
	bg->edges = list_edges;
	bg->next_index = calloc(n_pairs * 2, sizeof(int));
	bg->n_edges = n_pairs;
	bg->is_del = calloc(n_pairs * 2, sizeof(int));
	//todo bg->edge just only half of list edge because this graph is one direction
	for (int i = 0; i < n_pairs; i++) {
		int u = list_edges[i << 1];
		int v = list_edges[(i << 1) + 1];
		bg->next_index[i << 1] = bg->first_index[v];
		bg->first_index[v] = i << 1;
		bg->next_index[(i << 1) + 1] = bg->first_index[u];
		bg->first_index[u] = (i << 1) + 1;
	}
}

int have_edge(struct barcode_graph *bg, int u, int v)
{
	for (int i = bg->first_index[u]; i != -1; i = bg->next_index[i]) {
		if (bg->edges[i] == v)
			return 1;
	}
	return 0;
}

void del_edge(struct barcode_graph *bg, int index)
{
	log_trace("%d %d", bg->edges[index & (-2)], bg->edges[index | 1]);
	bg->is_del[index] = 1;
	bg->is_del[index ^ 1] = 1;
}

void del_pair_edge(struct barcode_graph *bg, int index)
{
	del_edge(bg, index);
	del_edge(bg, index ^ 2);
}

void del_four_edge(struct barcode_graph *bg, int index)
{
	del_pair_edge(bg, index);
	del_pair_edge(bg, index ^ 4);
}

void del_eight_edge(struct barcode_graph *bg, int index)
{
	del_four_edge(bg, index);
	del_four_edge(bg, index ^ 8);
}

void filter_bulge(struct barcode_graph *bg)
{
	int *near = calloc(1000, 4);
	for (int node = 0; node < bg->n_nodes; node++) {
		int deg_out = 0;
		for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
			if (bg->is_del[j] == 0) {
				if (j & 1) {
					near[deg_out] = bg->edges[j];
					deg_out++;
				}
			}
		if (deg_out == 2) {
			if (have_edge(bg, near[0], near[1])) {
				for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
					if (bg->edges[j] == near[1]) {
						del_pair_edge(bg, j);
						log_debug("Del bulge %d %d", node, near[1]);
					}
			} else if (have_edge(bg, near[1], near[0]))
				for (int j = bg->first_index[node];
				     j != -1; j = bg->next_index[j])
					if (bg->edges[j] == near[0]) {
						del_pair_edge(bg, j);
						log_debug("Del bulge %d %d", node, near[0]);
					}
		}
	}
	free(near);
}

void filter_by_deg(struct barcode_graph *bg, int thres_deg)
{
	for (int node = 0; node < bg->n_nodes; node++) {
		int deg_out = 0;
		for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
			if (bg->is_del[j] == 0) {
				if (j & 1)
					deg_out++;
			}
		if (deg_out > thres_deg) {
			for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
				if (bg->is_del[j] == 0 && j & 1) {
					log_debug("Del by deg %d %d", node, bg->edges[j]);
					del_pair_edge(bg, j);
				}
		}
	}

	for (int node = 0; node < bg->n_nodes; node++) {
		int deg_in = 0;
		for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
			if (bg->is_del[j] == 0) {
				if ((j & 1) == 0)
					deg_in++;
			}
		if (deg_in > thres_deg) {
			for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
				if (bg->is_del[j] == 0 && (j & 1) == 0) {
					log_debug("Del by deg %d %d", node, bg->edges[j]);
					del_pair_edge(bg, j);
				}
		}
	}
}

void filter_complex_barcode_graph(struct barcode_graph *bg)
{
//	filter_bulge(bg);
	filter_by_deg(bg, 4);
//	filter_bulge(bg);
//	filter_by_deg(bg, 1);
}

uint64_t get_share_read_pair(struct mini_hash_t *rp_table, int u, int v)
{
	//todo @huu
	uint64_t t = GET_CODE(u, v);
	uint64_t *val = mini_get(rp_table, t);
	if (val == (uint64_t *) EMPTY_BX)
		return 0;
	return *val;
}

int
check_read_pair(struct asm_graph_t *g, struct mini_hash_t *rp_table, struct shortest_path_info_t *r)
{
	int thres_share_read_pair;
	if (r->n_e >= LONG_PATH)
		thres_share_read_pair = MIN_PAIR_SUPPORT_PAIR_END;
	else
		thres_share_read_pair = MIN_PAIR_SUPPORT_PAIR_END_SOFT;

	int count_share_read_pair = 0;
	log_debug("Check from first edge: %d", r->path[0]);
	for (int i = 1; i < r->n_e; i++) {
		uint64_t t = get_share_read_pair(rp_table, r->path[0], g->edges[r->path[i]].rc_id);
		if (g->edges[r->path[i]].seq_len < VERY_SHORT_EDGE_LEN ||
		    g->edges[r->path[0]].seq_len < VERY_SHORT_EDGE_LEN) {
			log_debug("%d fail, length is short %d, thres: %d", r->path[i],
				  g->edges[r->path[i]].seq_len, VERY_SHORT_EDGE_LEN);
			continue;
		}
		if (t > MIN_READ_PAIR_COUNT) {
			count_share_read_pair++;
			log_debug("%d pass, n read-pair support %d, thres: %d", r->path[i], t,
				  MIN_READ_PAIR_COUNT);
		} else
			log_debug("%d fail, not enough read-pair support %d, thres: %d", r->path[i],
				  t, MIN_READ_PAIR_COUNT);
	}
	log_debug("Check to final edge: %d", r->path[r->n_e - 1]);
	for (int i = 0; i < r->n_e - 1; i++) {
		uint64_t t = get_share_read_pair(rp_table, r->path[i],
						 g->edges[r->path[r->n_e - 1]].rc_id);
		if (g->edges[r->path[i]].seq_len < VERY_SHORT_EDGE_LEN ||
		    g->edges[r->path[r->n_e - 1]].seq_len < VERY_SHORT_EDGE_LEN) {
			log_debug("%d fail, length is short %d, thres: %d", r->path[i],
				  g->edges[r->path[i]].seq_len, VERY_SHORT_EDGE_LEN);
			continue;
		}
		if (t > MIN_READ_PAIR_COUNT) {
			count_share_read_pair++;
			log_debug("%d pass, n read-pair support %d, thres: %d", r->path[i], t,
				  MIN_READ_PAIR_COUNT);
		} else
			log_debug("%d fail, not enough read-pair support %d, thres: %d", r->path[i],
				  t, MIN_READ_PAIR_COUNT);
	}
	if (count_share_read_pair > thres_share_read_pair) {
		log_debug("%d-%d passed: %d pairs of edges support. thres: %d", r->path[0],
			  r->path[r->n_e - 1], count_share_read_pair, thres_share_read_pair);
		return 1;
	}
	log_debug("%d-%d fail: %d pairs of edges support. thres: %d", r->path[0],
		  r->path[r->n_e - 1], count_share_read_pair, thres_share_read_pair);
	return 0;
}

void test_rp_table(struct mini_hash_t *rp_table, int n_pairs, int *list_edges)
{
	for (int i = 0; i < n_pairs; i++) {
		int u = list_edges[i * 2];
		int v = list_edges[i * 2 + 1];
		uint64_t res = get_share_read_pair(rp_table, u, v);
		for (int j = 0; j < 5; j++) {
			assert(get_share_read_pair(rp_table, u, v) == res);
		}
	}

}

void print_dot_graph(struct barcode_graph *bg, char *path)
{
	FILE *out = fopen(path, "w");
	fprintf(out, "digraph G {\n");
	for (int i = 0; i < bg->n_edges * 2; i += 2) {
		if (bg->is_del[i]) {
			assert(bg->is_del[i + 1]);
			continue;
		}
		assert(bg->is_del[i + 1] == 0);
		fprintf(out, "%d -> %d\n", bg->edges[i], bg->edges[i + 1]);
	}
	fprintf(out, "}\n");
	fclose(out);
}

void remove_tips_barcode_graph(struct asm_graph_t *g, struct barcode_graph *bg,
			       khash_t(long_spath) *stored)
{
	const int max_deg_out = 10;
	for (int node = 0; node < bg->n_nodes; node++) {
		int *out = calloc(max_deg_out, 4);
		int *out_id = calloc(max_deg_out, 4);
		int deg_out = 0;
		for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
			if (bg->is_del[j] == 0) {
				if (j & 1) {
					out[deg_out % max_deg_out] = bg->edges[j];
					out_id[deg_out % max_deg_out] = j;
					deg_out++;
				}
			}
		if (deg_out > max_deg_out)
			continue;
		int *flag = calloc(deg_out, 4);
		for (int i = 0; i < deg_out; i++)
			if (flag[i] == 0) {
				struct shortest_path_info_t *spath = get_shortest_path(g, node, out[i], stored);
				for (int j = 0; j < spath->n_e; ++j) {
					for (int l = 0; l < deg_out; l++)
						if (spath->path[j] == out[l] && l != i) {
							flag[l] = out[i] + 1;
						}
				}
			}
		for (int j = 0; j < deg_out; ++j) {
			if (flag[j]) {
				del_pair_edge(bg, out_id[j]);
				log_debug("Del remove tips hao %d %d, end %d", node, out[j],
					  flag[j] - 1);
			}
		}
		free(out);
		free(out_id);
		free(flag);
	}

	for (int node = 0; node < bg->n_nodes; node++) {
		int *out = calloc(max_deg_out, 4);
		int *out_id = calloc(max_deg_out, 4);
		int deg_out = 0;
		for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
			if (bg->is_del[j] == 0) {
				if ((j & 1) == 0) {
					out[deg_out % max_deg_out] = bg->edges[j];
					out_id[deg_out % max_deg_out] = j;
					deg_out++;
				}
			}
		if (deg_out > max_deg_out)
			continue;
		int *flag = calloc(deg_out, 4);
		for (int i = 0; i < deg_out; i++)
			if (flag[i] == 0) {
				struct shortest_path_info_t *spath = get_shortest_path(g, out[i], node, stored);
				for (int j = 0; j < spath->n_e; ++j) {
					for (int l = 0; l < deg_out; l++)
						if (spath->path[j] == out[l] && l != i) {
							flag[l] = out[i] + 1;
						}
				}
			}
		for (int j = 0; j < deg_out; ++j) {
			if (flag[j]) {
				del_pair_edge(bg, out_id[j]);
				log_debug("Del remove tips hao %d %d, end %d", out[j], node,
					  flag[j] - 1);
			}
		}
		free(out);
		free(out_id);
		free(flag);
	}
}

void filter_go_reverse_complement(struct asm_graph_t *g, struct barcode_graph *bg)
{
	for (int node = 0; node < bg->n_nodes; node++) {
		for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
			if (bg->is_del[j] == 0) {
				int u = bg->edges[j];
				if (g->edges[node].rc_id == u) {
					del_pair_edge(bg, j);
					log_debug("Del reverse complement %d %d", node, u);
				}
			}
	}
}

void print_del_barcode_graph(struct barcode_graph *bg)
{
	log_debug("log_barcode_graph");
	for (int i = 0; i < bg->n_edges * 2; i += 2) {
		log_debug("%d %d", bg->edges[i], bg->is_del[i]);
		log_debug("%d %d", bg->edges[i + 1], bg->is_del[i + 1]);
	}
}

void filter_shortest_path_and_readpair(struct asm_graph_t *g, struct barcode_graph *bg,
				       khash_t(long_spath) *stored,
				       struct mini_hash_t *rp_table)
{
	int count = 0;
	int del_by_shortest = 0, del_by_too_long = 0, del_by_read_pair = 0;
	for (int i = 0; i < bg->n_edges; i++) {
		if (bg->is_del[i << 1] == 0) {
			int u = bg->edges[i << 1];
			int v = bg->edges[(i << 1) + 1];
			struct shortest_path_info_t *r = get_shortest_path(g, u, v, stored);
			if (r == NULL) {
				del_pair_edge(bg, i << 1);
				log_debug("Del not found shortest path %d %d", u, v);
				del_by_shortest++;
				continue;
			}
			if (r->sum_seq > MAX_RADIUS) {
				del_pair_edge(bg, i << 1);
				log_debug("Del too long path %d %d", u, v);
				del_by_too_long++;
				continue;
			}
			if (!check_read_pair(g, rp_table, r)) {
				del_pair_edge(bg, i << 1);
				log_debug("Del read pair %d %d", u, v);
				del_by_read_pair++;
				continue;
			}
			count++;
		}
	}
	log_info("Del by shortest %d", del_by_shortest);
	log_info("Del by too long %d", del_by_too_long);
	log_info("Del by read pair %d", del_by_read_pair);
	log_info("n edges after filter shortest path and read_pair %d", count);
}

void filter_go_both_reverse_complement(struct asm_graph_t *g, struct barcode_graph *bg)
{
	int *list_adj = calloc(1000, 4);
	int *list_id = calloc(1000, 4);
	for (int node = 0; node < bg->n_nodes; node++) {
		int n_adj = 0;
		for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
			if (bg->is_del[j] == 0 && (j & 1)) {
				list_adj[n_adj] = bg->edges[j];
				list_id[n_adj] = j;
				n_adj++;
			}
		for (int j = 0; j < n_adj; j++) {
			int u = list_adj[j];
			for (int l = j + 1; l < n_adj; l++) {
				int v = list_adj[l];
				if (g->edges[u].rc_id == v) {
					log_debug("Del go both reverse pair %d %d", node, u);
					del_eight_edge(bg, list_id[j]);
				}
			}
		}
	}
	for (int node = 0; node < bg->n_nodes; node++) {
		int n_adj = 0;
		for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
			if (bg->is_del[j] == 0 && (j & 1) == 0) {
				list_adj[n_adj] = bg->edges[j];
				list_id[n_adj] = j;
				n_adj++;
			}
		for (int j = 0; j < n_adj; j++) {
			int u = list_adj[j];
			for (int l = j + 1; l < n_adj; l++) {
				int v = list_adj[l];
				if (g->edges[u].rc_id == v) {
					log_debug("Del go both reverse pair %d %d", u, node);
					del_eight_edge(bg, list_id[j]);
				}
			}
		}
	}
}

void filter_list_edge(struct opt_proc_t *opt, struct mini_hash_t *rp_table, struct asm_graph_t *g,
		      int n_pairs, int *list_edges, int *n_ret, int **list_ret)
{
	khash_t(long_spath) *stored = kh_init(long_spath);

	int n_edges = 0;
	for (int i = 0; i < n_pairs * 2; i++)
		if (n_edges < list_edges[i])
			n_edges = list_edges[i];

	struct barcode_graph *bg = init_for_barcode_graph(n_edges + 1);
	append_list_edge(bg, n_pairs, list_edges);
	//todo @huu test move remove tip to here
	filter_shortest_path_and_readpair(g, bg, stored, rp_table);
	print_dot_graph(bg, "after_filter_BCandpair.dot");

	filter_go_both_reverse_complement(g, bg);
	filter_go_reverse_complement(g, bg);
	print_dot_graph(bg, "after_filter_reverse.dot");

	filter_by_deg(bg, 4);
	remove_tips_barcode_graph(g, bg, stored);
	print_dot_graph(bg, "after_remove_tips.dot");
	long_spath_destroy(stored);

	filter_complex_barcode_graph(bg);
	print_dot_graph(bg, "after_filter_complex.dot");

	filter_bulge(bg);
	print_dot_graph(bg, "after_filter_bulge.dot");

	filter_by_deg(bg, 1);
	print_dot_graph(bg, "after_all.dot");
	print_del_barcode_graph(bg);

	khash_t(set_long) *mark_link = kh_init(set_long);
	int *list_res = NULL, n_res = 0;
	for (int i = 0; i < bg->n_edges * 2; i += 2) {
		if (bg->is_del[i]) {
			assert(bg->is_del[i + 1]);
			continue;
		}
		assert(bg->is_del[i + 1] == 0);
		list_res = realloc(list_res, ((n_res + 1) << 1) * sizeof(int));
		int v = bg->edges[i];
		int u = bg->edges[i + 1];
		uint64_t code = GET_CODE(v, u);
		if (!kh_set_long_exist(mark_link, code))
			kh_set_long_add(mark_link, code);
		else
			log_error("wtf");
		list_res[n_res << 1] = v;
		list_res[(n_res << 1) + 1] = u;
		++n_res;
		log_debug("Add edge to final list %d %d", v, u);
	}
	destroy_barcode_graph(bg);
	kh_destroy(set_long, mark_link);

	*list_ret = list_res;
	*n_ret = n_res;
	log_info("n edges after in deg and out deg: %d", *n_ret);
}

void print_bx_count(khash_t(long_int) *res, struct opt_proc_t *opt)
{
	char path[1024];
	sprintf(path, "%s/bc_hits_edge_pairs.txt", opt->out_dir);
	FILE *fp = fopen(path, "w");
	for (khiter_t k = kh_begin(res); k != kh_end(res); ++k) {
		if (kh_exist(res, k)) {
			fprintf(fp, "%lu %lu %d\n", kh_key(res, k) >> 32,
				kh_key(res, k) & 0x00000000ffffffff, kh_value(res, k));
		}
	}
	fclose(fp);
}

void append_edge(int *n_edges, int **list_edges, int u, int v)
{
	*list_edges = realloc(*list_edges, (((*n_edges) + 1) << 1) * sizeof(int));
	(*list_edges)[*n_edges << 1] = u;
	(*list_edges)[(*n_edges << 1) + 1] = v;
	(*n_edges)++;
}

void save_khash(khash_t(long_int) *h, char *path)
{
	FILE *f = xfopen(path, "wb");
	xfwrite(&h->n_buckets, sizeof(h->n_buckets), 1, f);
	xfwrite(&h->size, sizeof(h->size), 1, f);
	xfwrite(&h->n_occupied, sizeof(h->n_occupied), 1, f);
	xfwrite(&h->upper_bound, sizeof(h->upper_bound), 1, f);

	//flag
	xfwrite(h->flags, sizeof(khint32_t), __ac_fsize(h->n_buckets), f);
	//key
	xfwrite(h->keys, sizeof(uint64_t), h->n_buckets, f);
	//value
	xfwrite(h->vals, sizeof(int), h->n_buckets, f);

	fclose(f);
}

khash_t(long_int) *load_khash(char *path)
{
	khash_t(long_int) *h = kh_init(long_int);
	FILE *f = xfopen(path, "rb");
	xfread(&h->n_buckets, sizeof(h->n_buckets), 1, f);
	xfread(&h->size, sizeof(h->size), 1, f);
	xfread(&h->n_occupied, sizeof(h->n_occupied), 1, f);
	xfread(&h->upper_bound, sizeof(h->upper_bound), 1, f);

	//flag
	h->flags = calloc(__ac_fsize(h->n_buckets), sizeof(khint32_t));
	xfread(h->flags, sizeof(khint32_t), __ac_fsize(h->n_buckets), f);
	//key
	h->keys = calloc(h->n_buckets, sizeof(uint64_t));
	xfread(h->keys, sizeof(uint64_t), h->n_buckets, f);
	//value
	h->vals = calloc(h->n_buckets, sizeof(int));
	xfread(h->vals, sizeof(int), h->n_buckets, f);

	fclose(f);
	return h;
}

void compare_equal(khash_t(long_int) *a, khash_t(long_int) *b)
{
	assert(a->n_buckets == b->n_buckets);
	assert(a->size == b->size);
	assert(a->n_occupied == b->n_occupied);
	assert(a->upper_bound == b->upper_bound);

	for (int i = 0; i < __ac_fsize(a->n_buckets); i++) {
		assert(a->flags[i] == b->flags[i]);
	}
	for (int i = 0; i < a->n_buckets; i++) {
		assert(a->keys[i] == b->keys[i]);
		assert(a->vals[i] == b->vals[i]);
	}
}

void save_simple_minihash(struct mini_hash_t *simple, char *path)
{
	FILE *f = fopen(path, "wb");
	xfwrite(&simple->size, sizeof(uint64_t), 1, f);
	xfwrite(&simple->count, sizeof(uint64_t), 1, f);
	xfwrite(&simple->max_cnt, sizeof(uint64_t), 1, f);
	xfwrite(&simple->prime_index, sizeof(int), 1, f);
	//write key
	xfwrite(simple->key, sizeof(uint64_t), simple->size, f);
	//write h
	xfwrite(simple->h, sizeof(uint64_t), simple->size, f);
	fclose(f);
}

struct mini_hash_t *load_simple_minihash(char *path)
{
	FILE *f = fopen(path, "rb");
	struct mini_hash_t *simple = calloc(1, sizeof(struct mini_hash_t));

	xfread(&simple->size, sizeof(uint64_t), 1, f);
	xfread(&simple->count, sizeof(uint64_t), 1, f);
	xfread(&simple->max_cnt, sizeof(uint64_t), 1, f);
	xfread(&simple->prime_index, sizeof(int), 1, f);
	//read key
	simple->key = calloc(simple->size, sizeof(uint64_t));
	xfread(simple->key, sizeof(uint64_t), simple->size, f);
	//read h
	simple->h = calloc(simple->size, sizeof(uint64_t));
	xfread(simple->h, sizeof(uint64_t), simple->size, f);
	fclose(f);
	return simple;
}

void get_list_contig(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	log_info("get list contig");
	struct mini_hash_t *rp_table = NULL;
	khash_t(long_int) *all_count = NULL;
	if (opt->var == NULL) {
		log_info("save ");
		struct mm_bundle_t *t = mm_hit_all_barcodes(opt, g);
		struct mini_hash_t *bx_table = t->bx_table;
		rp_table = t->rp_table;
		save_simple_minihash(t->rp_table, str_concate(opt->out_dir, "/rptable.mini"));
		all_count = count_edge_link_shared_bc(g, bx_table);
		save_khash(all_count, str_concate(opt->out_dir, "/allcount.khash"));
	} else {
		all_count = load_khash(opt->var[0]);
		rp_table = load_simple_minihash(opt->var[1]);
	}

	print_bx_count(all_count, opt);
	int n_edges = 0;
	int *list_edges = NULL;
	double global_cov = get_genome_coverage_h(g);
	for (int i = kh_begin(all_count); i != kh_end(all_count); i++) {
		if (kh_exist(all_count, i)) {
			uint64_t key = kh_key(all_count, i);
			int u = (key >> 32) & (uint32_t) (-1);
			int v = key & (uint32_t) (-1);
			if (u == v || u == g->edges[u].rc_id) {
				continue;
			}
			double cov_u = __get_edge_cov(&g->edges[u], g->ksize);
			double cov_v = __get_edge_cov(&g->edges[v], g->ksize);
			//todo @huu in metagenomics case, consider relative between u and v
			if (global_cov * 1.5 < cov_u || global_cov * 1.5 < cov_v) {
				continue;
			}
			uint64_t val = kh_value(all_count, i);
			if (g->edges[u].seq_len < MIN_EDGE_LEN ||
			    g->edges[v].seq_len < MIN_EDGE_LEN)
				continue;
			int len_u = MIN(g->edges[u].seq_len, MOLECULE_DENSITY);
			int len_v = MIN(g->edges[v].seq_len, MOLECULE_DENSITY);
			if (val * 1.0 / (len_u + len_v) < MIN_SHARED_BARCODE_RATIO) {
//				log_debug("Edge %d %d, len_u %d, len_v %d, share only %d barcodes", u, v, len_u, len_v, val);
				continue;
			}
			int u_rc = g->edges[u].rc_id;
			int v_rc = g->edges[v].rc_id;
			append_edge(&n_edges, &list_edges, u, v);
			append_edge(&n_edges, &list_edges, v_rc, u_rc);

			append_edge(&n_edges, &list_edges, u, v_rc);
			append_edge(&n_edges, &list_edges, v, u_rc);

			append_edge(&n_edges, &list_edges, u_rc, v);
			append_edge(&n_edges, &list_edges, v_rc, u);

			append_edge(&n_edges, &list_edges, u_rc, v_rc);
			append_edge(&n_edges, &list_edges, v, u);
		}
	}
	kh_destroy(long_int, all_count);

	log_info("n pair pass barcode count filtering %d", n_edges);
	int n_res = 0, *list_res = NULL;
	filter_list_edge(opt, rp_table, g, n_edges, list_edges, &n_res, &list_res);
	free(list_edges);
	destroy_mini_hash(rp_table);

	create_barcode_molecules(opt, list_res, n_res * 2, g);
	free(list_res);

}
