//
// Created by che on 04/12/2019.
//

#include <stdlib.h>
#include <minimizers/count_barcodes.h>
#include "cluster_molecules.h"
#include "log.h"
#include "utils.h"

#define MIN_SHARE_BARCODE_COUNT 100
#define MIN_EDGE_LEN 500
#define MIN_READ_PAIR_COUNT 1
#define VERY_SHORT_EDGE_LEN 250
#define LONG_PATH 10
#define SHORT_PATH 2
#define MIN_PAIR_SUPPORT_PAIR_END 1
#define MIN_PAIR_SUPPORT_PAIR_END_SOFT 0

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
						bg->is_del[j] = 1;
						bg->is_del[j ^ 1] = 1;
						log_debug("Del bulge %d %d", node, near[1]);
					}
			} else if (have_edge(bg, near[1], near[0]))
				for (int j = bg->first_index[node];
				     j != -1; j = bg->next_index[j])
					if (bg->edges[j] == near[0]) {
						bg->is_del[j] = 1;
						bg->is_del[j ^ 1] = 1;
						log_debug("Del bulge %d %d", node, near[0]);
					}
		}
	}
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
					bg->is_del[j] = 1;
					bg->is_del[j ^ 1] = 1;
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
					bg->is_del[j] = 1;
					bg->is_del[j ^ 1] = 1;
				}
		}
	}
}

void filter_complex_barcode_graph(struct barcode_graph *bg)
{
	filter_bulge(bg);
//	filter_by_deg(bg, 3);
//	filter_bulge(bg);
//	filter_by_deg(bg, 1);
}

uint64_t get_share_read_pair(struct mini_hash_t *rp_table, int u, int v)
{
	//todo @huu
	uint64_t t = GET_CODE(u, v);
	uint64_t *val = mini_get(rp_table, t);
	if (val == (uint64_t *)EMPTY_BX)
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
		if (g->edges[r->path[i]].seq_len < VERY_SHORT_EDGE_LEN || g->edges[r->path[0]].seq_len < VERY_SHORT_EDGE_LEN) {
			log_debug("%d fail, length is short %d, thres: %d", r->path[i], g->edges[r->path[i]].seq_len, VERY_SHORT_EDGE_LEN);
			continue;
		}
		if (t > MIN_READ_PAIR_COUNT) {
			count_share_read_pair++;
			log_debug("%d pass, n read-pair support %d, thres: %d", r->path[i], t, MIN_READ_PAIR_COUNT);
		} else
			log_debug("%d fail, not enough read-pair support %d, thres: %d", r->path[i], t, MIN_READ_PAIR_COUNT);
	}
	log_debug("Check to final edge: %d", r->path[r->n_e - 1]);
	for (int i = 0; i < r->n_e - 1; i++) {
		uint64_t t = get_share_read_pair(rp_table, r->path[i],
						 g->edges[r->path[r->n_e - 1]].rc_id);
		if (g->edges[r->path[i]].seq_len < VERY_SHORT_EDGE_LEN ||
		    g->edges[r->path[r->n_e - 1]].seq_len < VERY_SHORT_EDGE_LEN) {
			log_debug("%d fail, length is short %d, thres: %d", r->path[i], g->edges[r->path[i]].seq_len, VERY_SHORT_EDGE_LEN);
			continue;
		}
		if (t > MIN_READ_PAIR_COUNT) {
			count_share_read_pair++;
			log_debug("%d pass, n read-pair support %d, thres: %d", r->path[i], t, MIN_READ_PAIR_COUNT);
		} else
			log_debug("%d fail, not enough read-pair support %d, thres: %d", r->path[i], t, MIN_READ_PAIR_COUNT);
	}
	if (count_share_read_pair > thres_share_read_pair) {
		log_debug("%d-%d passed: %d pairs of edges support. thres: %d", r->path[0],r->path[r->n_e - 1], count_share_read_pair, thres_share_read_pair);
		return 1;
	}
	log_debug("%d-%d fail: %d pairs of edges support. thres: %d", r->path[0],r->path[r->n_e- 1], count_share_read_pair, thres_share_read_pair);
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
	for (int node = 0; node < bg->n_nodes; node++) {
		int out[2];
		int out_id[2];
		int deg_out = 0;
		for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
			if (bg->is_del[j] == 0) {
				if (j & 1) {
					out[deg_out & 1] = bg->edges[j];
					out_id[deg_out & 1] = j;
					deg_out++;
				}
			}

		if (deg_out != 2)
			continue;
		struct shortest_path_info_t *spath = get_shortest_path(g, node,
								       out[0], stored);
		int flag[2] = {};
		for (int i = 0; i < spath->n_e; ++i) {
			if (spath->path[i] == out[1]) {
				flag[1] = 1;
				break;
			}
		}

		spath = get_shortest_path(g, node, out[1], stored);
		for (int i = 0; i < spath->n_e; ++i) {
			if (spath->path[i] == out[0]) {
				flag[0] = 1;
				break;
			}
		}
		if (!flag[0] && !flag[1])
			flag[0] = flag[1] = 1;
		for (int i = 0; i < 2; ++i) {
			if (flag[i]) {
				bg->is_del[out_id[i]] = 1;
				bg->is_del[out_id[i] ^ 1] = 1;
				log_debug("Del remove tips hao %d %d", node, out[i]);
			}
		}
	}


	for (int node = 0; node < bg->n_nodes; ++node) {
		int in[2];
		int in_id[2];
		int deg_in = 0;
		for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
			if (bg->is_del[j] == 0) {
				if ((j & 1) == 0) {
					in[deg_in & 1] = bg->edges[j];
					in_id[deg_in & 1] = j;
					deg_in++;
				}
			}

		if (deg_in != 2)
			continue;

		struct shortest_path_info_t *spath = get_shortest_path(g, in[0],
								       node, stored);
		int flag[2] = {};
		for (int i = 0; i < spath->n_e; ++i) {
			if (spath->path[i] == in[1]) {
				flag[1] = 1;
				break;
			}
		}

		spath = get_shortest_path(g, in[1], node, stored);
		for (int i = 0; i < spath->n_e; ++i) {
			if (spath->path[i] == in[0]) {
				flag[0] = 1;
				break;
			}
		}
		if (!flag[0] && !flag[1])
			flag[0] = flag[1] = 1;
		for (int i = 0; i < 2; ++i) {
			if (flag[i]) {
				bg->is_del[in_id[i]] = 1;
				bg->is_del[in_id[i] ^ 1] = 1;
				log_debug("Del remove tips hao %d %d", in[i], node);
			}
		}
	}
}

void filter_graph_reverse_complement(struct asm_graph_t *g, struct barcode_graph *bg)
{
	for (int node = 0; node < bg->n_nodes; node++) {
		for (int j = bg->first_index[node]; j != -1; j = bg->next_index[j])
			if (bg->is_del[j] == 0) {
				int u = bg->edges[j];
				if (g->edges[node].rc_id == u) {
					bg->is_del[j] = 1;
					bg->is_del[j ^ 1] = 1;
					log_debug("Del reverse complement %d %d", node, u);
				}
			}
	}
}

void filter_list_edge(struct opt_proc_t *opt, struct mini_hash_t *rp_table, struct asm_graph_t *g,
		      int n_pairs, int *list_edges, int *n_ret, int **list_ret)
{
	test_rp_table(rp_table, n_pairs, list_edges);
	khash_t(long_spath) *stored = kh_init(long_spath);

	int n_filted = 0, *list_filted = NULL;
	for (int i = 0; i < n_pairs; i++) {
		int u = list_edges[i << 1];
		int v = list_edges[(i << 1) + 1];
		struct shortest_path_info_t *r = get_shortest_path(g, u, v, stored);
		if (r == NULL) {
			log_debug("Del not found shortest path %d %d", u, v);
			continue;
		}
		if (r->sum_seq > MAX_RADIUS) {
			log_debug("Del too long path %d %d", u, v);
			continue;
		}
		if (!check_read_pair(g, rp_table, r)) {
			log_debug("Del read pair %d %d", u, v);
			continue;
		}
		list_filted = realloc(list_filted, (n_filted + 1) * 2 * sizeof(int));
		list_filted[n_filted << 1] = u;
		list_filted[(n_filted << 1) + 1] = v;
		n_filted++;
	}

	log_info("n edges after filter shortest path and read_pair %d", n_filted);
	int n_edges = 0;
	for (int i = 0; i < n_filted; i++)
		if (n_edges < list_filted[i])
			n_edges = list_filted[i];

	struct barcode_graph *bg = init_for_barcode_graph(n_edges + 1);
	append_list_edge(bg, n_filted, list_filted);
	print_dot_graph(bg, "after_filter_BCandpair.dot");
	filter_graph_reverse_complement(g, bg);
	remove_tips_barcode_graph(g, bg, stored);
	print_dot_graph(bg, "after_remove_tips.dot");
	filter_complex_barcode_graph(bg);

	khash_t(set_long) *mark_link = kh_init(long_int);
	print_dot_graph(bg, "after_filter_complex.dot");
	int *list_res = NULL, n_res = 0;
	for (int i = 0; i < bg->n_edges * 2; i += 2) {
		if (bg->is_del[i]) {
			assert(bg->is_del[i + 1]);
			continue;
		}
		assert(bg->is_del[i + 1] == 0);
		list_res = realloc(list_res, ((n_res + 2) << 1) * sizeof(int));
		int v = bg->edges[i];
		int u = bg->edges[i + 1];
		int v_rc = g->edges[v].rc_id;
		int u_rc = g->edges[u].rc_id;
		uint64_t code = GET_CODE(v, u);
		uint64_t code_rc = GET_CODE(u_rc, v_rc);
		if (kh_set_long_exist(mark_link, code)
				|| kh_set_long_exist(mark_link, code_rc))
			continue;
		kh_set_long_add(mark_link, code);
		kh_set_long_add(mark_link, code_rc);
		list_res[n_res << 1] = v;
		list_res[(n_res << 1) + 1] = u;
		++n_res;
		list_res[n_res << 1] = u_rc;
		list_res[(n_res << 1) + 1] = v_rc;
		++n_res;
		log_debug("%d %d", v, u);
		log_debug("%d %d", u_rc, v_rc);
	}
	kh_destroy(set_long, mark_link);

	*list_ret = list_res;
	*n_ret = n_res;
	log_info("n edges after in deg and out deg: %d", *n_ret);
	destroy_barcode_graph(bg);
}

void print_bx_count(khash_t(long_int) *res, struct opt_proc_t *opt)
{
	char path[1024];
	sprintf(path, "%s/bc_hits_edge_pairs.txt", opt->out_dir);
	FILE *fp = fopen(path, "w");
	for (khiter_t k = kh_begin(res); k != kh_end(res); ++k) {
		if (kh_exist(res, k)) {
			fprintf(fp, "%lu %lu %d\n", kh_key(res, k) >> 32, kh_key(res, k) & 0x00000000ffffffff, kh_value(res, k));
		}
	}
	fclose(fp);
}

void get_list_contig(struct opt_proc_t *opt, struct asm_graph_t *g)
{
	struct mm_bundle_t *t = mm_hit_all_barcodes(opt);
	struct mini_hash_t *bx_table = t->bx_table;
	struct mini_hash_t *rp_table = t->rp_table;
	khash_t(long_int) *all_count = count_edge_link_shared_bc(g, bx_table);

	print_bx_count(all_count, opt);

	int n_edges = 0;
	int *list_edges = NULL;
	for (int i = kh_begin(all_count); i != kh_end(all_count); i++) {
		if (kh_exist(all_count, i)) {
			uint64_t key = kh_key(all_count, i);
			int u = (key >> 32) & (uint32_t) (-1);
			int v = key & (uint32_t) (-1);
			uint64_t val = kh_value(all_count, i);
			if (g->edges[u].seq_len < MIN_EDGE_LEN || g->edges[v].seq_len < MIN_EDGE_LEN)
				continue;
			int len_u = MIN(g->edges[u].seq_len, 3000);
			int len_v = MIN(g->edges[v].seq_len, 3000);
			if (val * 1.0 / (len_u + len_v) < 0.033) {
				continue;
			}
			list_edges = realloc(list_edges, ((n_edges + 1) << 1) * sizeof(int));
			list_edges[n_edges << 1] = u;
			list_edges[(n_edges << 1) + 1] = v;
			n_edges++;
			list_edges = realloc(list_edges, ((n_edges + 1) << 1) * sizeof(int));
			list_edges[n_edges << 1] = g->edges[v].rc_id;
			list_edges[(n_edges << 1) + 1] = g->edges[u].rc_id;
			n_edges++;
		}
	}
	log_info("n pair share 100 bc %d", n_edges);
	int n_res = 0, *list_res = NULL;
	filter_list_edge(opt, rp_table, g, n_edges, list_edges, &n_res, &list_res);
	create_barcode_molecules(opt, list_res, n_res * 2, g);
}
