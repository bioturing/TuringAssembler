#include "scaffolding/scaffolding.h"
#include "attribute.h"
#include "cluster_molecules.h"
#include "helper.h"
#include "verbose.h"
#include "assembly_graph.h"
#include "minimizers/count_barcodes.h"
#include "minimizers/smart_load.h"
#include "minimizers/minimizers.h"
#include "log.h"
#include "process.h"
#include "utils.h"

//int get_pair_distance(int v, int u, khash_t(long_spath) *spath_info)
//{
//	uint64_t code = GET_CODE(v, u);
//	if (kh_long_spath_exist(spath_info, code) == 0)
//		return -1;
//
//	return kh_long_spath_get(spath_info, code)->len;
//}

//int check_connected(struct asm_graph_t *g, int v, int u,
//		khash_t(long_spath) *spath_info)
//{
//	int tg = g->edges[v].target;
//	int sr = g->nodes[g->edges[u].source].rc_id;
//	for (int i = 0; i < g->nodes[tg].deg; ++i){
//		for (int j = 0; j < g->nodes[sr].deg; ++j){
//			int w = g->nodes[tg].adj[i];
//			int t = g->edges[g->nodes[sr].adj[j]].rc_id;
//			if (w == u || t == v)
//				return 1;
//			int d = get_pair_distance(w, t, spath_info);
//			if (d == -1)
//				continue;
//			if(d > MAX_RADIUS)
//				log_error("Something went wrong, probably Dijkstra is incorrect");
//			return 1;
//		}
//	}
//	return 0;
//}

//void get_edge_links_by_distance(struct asm_graph_t *g, int *edges, int n_e,
//		khash_t(long_spath) *spath_info, khash_t(long_int) *is_connected,
//		khash_t(long_int) *count_link)
//{
//	for (int i = 0; i < n_e; ++i){
//		for (int j = 0; j < n_e; ++j){
//			if (i == j)
//				continue;
//			int v = edges[i];
//			int u = edges[j];
//			uint64_t code = GET_CODE(v, u);
//			int ok;
//			if (kh_long_int_exist(is_connected, code)){
//				ok = kh_long_int_get(is_connected, code);
//			} else {
//				ok = check_connected(g, v, u, spath_info);
//				kh_long_int_set(is_connected, code, ok);
//			}
//			if (!ok)
//				continue;
//			int cur_count = 0;
//			if (kh_long_int_exist(count_link, code))
//				cur_count = kh_long_int_get(count_link, code);
//			kh_long_int_set(count_link, code, cur_count + 1);
//		}
//	}
//}

//void count_edge_links_bc(struct opt_proc_t *opt)
//{
//	khash_t(long_int) *link_count = kh_init(long_int);
//
//	struct bc_hit_bundle_t bc_hit_bundle;
//	get_bc_hit_bundle(opt, &bc_hit_bundle);
//	struct asm_graph_t *g = bc_hit_bundle.g;
//
//	khash_t(long_spath) *spath_info = kh_init(long_spath);
//	khash_t(long_int) *is_connected = kh_init(long_int);
//	get_all_shortest_paths_dp(bc_hit_bundle.g, spath_info);
//
//	struct barcode_list_t blist;
//	get_barcode_list(opt->bx_str, &blist);
//
//	FILE *bc_log = fopen("bc_log.txt", "w");
//	for (int i = 0; i < blist.n_bc; ++i){
//		if ((i + 1) % 10000 == 0)
//			log_debug("%d/%d barcodes processed", i + 1, blist.n_bc);
//		if (blist.read_count[i] < MIN_BC_READ_COUNT
//			|| blist.read_count[i] > MAX_BC_READ_COUNT)
//			continue;
//
//		struct mm_hits_t *hits = get_hits_from_barcode(blist.bc_list[i],
//				&bc_hit_bundle);
//
//		int *edges;
//		int n_e;
//		hits_to_edges(g, hits, &edges, &n_e);
//		mm_hits_destroy(hits);
//
//		get_edge_links_by_distance(bc_hit_bundle.g, edges, n_e, spath_info,
//				is_connected, link_count);
//
//		fprintf(bc_log, "%s: ", blist.bc_list[i]);
//		for (int i = 0; i < n_e; ++i)
//			fprintf(bc_log, "%d,", edges[i]);
//		fprintf(bc_log, "\n");
//		free(edges);
//	}
//	fclose(bc_log);
//
//	log_info("Writing all shortest paths");
//	FILE *all_paths = fopen("all_shortest_paths.txt", "w");
//	khash_t(set_long) *mark = kh_init(set_long);
//	for (int v = 0; v < g->n_e; ++v){
//		int tg = g->edges[v].target;
//		for (int i = 0; i < g->nodes[tg].deg; ++i){
//			int u = g->nodes[tg].adj[i];
//			kh_set_long_insert(mark, GET_CODE(v, u));
//			fprintf(all_paths, "%d to %d: %d,%d\n", v, u, v, u);
//		}
//	}
//	for (khiter_t it = kh_begin(spath_info); it != kh_end(spath_info); ++it){
//		if (!kh_exist(spath_info, it))
//			continue;
//		uint64_t code = kh_key(spath_info, it);
//		struct shortest_path_info_t *wrapper = kh_val(spath_info, it);
//		int v = code >> 32;
//		int u = code & ((uint32_t) -1);
//
//		int source = g->nodes[g->edges[v].source].rc_id;
//		int target = g->edges[u].target;
//		for (int i = 0; i < g->nodes[source].deg; ++i){
//			int w = g->edges[g->nodes[source].adj[i]].rc_id;
//			for (int j = 0; j < g->nodes[target].deg; ++j){
//				int t = g->nodes[target].adj[j];
//				uint64_t new_code = GET_CODE(w, t);
//				if (kh_set_long_exist(mark, new_code))
//					continue;
//				kh_set_long_insert(mark, new_code);
//				fprintf(all_paths, "%d to %d: ", w, t);
//				fprintf(all_paths, "%d,", w);
//				int p = v;
//				while (p != -1){
//					fprintf(all_paths, "%d,", p);
//					uint64_t new_code = GET_CODE(p, u);
//					p = kh_long_spath_get(spath_info, new_code)->trace;
//				}
//				fprintf(all_paths, "%d\n", t);
//
//			}
//		}
//	}
//	kh_destroy(set_long, mark);
//	fclose(all_paths);
//
//	barcode_list_destroy(&blist);
//	bc_hit_bundle_destroy(&bc_hit_bundle);
//	kh_destroy(long_spath, spath_info);
//	kh_destroy(long_int, is_connected);
//
//
//
//	FILE *f = fopen(opt->lc, "w");
//	for (khiter_t it = kh_begin(link_count); it != kh_end(link_count); ++it){
//		if (!kh_exist(link_count, it))
//			continue;
//		uint64_t code = kh_key(link_count, it);
//		int count = kh_val(link_count, it);
//		int u = code >> 32;
//		int v = code & ((((uint64_t) 1) << 32) - 1);
//		if (count != -1)
//			fprintf(f, "%d %d %d\n", u, v, count);
//	}
//	fclose(f);
//	kh_destroy(long_int, link_count);
//}

void load_pair_edge_count(char *path, khash_t(long_int) *h_all)
{
	FILE *f = fopen(path, "r");
	int u, v, c;
	while (fscanf(f, "%d %d %d\n", &v, &u, &c) == 3) {
		uint64_t code = GET_CODE(v, u);
		kh_long_int_set(h_all, code, c);
	}
	fclose(f);
}

void print_barcode_graph(struct opt_proc_t *opt)
{
	khash_t(long_int) *h_all = kh_init(long_int);
	load_pair_edge_count(opt->in_fasta, h_all);

	struct bc_hit_bundle_t bc_hit_bundle;
	get_bc_hit_bundle(opt, &bc_hit_bundle);
	struct asm_graph_t *g = bc_hit_bundle.g;

	struct mm_hits_t *hits = get_hits_from_barcode(opt->bx_str, &bc_hit_bundle);

	int *edges;
	int n_e;
	hits_to_edges(g, hits, &edges, &n_e);

	FILE *f = fopen(opt->lc, "w");
	fprintf(f, "digraph %s{\n", opt->bx_str);
	float unit_cov = get_genome_coverage(g);
	for (int i = 0; i < n_e; ++i) {
		for (int j = 0; j < n_e; ++j) {
			if (i == j)
				continue;
			int v = edges[i];
			int u = edges[j];
			if (v == g->edges[u].rc_id && is_repeat(g, v) == 0)
				continue;
			uint64_t code = GET_CODE(v, u);
			int val = kh_long_int_try_get(h_all, code, 0);
			if (val < opt->thresh)
				continue;
			fprintf(f, "\t%d -> %d [label=\"%d\"]\n", v, u, val);
		}
	}

	for (int i = 0; i < n_e; ++i) {
		int v = edges[i];
		float cov = __get_edge_cov(g->edges + v, g->ksize);
		float ratio = cov / unit_cov;
		if (ratio >= 0.8 && ratio <= 1.2)
			fprintf(f, "%d [style=\"filled\",fillcolor=green]\n",
				v);
		else
			fprintf(f, "%d [style=\"filled\",fillcolor=violet]\n",
				v);
	}
	fprintf(f, "}");
	fclose(f);
	free(edges);
}

void get_barcode_list(char *bc_count_path, struct barcode_list_t *blist)
{
	int n = 0;
	int m = 1;
	char **bc_list = calloc(1, sizeof(char *));
	int *read_count = calloc(1, sizeof(int));
	char bc[19];
	int count;
	FILE *f = fopen(bc_count_path, "r");
	while (fscanf(f, "%s\t%d\n", bc, &count) == 2) {
		if (n == m) {
			m <<= 1;
			bc_list = realloc(bc_list, sizeof(char *) * m);
			read_count = realloc(read_count, sizeof(int) * m);
		}
		bc_list[n] = calloc(19, sizeof(char));
		memcpy(bc_list[n], bc, sizeof(char) * 19);
		read_count[n] = count;
		++n;
	}
	fclose(f);
	bc_list = realloc(bc_list, sizeof(char *) * n);
	read_count = realloc(read_count, sizeof(int) * n);
	blist->n_bc = n;
	blist->bc_list = bc_list;
	blist->read_count = read_count;
}

void add_simple_node(struct simple_graph_t *sg, int v)
{
	if (kh_int_node_exist(sg->nodes, v) == 0) {
		struct simple_node_t *snode = calloc(1,
						     sizeof(struct simple_node_t));
		kh_int_node_set(sg->nodes, v, snode);
	}
}

void add_simple_edge(struct simple_graph_t *sg, int v, int u)
{
	struct simple_node_t *snode = kh_int_node_get(sg->nodes, v);
	snode->adj = realloc(snode->adj, sizeof(int) * (snode->deg + 1));
	snode->adj[snode->deg] = u;
	++snode->deg;

	snode = kh_int_node_get(sg->nodes, u);
	snode->rv_adj = realloc(snode->rv_adj, sizeof(int) * (snode->rv_deg + 1));
	snode->rv_adj[snode->rv_deg] = v;
	++snode->rv_deg;
}

void init_simple_graph(struct asm_graph_t *g, struct simple_graph_t *sg)
{
	sg->g = g;
	sg->nodes = kh_init(int_node);
	sg->is_loop = kh_init(set_int);
	sg->is_complex = kh_init(set_int);
	sg->path_len = kh_init(int_int);
	sg->next = kh_init(int_int);
}

void build_simple_graph(int *edges, int n_e, khash_t(long_int) *all_bc,
			struct simple_graph_t *sg)
{
	struct asm_graph_t *g = sg->g;
	for (int i = 0; i < n_e; ++i)
		add_simple_node(sg, edges[i]);
	for (int i = 0; i < n_e; ++i) {
		for (int j = 0; j < n_e; ++j) {
			if (i == j)
				continue;
			int v = edges[i];
			int u = edges[j];
			uint64_t code = GET_CODE(v, u);
			uint64_t code_rc = GET_CODE(g->edges[u].rc_id,
						    g->edges[v].rc_id);
			int val = 0;
			val = max(val, kh_long_int_try_get(all_bc, code, 0));
			val = max(val, kh_long_int_try_get(all_bc, code_rc, 0));

			if (val < MIN_BARCODE_EDGE_COUNT)
				continue;
			add_simple_edge(sg, v, u);
		}
	}
}

void simple_graph_destroy(struct simple_graph_t *sg)
{
	for (khiter_t it = kh_begin(sg->nodes); it != kh_end(sg->nodes); ++it) {
		if (!kh_exist(sg->nodes, it))
			continue;
		struct simple_node_t *snode = kh_val(sg->nodes, it);
		free(snode->adj);
		free(snode->rv_adj);
		free(snode);
	}
	kh_destroy(int_node, sg->nodes);
	kh_destroy(set_int, sg->is_loop);
	kh_destroy(int_int, sg->path_len);
	kh_destroy(int_int, sg->next);
}

void check_loop_dfs(struct simple_graph_t *sg, int v, khash_t(set_int) *visited,
		    khash_t(set_int) *in_dfs)
{
	if (kh_set_int_exist(in_dfs, v)) {
		kh_set_int_insert(sg->is_loop, v);
		return;
	}
	if (kh_set_int_exist(visited, v))
		return;
	kh_set_int_insert(visited, v);
	kh_set_int_insert(in_dfs, v);
	struct simple_node_t *snode = kh_int_node_get(sg->nodes, v);
	for (int i = 0; i < snode->deg; ++i) {
		int u = snode->adj[i];
		check_loop_dfs(sg, u, visited, in_dfs);
	}
	kh_set_int_erase(in_dfs, v);
}

void find_DAG(struct simple_graph_t *sg)
{
	struct asm_graph_t *g = sg->g;
	khash_t(int_node) *nodes = sg->nodes;
	khash_t(set_int) *visited = kh_init(set_int);
	khash_t(set_int) *in_dfs = kh_init(set_int);
	for (khiter_t it = kh_begin(nodes); it != kh_end(nodes); ++it) {
		if (!kh_exist(nodes, it))
			continue;
		int v = kh_key(nodes, it);
		check_loop_dfs(sg, v, visited, in_dfs);
	}
	kh_destroy(set_int, visited);
	kh_destroy(set_int, in_dfs);
}

void get_longest_path_dfs(struct simple_graph_t *sg, int v,
			  khash_t(set_int) *done_dfs)
{
	if (kh_set_int_exist(done_dfs, v))
		return;
	struct simple_node_t *snode = kh_int_node_get(sg->nodes, v);
	int max_len = 0;
	int next = -1;
	for (int i = 0; i < snode->deg; ++i) {
		int u = snode->adj[i];
		get_longest_path_dfs(sg, u, done_dfs);
		int next_len = kh_int_int_get(sg->path_len, u);
		if (max_len < next_len) {
			max_len = next_len;
			next = u;
		}
	}
	kh_int_int_set(sg->path_len, v, max_len + 1);
	kh_int_int_set(sg->next, v, next);
	kh_set_int_insert(done_dfs, v);
}

void get_longest_path(struct simple_graph_t *sg)
{
	khash_t(set_int) *done_dfs = kh_init(set_int);
	for (khiter_t it = kh_begin(sg->nodes); it != kh_end(sg->nodes); ++it) {
		if (!kh_exist(sg->nodes, it))
			continue;
		int v = kh_key(sg->nodes, it);
		if (kh_set_int_exist(sg->is_complex, v))
			continue;
		get_longest_path_dfs(sg, v, done_dfs);
	}
	kh_destroy(set_int, done_dfs);
}

void filter_complex_regions(struct simple_graph_t *sg)
{
	struct asm_graph_t *g = sg->g;
	khash_t(int_node) *nodes = sg->nodes;
	khash_t(set_int) *visited = kh_init(set_int);
	int n_total = 0;
	int n_complex = 0;
	int n_simple_edges = 0;
	for (khiter_t it = kh_begin(nodes); it != kh_end(nodes); ++it){
		if (!kh_exist(nodes, it))
			continue;
		int v = kh_key(nodes, it);
		if (kh_set_int_exist(visited, v))
			continue;

		struct queue_t q;
		init_queue(&q, 1024);
		push_queue(&q, pointerize(&v, sizeof(int)));
		kh_set_int_insert(visited, v);

		khash_t(set_int) *component = kh_init(set_int);
		int has_rc = 0;
		int has_loop = 0;
		int n_source = 0;
		int n_sink = 0;
		while (!is_queue_empty(&q)) {
			int v = *(int *) get_queue(&q);
			free(get_queue(&q));
			pop_queue(&q);

			struct simple_node_t *snode = kh_int_node_get(sg->nodes,
								      v);
			if (snode->deg == 0)
				++n_sink;
			if (snode->rv_deg == 0)
				++n_source;
			if (kh_set_int_exist(component, g->edges[v].rc_id))
				has_rc = 1;
			if (kh_set_int_exist(sg->is_loop, v))
				has_loop = 1;
			kh_set_int_insert(component, v);

			for (int i = 0; i < snode->deg + snode->rv_deg; ++i) {
				int u = i < snode->deg ? snode->adj[i] :
					snode->rv_adj[i - snode->deg];
				if (kh_set_int_exist(visited, u))
					continue;
				kh_set_int_insert(visited, u);
				push_queue(&q, pointerize(&u, sizeof(int)));
			}
		}
		destroy_queue(&q);
		++n_total;

		if (!has_rc && !has_loop && n_source == 1 && n_sink == 1
				&& kh_size(component) > 1){
			n_simple_edges += kh_size(component);
			goto free_stage_1;
		}
		++n_complex;
		for (khiter_t it = kh_begin(component); it != kh_end(component);
		     ++it) {
			if (!kh_exist(component, it))
				continue;
			int v = kh_key(component, it);
			kh_set_int_insert(sg->is_complex, v);
		}
free_stage_1:
		kh_destroy(set_int, component);
	}
	kh_destroy(set_int, visited);
	log_info("Number of all regions: %d", n_total);
	log_info("Number of simple regions: %d", n_total - n_complex);
	log_info("Number of edges in simple regions: %d", n_simple_edges);
}

void print_simple_graph(struct simple_graph_t *sg, int *edges, int n_e, FILE *f)
{
	for (int i = 0; i < n_e; ++i) {
		int v = edges[i];
		if (kh_set_int_exist(sg->is_complex, v))
			continue;
		struct simple_node_t *snode = kh_int_node_get(sg->nodes, v);
		for (int j = 0; j < snode->deg; ++j) {
			int u = snode->adj[j];
			if (kh_set_int_exist(sg->is_complex, u))
				continue;
			fprintf(f, "%d %d\n", v, u);
		}
	}
}

/**
 * Get the shortest path that goes from v to u such that:
 * 	+ The total length of the sequences between v and u (not inclusive) <= MAX_RADIUS
 * 	+ The number of edges between v and u (not inclusive) <= MAX_PATH_LEN
 * @param g: The assembly graph
 * @param spath: The result structure from the function get_all_shortest_paths_dp
 * @param v, u: The query
 * @param path: The result path
 * @param n_v: path len
 */
//int extract_shortest_path(struct asm_graph_t *g, khash_t(long_spath) *spath,
//		int v, int u, int **path, int *n_v)
//{
//	int v_tg = g->edges[v].target;
//	for (int i = 0; i < g->nodes[v_tg].deg; ++i){
//		int w = g->nodes[v_tg].adj[i];
//		if (w == u){
//			*n_v = 2;
//			*path = calloc(2, sizeof(int));
//			(*path)[0] = v;
//			(*path)[1] = u;
//			return 0;
//		}
//	}
//
//	int u_sr = g->nodes[g->edges[u].source].rc_id;
//	int w = -1;
//	int t = -1;
//	int min_len = 1e9;
//
//	for (int i = 0; i < g->nodes[v_tg].deg; ++i){
//		int tmpw = g->nodes[v_tg].adj[i];
//		for (int j = 0; j < g->nodes[u_sr].deg; ++j){
//			int tmpt = g->edges[g->nodes[u_sr].adj[j]].rc_id;
//			uint64_t code = GET_CODE(tmpw, tmpt);
//			if (kh_long_spath_exist(spath, code) == 0)
//				continue;
//			int len = kh_long_spath_get(spath, code)->len;
//			if (min_len > len){
//				min_len = len;
//				w = tmpw;
//				t = tmpt;
//			}
//		}
//	}
//
//	if (w == -1){
//		*path = NULL;
//		*n_v = 0;
//		return -1;
//	}
//	*path = calloc(MAX_PATH_LEN + 2, sizeof(int));
//	*n_v = 1;
//	(*path)[0] = v;
//	int i = w;
//	while (i != -1){
//		(*path)[(*n_v)++] = i;
//		uint64_t code = GET_CODE(i, t);
//		i = kh_long_spath_get(spath, code)->trace;
//	}
//	(*path)[(*n_v)++] = u;
//	return min_len;
//}

//void fill_gap(struct asm_graph_t *g, int v, int u, khash_t(long_spath) *spath,
//		struct simple_graph_t *sg, char **seq)
//{
//	int len = strlen(*seq);
//	int *path;
//	int n_v;
//	extract_shortest_path(g, spath, v, u, &path, &n_v);
//	for (int i = 1; i < n_v; ++i){
//		int w = path[i];
//		char *tmp;
//		decode_seq(&tmp, g->edges[w].seq, g->edges[w].seq_len);
//		*seq = realloc(*seq, len + g->edges[w].seq_len - g->ksize + 1);
//		strcpy(*seq + len, tmp + g->ksize);
//		len += g->edges[w].seq_len - g->ksize;
//		free(tmp);
//	}
//	free(path);
//}

int check_ignore_path(struct asm_graph_t *g, double global_avg_cov,
		int *contig_path, int n_contig_path, double *local_cov)
{
	int ksize = g->ksize;
	double sum_cov = 0, sum_len = 0;
	for (int i = 0; i < n_contig_path; i++) {
		int i_e = contig_path[i];
		double cov = __get_edge_cov(&g->edges[i_e], g->ksize);
		assert(cov >= 0);
		if (cov < MIN_COVERAGE_TO_BE_IGNORE * global_avg_cov){
			log_debug("Ignore path cov: %.3f", cov);
			return -1;
		}
		if (cov > COVERAGE_RATIO_TO_BE_REPEAT * global_avg_cov)
			continue;
		sum_cov += (g->edges[i_e].seq_len-ksize) *cov;
		sum_len += (g->edges[i_e].seq_len-ksize);
	}
	if (sum_len == 0) {
//		assert(sum_len != 0 && "all edge in path is repeat");
		*local_cov = global_avg_cov;
	} else
		*local_cov = sum_cov/sum_len;

	if (*local_cov < MIN_COVERAGE_TO_BE_IGNORE * global_avg_cov) {
		log_debug("Ignore path because seem to be already consider");
		log_debug("local_cov: %.2f, ratio: %.2f, avg_cov: %.2f",
				*local_cov, MIN_COVERAGE_TO_BE_IGNORE,
				global_avg_cov);
		return -1;
	}
	return 0;
}

void concate_edges_fill_N(struct asm_graph_t *g, int *path, int n,
		khash_t(long_spath) *stored, struct asm_edge_t *e)
{
	char *seq;
	decode_seq(&seq, g->edges[path[0]].seq, g->edges[path[0]].seq_len);
	int total_len = strlen(seq);

	int n_holes = 0;
	uint32_t *p_holes = NULL;
	uint32_t *l_holes = NULL;
	int total_count = g->edges[path[0]].count;
	for(int i = 1; i < n; i++) {
		int a = path[i-1], b = path[i];
		struct shortest_path_info_t *res = get_shortest_path(g, a, b, stored);
		if (res == NULL)
			log_error("Can not find shortest path between %d and %d!", a, b);
		if (res->n_e > 2){
			p_holes = realloc(p_holes, (n_holes + 1) * sizeof(uint32_t));
			l_holes = realloc(l_holes, (n_holes + 1) * sizeof(uint32_t));
			int N_len = 0;
			for (int j = 1; j < res->n_e - 1; ++j){
				int e = res->path[j];
				N_len += g->edges[e].seq_len - g->ksize;
			}
			p_holes[n_holes] = total_len;
			l_holes[n_holes] = N_len;
			++n_holes;
		}
		int b_len = g->edges[b].seq_len;
		double tmp_cov = __get_edge_cov(&g->edges[b], g->ksize);
		log_debug("cov %d %d %lf", b, g->edges[b].count, tmp_cov);
		total_count += g->edges[b].count;

		seq = realloc(seq, (total_len + b_len - g->ksize + 1) * sizeof(char));
		char *tmp_seq;
		decode_seq(&tmp_seq, g->edges[b].seq, b_len);
		strcat(seq, tmp_seq + g->ksize);
		total_len += b_len - g->ksize;
		free(tmp_seq);
	}

	uint32_t *seq_encode;
	encode_seq(&seq_encode, seq);
	free(seq);

	e->seq = seq_encode;
	e->seq_len = total_len;
	e->count = total_count;
	log_debug("total_count %d %lf", total_count, __get_edge_cov(e, g->ksize));
	e->n_holes = n_holes;
	e->l_holes = l_holes;
	e->p_holes = p_holes;
}

void add_path_to_edges(struct asm_graph_t *g, struct asm_graph_t *g_new,
		khash_t(long_spath) *stored, int *contig_path, int n_contig_path,
		double avg_cov, int *visited)
{
	double local_cov;
	int ret = check_ignore_path(g, avg_cov, contig_path, n_contig_path,
			&local_cov);
	if (ret == -1)
		return;

	/*
	 * Print the path
	 */
	char path[8192];
	char ctg[16];
	path[0] = '\0';
	for (int i = 0; i < n_contig_path; ++i) {
		sprintf(ctg, "%d ", contig_path[i]);
		strcat(path, ctg);
	}
	log_debug("Path: %s", path);

	struct asm_edge_t new_edge;
	concate_edges_fill_N(g, contig_path, n_contig_path, stored, &new_edge);

	for (int i = 0; i < n_contig_path; i++) {
		int i_e  = contig_path[i];
		int i_e_rc = g->edges[i_e].rc_id;
		visited[i_e]++;
		visited[i_e_rc]++;
		int tmp = MIN((g->edges[i_e].seq_len - g->ksize) * local_cov, g->edges[i_e].count);
		log_debug("Visited and decrease count of %d: from %d to %d", i_e, g->edges[i_e].count, g->edges[i_e].count - tmp);
		g->edges[i_e].count -= tmp;
		log_debug("Visited and decrease count of %d: from %d to %d", i_e_rc, g->edges[i_e_rc].count, g->edges[i_e_rc].count - tmp);
		g->edges[i_e_rc].count -= tmp;
	}


	g_new->edges = realloc(g_new->edges, (g_new->n_e+2) * sizeof(struct asm_edge_t));
	g_new->edges[g_new->n_e] = new_edge;
	asm_clone_seq_reverse(g_new->edges + g_new->n_e + 1,
			g_new->edges + g_new->n_e);
	assert(g_new->edges[g_new->n_e].count == g_new->edges[g_new->n_e+1].count);
	g_new->n_e += 2;
}

void create_barcode_molecules(struct opt_proc_t *opt, int *edges, int n_e,
		struct asm_graph_t *g)
{
	struct asm_graph_t *g_new = calloc(1, sizeof(struct asm_graph_t));
	g_new->ksize = g->ksize;

	struct paths_bundle_t *paths_bundle = calloc(1, sizeof(struct paths_bundle_t));
	get_all_longest_paths(edges, n_e, g, paths_bundle);
	double global_cov = get_genome_coverage_h(g);
	log_info("Global genome coverage %.2f", global_cov);
	int *visited = calloc(g->n_e, sizeof(int));
	log_info("Init visited hash table");
	khash_t(long_spath) *stored = kh_init(long_spath);
	for (int i = 0; i < paths_bundle->n_paths; ++i) {
		add_path_to_edges(g, g_new, stored, paths_bundle->paths[i].edges,
				paths_bundle->paths[i].n_e, global_cov, visited);
	}
	long_spath_destroy(stored);
	paths_bundle_destroy(paths_bundle);

	for (int i = 0; i < g->n_e; ++i){
		int e = i;
		int e_rc = g->edges[e].rc_id;
		int been_touch = visited[e] + visited[e_rc];
		if (__get_edge_cov(&g->edges[i], g->ksize) <= MIN_COVERAGE_TO_BE_IGNORE * global_cov && been_touch != 0) {
			log_debug("Ignore edge %d, coverage %.2f, ration threshold %.2f, coverage threshold %.2f", i, __get_edge_cov(&g->edges[i], g->ksize), MIN_COVERAGE_TO_BE_IGNORE, global_cov);
			continue;
		}
		if (g->edges[i].count == 0)
			continue;
		assert(g->edges[i].count > 0);
		if (e > e_rc) {
			continue;
		}
		g_new->edges = realloc(g_new->edges, (g_new->n_e+2) * sizeof(struct asm_edge_t));
		g_new->edges[g_new->n_e].seq = g->edges[i].seq;
		g_new->edges[g_new->n_e].seq_len = g->edges[i].seq_len;
		g_new->edges[g_new->n_e].count = g->edges[i].count;
		g_new->edges[g_new->n_e].n_holes = 0;
		asm_clone_seq_reverse(g_new->edges + g_new->n_e + 1,
				g_new->edges + g_new->n_e);
		assert(g_new->edges[g_new->n_e].count == g_new->edges[g_new->n_e+1].count);
		g_new->n_e += 2;
	}
	free(visited);

	for (int i = 0; i < g_new->n_e; ++i){
		g_new->edges[i].source = i << 1;
		g_new->edges[i].target = (i << 1) + 1;
		g_new->edges[i].rc_id = i ^ 1;
	}
	g_new->n_v = g_new->n_e * 2;
	g_new->nodes = calloc(g_new->n_e * 2, sizeof(struct asm_node_t));
	for (int e = 0; e < g_new->n_e; ++e){
		int e_rc = g_new->edges[e].rc_id;
		int v = g_new->edges[e].source;
		int u = g_new->edges[e].target;
		int v_rc = g_new->edges[e_rc].source;
		int u_rc = g_new->edges[e_rc].target;

		g_new->nodes[v].rc_id = u_rc;
		g_new->nodes[v].deg = 1;
		g_new->nodes[v].adj = calloc(1, sizeof(gint_t));
		g_new->nodes[v].adj[0] = e;

		g_new->nodes[u].rc_id = v_rc;
	}
	for (int i = 0; i < g_new->n_e; ++i){
		log_debug("Edge %d len %d count %d", i, g_new->edges[i].seq_len,
				g_new->edges[i].count);
	}

	save_graph_info(opt->out_dir, g_new, "level_3");
	asm_graph_destroy(g_new);
	free(g_new);
}

void barcode_list_destroy(struct barcode_list_t *blist)
{
	for (int i = 0; i < blist->n_bc; ++i)
		free(blist->bc_list[i]);
	free(blist->bc_list);
	free(blist->read_count);
}

int is_repeat(struct asm_graph_t *g, int e)
{
	float unit_cov = get_genome_coverage(g);
	float cov = __get_edge_cov(g->edges + e, g->ksize);
	float ratio = cov / unit_cov;
	if (ratio > 1.2)
		return 1;
	return 0;
}

void print_graph_component(struct simple_graph_t *sg, char *bc, FILE *f)
{
	struct asm_graph_t *g = sg->g;
	khash_t(set_int) *visited = kh_init(set_int);
	for (khiter_t it = kh_begin(sg->nodes); it != kh_end(sg->nodes); ++it) {
		if (!kh_exist(sg->nodes, it))
			continue;
		int s = kh_key(sg->nodes, it);
		if (kh_set_int_exist(visited, s))
			continue;

		int has_rc = 0;
		int has_loop = 0;
		khash_t(set_int) *component = kh_init(set_int);

		struct queue_t q;
		init_queue(&q, 1024);
		push_queue(&q, pointerize(&s, sizeof(int)));
		kh_set_int_insert(visited, s);
		while (!is_queue_empty(&q)) {
			int v = *(int *) get_queue(&q);
			free(get_queue(&q));
			pop_queue(&q);

			if (kh_set_int_exist(component, g->edges[v].rc_id))
				has_rc = 1;
			if (kh_set_int_exist(sg->is_loop, v))
				has_loop = 1;
			kh_set_int_insert(component, v);

			struct simple_node_t *node = kh_int_node_get(sg->nodes,
								     v);
			for (int i = 0; i < node->deg; ++i) {
				int u = node->adj[i];
				if (kh_set_int_exist(visited, u))
					continue;
				kh_set_int_insert(visited, u);
				push_queue(&q, pointerize(&u, sizeof(int)));
			}

			node = kh_int_node_get(sg->nodes, g->edges[v].rc_id);
			for (int i = 0; i < node->deg; ++i) {
				int u = g->edges[node->adj[i]].rc_id;
				if (kh_int_node_exist(sg->nodes, u) == 0)
					continue;
				if (kh_set_int_exist(visited, u))
					continue;
				kh_set_int_insert(visited, u);
				push_queue(&q, pointerize(&u, sizeof(int)));
			}
		}
		if (!has_loop && !has_rc && kh_size(component) > 1) {
			fprintf(f, "digraph %s{\n", bc);
			for (khiter_t it = kh_begin(component); it != kh_end(component);
			     ++it) {
				if (!kh_exist(component, it))
					continue;
				int s = kh_key(component, it);
				struct simple_node_t *snode = kh_int_node_get(
					sg->nodes, s);
				for (int i = 0; i < snode->deg; ++i)
					fprintf(f, "\t%d -> %d\n", s, snode->adj[i]);
				float unit_cov = get_genome_coverage(sg->g);
				float cov = __get_edge_cov(g->edges + s, g->ksize);
				float ratio = cov / unit_cov;
				if (ratio >= 0.8 && ratio <= 1.2)
					fprintf(f, "\t%d [style=\"filled\",fillcolor=green]\n",
						s);
				else
					fprintf(f, "\t%d [style=\"filled\",fillcolor=violet]\n",
						s);
			}
			fprintf(f, "}\n");
		}
		kh_destroy(set_int, component);
		destroy_queue(&q);
	}
	kh_destroy(set_int, visited);
}

void hits_to_edges(struct asm_graph_t *g, struct mm_hits_t *hits, int **edges,
		   int *n_e)
{
	khash_t(set_int) *tmp = kh_init(set_int);
	for (khiter_t it = kh_begin(hits->edges); it != kh_end(hits->edges); ++it) {
		if (!kh_exist(hits->edges, it))
			continue;
		int e = kh_key(hits->edges, it);
		int rc = g->edges[e].rc_id;
		kh_set_int_insert(tmp, e);
		kh_set_int_insert(tmp, rc);
	}

	*edges = calloc(kh_size(tmp), sizeof(int));
	*n_e = 0;
	for (khiter_t it = kh_begin(tmp); it != kh_end(tmp); ++it) {
		if (!kh_exist(tmp, it))
			continue;
		int e = kh_key(tmp, it);
		(*edges)[(*n_e)++] = e;
	}
	kh_destroy(set_int, tmp);
}

/**
 * Get all shortest paths between all the edges in the assembly graph.
 * 	Every path p goes from v to u satisfies these two conditions:
 * 		+ The sum of sequences between v and u (inclusive) is <= MAX_RADIUS
 * 		+ The number of edges between v and u (inclusive) is <= MAX_PATH_LEN
 * 	The result is stored in the structure khash_t(long_spath) *spath_info, inwhich
 * 		keys are uint64_t integers that encode v and u (v = key >> 32, u = key & -1)
 * 		and values are pairs of values (len, trace) where len is the length of the
 * 		shortest path from v to u and trace is the next edge on that path.
 * 		If trace is -1 indicates the path's end.
 * @param g: the assembly graph
 * @param spath_info the result structure
 */
//void get_all_shortest_paths_dp(struct asm_graph_t *g, khash_t(long_spath) *spath_info)
//{
//	khash_t(long_int) *L_pre = kh_init(long_int);
//	for (int i = 0; i < g->n_e; ++i){
//		int v = i;
//		if (g->edges[v].seq_len > MAX_RADIUS)
//			continue;
//		kh_long_int_set(L_pre, GET_CODE(v, v), g->edges[v].seq_len);
//
//		struct shortest_path_info_t *wrapper = calloc(1,
//				sizeof(struct shortest_path_info_t));
//		wrapper->len = g->edges[v].seq_len;
//		wrapper->trace = -1;
//
//		kh_long_spath_set(spath_info, GET_CODE(v, v), wrapper);
//	}
//	for (int i = 2; i <= MAX_PATH_LEN; ++i){
//		khash_t(long_int) *L_cur = kh_init(long_int);
//		for (int j = 0; j < g->n_e; ++j){
//			if ((j + 1) % 10000 == 0)
//				log_debug("%d-th iteration: %d edges processed",
//						i, j + 1);
//			int v = j;
//			if (g->edges[v].seq_len > MAX_RADIUS)
//				continue;
//			int *edges;
//			int n_e;
//			bfs_nearby(g, v, i, &edges, &n_e);
//			for (int k = 0; k < n_e; ++k){
//				int u = edges[k];
//				if (g->edges[u].seq_len > MAX_RADIUS)
//					continue;
//				uint64_t code = GET_CODE(v, u);
//				int tg = g->edges[v].target;
//				int min_len = 1e9;
//				int next = -1;
//				for (int h = 0; h < g->nodes[tg].deg; ++h){
//					int w = g->nodes[tg].adj[h];
//					uint64_t code = GET_CODE(w, u);
//					if (kh_long_int_exist(L_pre, code)
//						== 0)
//						continue;
//					int new_len = kh_long_int_get(L_pre, code)
//						+ g->edges[v].seq_len - g->ksize;
//					if (min_len > new_len){
//						min_len = new_len;
//						next = w;
//					}
//				}
//				if (min_len > MAX_RADIUS)
//					continue;
//				kh_long_int_set(L_cur, code, min_len);
//
//				struct shortest_path_info_t *wrapper;
//				if (kh_long_spath_exist(spath_info, code)){
//					wrapper = kh_long_spath_get(spath_info, code);
//				} else {
//					wrapper = calloc(1, sizeof(struct shortest_path_info_t));
//					wrapper->len = 1e9;
//					wrapper->trace = 0;
//					kh_long_spath_set(spath_info, code, wrapper);
//				}
//				if (wrapper->len > min_len){
//					wrapper->len = min_len;
//					wrapper->trace = next;
//				}
//			}
//			free(edges);
//		}
//		kh_destroy(long_int, L_pre);
//		L_pre = L_cur;
//	}
//	kh_destroy(long_int, L_pre);
//}

void bfs_nearby(struct asm_graph_t *g, int source, int radius, int **edges, int *n_e)
{
	khash_t(int_int) *L = kh_init(int_int);
	struct queue_t q;
	init_queue(&q, 1024);
	push_queue(&q, pointerize(&source, sizeof(int)));
	kh_int_int_set(L, source, 1);
	while (!is_queue_empty(&q)) {
		int v = *(int *) get_queue(&q);
		free(get_queue(&q));
		pop_queue(&q);

		int len = kh_int_int_get(L, v);
		if (len > radius)
			break;
		int tg = g->edges[v].target;
		for (int i = 0; i < g->nodes[tg].deg; ++i) {
			int u = g->nodes[tg].adj[i];
			if (kh_int_int_exist(L, u))
				continue;
			push_queue(&q, pointerize(&u, sizeof(int)));
			kh_int_int_set(L, u, len + 1);
		}
	}
	free_queue_content(&q);
	destroy_queue(&q);

	*edges = calloc(kh_size(L), sizeof(int));
	*n_e = 0;
	for (khiter_t it = kh_begin(L); it != kh_end(L); ++it) {
		if (!kh_exist(L, it))
			continue;
		(*edges)[(*n_e)++] = kh_key(L, it);
	}
	kh_destroy(int_int, L);
}

//void get_all_pair_edges(struct asm_graph_t *g, khash_t(long_int) *pair_edges)
//{
//	for (int v = 0; v < g->n_e; ++v){
//		int target = g->edges[v].target;
//		for (int i = 0; i < g->nodes[target].deg; ++i){
//			int u = g->nodes[target].adj[i];
//			uint64_t code = GET_CODE(v, u);
//			kh_long_int_set(pair_edges, code, 0);
//		}
//	}
//
//	khash_t(long_spath) *spath = kh_init(long_spath);
//	get_all_shortest_paths_dp(g, spath);
//	for (khiter_t it = kh_begin(spath); it != kh_end(spath); ++it){
//		if (!kh_exist(spath, it))
//			continue;
//		uint64_t code = kh_key(spath, it);
//		int len = kh_val(spath, it)->len;
//		if (len > MAX_RADIUS)
//			log_error("Something went wrong when finding all the sohrtest paths");
//		int v = code >> 32;
//		int u = (uint32_t) code & -1;
//		int source = g->nodes[g->edges[v].source].rc_id;
//		int target = g->edges[u].target;
//		for (int i = 0; i < g->nodes[source].deg; ++i){
//			int w = g->edges[g->nodes[source].adj[i]].rc_id;
//			for (int j = 0; j < g->nodes[target].deg; ++j){
//				int t = g->nodes[target].adj[j];
//				uint64_t old_code = GET_CODE(w, t);
//				int old_len = kh_long_int_try_get(pair_edges,
//						old_code, 1e9);
//				if (old_len > len)
//					kh_long_int_set(pair_edges, old_code,
//							len);
//			}
//		}
//	}
//	kh_destroy(long_spath, spath);
//}

void build_graph_from_edges_list(int *edges, int n_e, struct asm_graph_t *g,
				 struct simple_graph_t *sg)
{
	init_simple_graph(g, sg);
	for (int i = 0; i < n_e; i += 2) {
		int v = edges[i];
		int u = edges[i + 1];
		add_simple_node(sg, v);
		add_simple_node(sg, u);
		add_simple_edge(sg, v, u);
	}
}

void get_all_longest_paths(int *edges, int n_e, struct asm_graph_t *g,
			   struct paths_bundle_t *path_bundle)
{
	struct simple_graph_t sg;
	build_graph_from_edges_list(edges, n_e, g, &sg);
	find_DAG(&sg);
	filter_complex_regions(&sg);
	get_longest_path(&sg);


	float unit_cov = get_genome_coverage(g);
	for (khiter_t it = kh_begin(sg.nodes); it != kh_end(sg.nodes); ++it) {
		if (!kh_exist(sg.nodes, it))
			continue;
		int source = kh_key(sg.nodes, it);
		struct simple_node_t *snode = kh_val(sg.nodes, it);
		if (kh_set_int_exist(sg.is_complex, source))
			continue;
		if (snode->rv_deg != 0)
			continue;
		if (__get_edge_cov(&g->edges[source], g->ksize) <= 0.5 * unit_cov)
			continue;

		int *path = calloc(1, sizeof(int));
		int n_e = 1;
		path[0] = source;
		for (int v = kh_int_int_get(sg.next, source); v != -1;
		     v = kh_int_int_get(sg.next, v)) {
			path = realloc(path, (n_e + 1) * sizeof(int));
			path[n_e] = v;
			n_e++;
		}
		path_bundle->paths = realloc(path_bundle->paths,
					     sizeof(struct simple_path_t) *
					     (path_bundle->n_paths + 1));
		struct simple_path_t new_path = {
			.edges = path,
			.n_e = n_e
		};
		path_bundle->paths[path_bundle->n_paths] = new_path;
		++path_bundle->n_paths;
	}
}

void paths_bundle_destroy(struct paths_bundle_t *paths_bundle)
{
	for (int i = 0; i < paths_bundle->n_paths; ++i)
		free(paths_bundle->paths[i].edges);
	free(paths_bundle->paths);
}

int check_adj_edges(struct asm_graph_t *g, int v, int u)
{
	int v_tg = g->edges[v].target;
	for (int i = 0; i < g->nodes[v_tg].deg; ++i) {
		int w = g->nodes[v_tg].adj[i];
		if (w == u)
			return 1;
	}
	return 0;
}

struct shortest_path_info_t *get_shortest_path(struct asm_graph_t *g, int s,
					       int t, khash_t(long_spath) *stored)
{
	if (kh_long_spath_exist(stored, GET_CODE(s, t)))
		return kh_long_spath_get(stored, GET_CODE(s, t));

	struct shortest_path_info_t *res = NULL;
	if (check_adj_edges(g, s, t)) {
		res = calloc(1, sizeof(struct shortest_path_info_t));
		res->sum_seq = 0;
		res->n_e = 2;
		res->path = calloc(2, sizeof(int));
		res->path[0] = s;
		res->path[1] = t;
		kh_long_spath_set(stored, GET_CODE(s, t), res);
		return res;
	}

	struct queue_t q;
	init_queue(&q, 1024);
	khash_t(int_int) *best_L = kh_init(int_int);
	khash_t(long_int) *best_P = kh_init(long_int);
	khash_t(int_int) *best_D = kh_init(int_int);

	int s_tg = g->edges[s].target;
	struct len_info_t wrapper;
	for (int i = 0; i < g->nodes[s_tg].deg; ++i) {
		int v = g->nodes[s_tg].adj[i];
		wrapper.v = v;
		wrapper.sum_seq = g->edges[v].seq_len;
		wrapper.len = 1;
		push_queue(&q, pointerize(&wrapper, sizeof(struct len_info_t)));
		kh_long_int_set(best_P, GET_CODE(v, 1), -1);
		kh_int_int_set(best_L, v, 0);
		kh_int_int_set(best_D, v, 1);
	}

	while (!is_queue_empty(&q)) {
		struct len_info_t *wrapper = get_queue(&q);
		int v = wrapper->v;
		int sum_seq = wrapper->sum_seq;
		int len = wrapper->len;
		free(wrapper);
		pop_queue(&q);
		if (len == MAX_PATH_LEN)
			continue;

		if (g->edges[v].seq_len > MIN_EDGE_LEN)
			continue;
		int v_tg = g->edges[v].target;
		for (int i = 0; i < g->nodes[v_tg].deg; ++i) {
			int u = g->nodes[v_tg].adj[i];
			int new_sum_seq = sum_seq + g->edges[u].seq_len - g->ksize;
			int new_len = len + 1;
			if (new_sum_seq < kh_int_int_try_get(best_L, u, 1e9)) {
				kh_int_int_set(best_L, u, new_sum_seq);
				kh_long_int_set(best_P, GET_CODE(u, new_len), v);
				kh_int_int_set(best_D, u, new_len);

				struct len_info_t new_wrapper = {
					.v = u,
					.sum_seq = new_sum_seq,
					.len = new_len
				};
				push_queue(&q, pointerize(&new_wrapper,
							  sizeof(struct len_info_t)));
			}
		}
	}
	destroy_queue(&q);

	int best_w = -1;
	int best_sum_seq = 1e9;
	for (khiter_t it = kh_begin(best_L); it != kh_end(best_L); ++it) {
		if (!kh_exist(best_L, it))
			continue;
		int w = kh_key(best_L, it);
		int sum_seq = kh_val(best_L, it);
		if (sum_seq > best_sum_seq)
			continue;

		int w_tg = g->edges[w].target;
		for (int i = 0; i < g->nodes[w_tg].deg; ++i) {
			if (t == g->nodes[w_tg].adj[i]) {
				best_w = w;
				best_sum_seq = sum_seq;
				break;
			}
		}
	}

	if (best_w == -1)
		return NULL;

	int *path = calloc(MAX_PATH_LEN + 2, sizeof(int));
	int n_e = 1;
	path[0] = t;

	int d = kh_int_int_get(best_D, best_w);
	int v = best_w;
	while (v != -1) {
		path[n_e++] = v;
		v = kh_long_int_get(best_P, GET_CODE(v, d));
		--d;
	}
	path[n_e++] = s;
	for (int i = 0, j = n_e - 1; i < j; ++i, --j) {
		int tmp = path[i];
		path[i] = path[j];
		path[j] = tmp;
	}

	res = calloc(1, sizeof(struct shortest_path_info_t));
	res->sum_seq = best_sum_seq;
	res->n_e = n_e;
	res->path = path;
	kh_long_spath_set(stored, GET_CODE(s, t), res);

	kh_destroy(int_int, best_L);
	kh_destroy(long_int, best_P);
	kh_destroy(int_int, best_D);
	return res;
}

void long_spath_destroy(khash_t(long_spath) *stored)
{
	for (khiter_t it = kh_begin(stored); it != kh_end(stored); ++it){
		if (!kh_exist(stored, it))
			continue;
		struct shortest_path_info_t *spath = kh_val(stored, it);
		free(spath->path);
	}
	kh_destroy(long_spath, stored);
}
