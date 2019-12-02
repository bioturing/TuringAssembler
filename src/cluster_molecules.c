#include "attribute.h"
#include "cluster_molecules.h"
#include "helper.h"
#include "verbose.h"
#include "assembly_graph.h"
#include "minimizers/count_barcodes.h"
#include "minimizers/smart_load.h"
#include "minimizers/minimizers.h"
#include "log.h"

#define MAX_RADIUS 4000
#define MAX_PATH_LEN 30
#define MIN_BC_READ_COUNT 10
#define MAX_BC_READ_COUNT 88
#define MIN_BARCODE_EDGE_COUNT 100
#define GET_CODE(a, b) ((((uint64_t) (a)) << 32) | (b))

int get_pair_distance(int v, int u, khash_t(long_spath) *spath_info)
{
	uint64_t code = GET_CODE(v, u);
	if (kh_long_spath_exist(spath_info, code) == 0)
		return -1;
	return kh_long_spath_get(spath_info, code)->len;
}

int check_connected(struct asm_graph_t *g, int v, int u,
		khash_t(long_spath) *spath_info)
{
	int tg = g->edges[v].target;
	int sr = g->nodes[g->edges[u].source].rc_id;
	for (int i = 0; i < g->nodes[tg].deg; ++i){
		for (int j = 0; j < g->nodes[sr].deg; ++j){
			int w = g->nodes[tg].adj[i];
			int t = g->edges[g->nodes[sr].adj[j]].rc_id;
			if (w == u || t == v)
				return 1;
			int d = get_pair_distance(w, t, spath_info);
			if (d == -1)
				continue;
			if(d > MAX_RADIUS)
				log_error("Something went wrong, probably Dijkstra is incorrect");
			return 1;
		}
	}
	return 0;
}

void get_edge_links_by_distance(struct asm_graph_t *g, int *edges, int n_e,
		khash_t(long_spath) *spath_info, khash_t(long_int) *is_connected,
		khash_t(long_int) *count_link)
{
	for (int i = 0; i < n_e; ++i){
		for (int j = 0; j < n_e; ++j){
			if (i == j)
				continue;
			int v = edges[i];
			int u = edges[j];
			uint64_t code = GET_CODE(v, u);
			int ok;
			if (kh_long_int_exist(is_connected, code)){
				ok = kh_long_int_get(is_connected, code);
			} else {
				ok = check_connected(g, v, u, spath_info);
				kh_long_int_set(is_connected, code, ok);
			}
			if (!ok)
				continue;
			int cur_count = 0;
			if (kh_long_int_exist(count_link, code))
				cur_count = kh_long_int_get(count_link, code);
			kh_long_int_set(count_link, code, cur_count + 1);
		}
	}
}

void count_edge_links_bc(struct opt_proc_t *opt)
{
	khash_t(long_int) *link_count = kh_init(long_int);

	struct bc_hit_bundle_t bc_hit_bundle;
	get_bc_hit_bundle(opt, &bc_hit_bundle);
	struct asm_graph_t *g = bc_hit_bundle.g;

	khash_t(long_spath) *spath_info = kh_init(long_spath);
	khash_t(long_int) *is_connected = kh_init(long_int);
	get_all_shortest_paths_dp(bc_hit_bundle.g, spath_info);

	struct barcode_list_t blist;
	get_barcode_list(opt->bx_str, &blist);

	FILE *bc_log = fopen("bc_log.txt", "w");
	for (int i = 0; i < blist.n_bc; ++i){
		if ((i + 1) % 10000 == 0)
			log_debug("%d/%d barcodes processed", i + 1, blist.n_bc);
		if (blist.read_count[i] < MIN_BC_READ_COUNT
			|| blist.read_count[i] > MAX_BC_READ_COUNT)
			continue;

		struct mm_hits_t *hits = get_hits_from_barcode(blist.bc_list[i],
				&bc_hit_bundle);

		int *edges;
		int n_e;
		hits_to_edges(g, hits, &edges, &n_e);
		mm_hits_destroy(hits);

		/*get_edge_links_by_distance(bc_hit_bundle.g, edges, n_e, spath_info,
				is_connected, link_count);*/

		fprintf(bc_log, "%s: ", blist.bc_list[i]);
		for (int i = 0; i < n_e; ++i)
			fprintf(bc_log, "%d%c", edges[i], i + 1 == n_e ? '\n' : ',');
		free(edges);
	}
	fclose(bc_log);

	log_info("Writing all shortest paths");
	FILE *all_paths = fopen("all_shortest_paths.txt", "w");
	khash_t(set_long) *mark = kh_init(set_long);
	for (int v = 0; v < g->n_e; ++v){
		int tg = g->edges[v].target;
		for (int i = 0; i < g->nodes[tg].deg; ++i){
			int u = g->nodes[tg].adj[i];
			kh_set_long_add(mark, GET_CODE(v, u));
			fprintf(all_paths, "%d to %d: %d,%d\n", v, u, v, u);
		}
	}
	for (khiter_t it = kh_begin(spath_info); it != kh_end(spath_info); ++it){
		if (!kh_exist(spath_info, it))
			continue;
		uint64_t code = kh_key(spath_info, it);
		struct shortest_path_info_t *wrapper = kh_val(spath_info, it);
		int v = code >> 32;
		int u = code & ((uint32_t) -1);

		int source = g->nodes[g->edges[v].source].rc_id;
		int target = g->edges[u].target;
		for (int i = 0; i < g->nodes[source].deg; ++i){
			int w = g->edges[g->nodes[source].adj[i]].rc_id;
			for (int j = 0; j < g->nodes[target].deg; ++j){
				int t = g->nodes[target].adj[j];
				uint64_t new_code = GET_CODE(w, t);
				if (kh_set_long_exist(mark, new_code))
					continue;
				kh_set_long_add(mark, new_code);
				fprintf(all_paths, "%d to %d: ", w, t);
				fprintf(all_paths, "%d,", w);
				int p = v;
				while (p != -1){
					fprintf(all_paths, "%d,", p);
					uint64_t new_code = GET_CODE(p, u);
					p = kh_long_spath_get(spath_info, new_code)->trace;
				}
				fprintf(all_paths, "%d\n", t);

			}
		}
	}
	kh_destroy(set_long, mark);
	fclose(all_paths);

	barcode_list_destroy(&blist);
	bc_hit_bundle_destroy(&bc_hit_bundle);
	kh_destroy(long_spath, spath_info);
	kh_destroy(long_int, is_connected);



	FILE *f = fopen(opt->lc, "w");
	for (khiter_t it = kh_begin(link_count); it != kh_end(link_count); ++it){
		if (!kh_exist(link_count, it))
			continue;
		uint64_t code = kh_key(link_count, it);
		int count = kh_val(link_count, it);
		int u = code >> 32;
		int v = code & ((((uint64_t) 1) << 32) - 1);
		if (count != -1)
			fprintf(f, "%d %d %d\n", u, v, count);
	}
	fclose(f);
	kh_destroy(long_int, link_count);
}

void load_pair_edge_count(char *path, khash_t(long_int) *h_all)
{
	FILE *f = fopen(path, "r");
	int u, v, c;
	while (fscanf(f, "%d %d %d\n", &v, &u, &c) == 3){
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
	for (int i = 0; i < n_e; ++i){
		for (int j = 0; j < n_e; ++j){
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

	for (int i = 0; i < n_e; ++i){
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
	while (fscanf(f, "%s\t%d\n", bc, &count) == 2){
		if (n == m){
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
	if (kh_int_node_exist(sg->nodes, v) == 0){
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
	for (int i = 0; i < n_e; ++i){
		for (int j = 0; j < n_e; ++j){
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
	for (khiter_t it = kh_begin(sg->nodes); it != kh_end(sg->nodes); ++it){
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
	if (kh_set_int_exist(in_dfs, v)){
		kh_set_int_add(sg->is_loop, v);
		return;
	}
	if (kh_set_int_exist(visited, v))
		return;
	kh_set_int_add(visited, v);
	kh_set_int_add(in_dfs, v);
	struct simple_node_t *snode = kh_int_node_get(sg->nodes, v);
	for (int i = 0; i < snode->deg; ++i){
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
	for (khiter_t it = kh_begin(nodes); it != kh_end(nodes); ++it){
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
	for (int i = 0; i < snode->deg; ++i){
		int u = snode->adj[i];
		get_longest_path_dfs(sg, u, done_dfs);
		int next_len = kh_int_int_get(sg->path_len, u);
		if (max_len < next_len){
			max_len = next_len;
			next = u;
		}
	}
	kh_int_int_set(sg->path_len, v, max_len + 1);
	kh_int_int_set(sg->next, v, next);
	kh_set_int_add(done_dfs, v);
}

void get_longest_path(struct simple_graph_t *sg)
{
	khash_t(set_int) *done_dfs = kh_init(set_int);
	for (khiter_t it = kh_begin(sg->nodes); it != kh_end(sg->nodes); ++it){
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
	for (khiter_t it = kh_begin(nodes); it != kh_end(nodes); ++it){
		if (!kh_exist(nodes, it))
			continue;
		int v = kh_key(nodes, it);
		if (kh_set_int_exist(visited, v))
			continue;

		struct queue_t q;
		init_queue(&q, 1024);
		push_queue(&q, pointerize(&v, sizeof(int)));
		kh_set_int_add(visited, v);

		khash_t(set_int) *component = kh_init(set_int);
		int has_rc = 0;
		int has_loop = 0;
		int n_source = 0;
		int n_sink = 0;
		while (!is_queue_empty(&q)){
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
			kh_set_int_add(component, v);

			for (int i = 0; i < snode->deg + snode->rv_deg; ++i){
				int u = i < snode->deg ? snode->adj[i] :
					snode->rv_adj[i - snode->deg];
				if (kh_set_int_exist(visited, u))
					continue;
				kh_set_int_add(visited, u);
				push_queue(&q, pointerize(&u, sizeof(int)));
			}
		}
		if (!has_rc && !has_loop && n_source == 1 && n_sink == 1
				&& kh_size(component) > 1)
			continue;
		for (khiter_t it = kh_begin(component); it != kh_end(component);
				++it){
			if (!kh_exist(component, it))
				continue;
			int v = kh_key(component, it);
			kh_set_int_add(sg->is_complex, v);
		}
		kh_destroy(set_int, component);
	}
	kh_destroy(set_int, visited);
}

void print_simple_graph(struct simple_graph_t *sg, int *edges, int n_e, FILE *f)
{
	for (int i = 0; i < n_e; ++i){
		int v = edges[i];
		if (kh_set_int_exist(sg->is_complex, v))
			continue;
		struct simple_node_t *snode = kh_int_node_get(sg->nodes, v);
		for (int j = 0; j < snode->deg; ++j){
			int u = snode->adj[j];
			if (kh_set_int_exist(sg->is_complex, u))
				continue;
			fprintf(f, "%d %d\n", v, u);
		}
	}
}

void create_barcode_molecules(struct opt_proc_t *opt)
{
	struct bc_hit_bundle_t bc_hit_bundle;
	get_bc_hit_bundle(opt, &bc_hit_bundle);
	struct asm_graph_t *g = bc_hit_bundle.g;

	struct barcode_list_t blist;
	get_barcode_list(opt->bx_str, &blist);

	khash_t(long_int) *all_pairs = kh_init(long_int);
	load_pair_edge_count(opt->in_fasta, all_pairs);

	int *edges = calloc(g->n_e, sizeof(int));
	int n_e = g->n_e;
	for (int i = 0; i < n_e; ++i)
		edges[i] = i;

	struct simple_graph_t sg;
	init_simple_graph(g, &sg);
	build_simple_graph(edges, n_e, all_pairs, &sg);
	find_DAG(&sg);

	filter_complex_regions(&sg);
	//get_longest_path(&sg);

	FILE *f = fopen(opt->lc, "w");
	print_simple_graph(&sg, edges, n_e, f);
	fclose(f);

	simple_graph_destroy(&sg);


	free(edges);

	kh_destroy(long_int, all_pairs);
	barcode_list_destroy(&blist);
	bc_hit_bundle_destroy(&bc_hit_bundle);

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
	for (khiter_t it = kh_begin(sg->nodes); it != kh_end(sg->nodes); ++it){
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
		kh_set_int_add(visited, s);
		while (!is_queue_empty(&q)){
			int v = *(int *) get_queue(&q);
			free(get_queue(&q));
			pop_queue(&q);

			if (kh_set_int_exist(component, g->edges[v].rc_id))
				has_rc = 1;
			if (kh_set_int_exist(sg->is_loop, v))
				has_loop = 1;
			kh_set_int_add(component, v);

			struct simple_node_t *node = kh_int_node_get(sg->nodes,
					v);
			for (int i = 0; i < node->deg; ++i){
				int u = node->adj[i];
				if (kh_set_int_exist(visited, u))
					continue;
				kh_set_int_add(visited, u);
				push_queue(&q, pointerize(&u, sizeof(int)));
			}

			node = kh_int_node_get(sg->nodes, g->edges[v].rc_id);
			for (int i = 0; i < node->deg; ++i){
				int u = g->edges[node->adj[i]].rc_id;
				if (kh_int_node_exist(sg->nodes, u) == 0)
					continue;
				if (kh_set_int_exist(visited, u))
					continue;
				kh_set_int_add(visited, u);
				push_queue(&q, pointerize(&u, sizeof(int)));
			}
		}
		if (!has_loop && !has_rc && kh_size(component) > 1){
			fprintf(f, "digraph %s{\n", bc);
			for (khiter_t it = kh_begin(component); it != kh_end(component);
					++it){
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
	for (khiter_t it = kh_begin(hits->edges); it != kh_end(hits->edges); ++it){
		if (!kh_exist(hits->edges, it))
			continue;
		int e = kh_key(hits->edges, it);
		int rc = g->edges[e].rc_id;
		kh_set_int_add(tmp, e);
		kh_set_int_add(tmp, rc);
	}

	*edges = calloc(kh_size(tmp), sizeof(int));
	*n_e = 0;
	for (khiter_t it = kh_begin(tmp); it != kh_end(tmp); ++it){
		if (!kh_exist(tmp, it))
			continue;
		int e = kh_key(tmp, it);
		(*edges)[(*n_e)++] = e;
	}
	kh_destroy(set_int, tmp);
}

void get_all_shortest_paths_dp(struct asm_graph_t *g, khash_t(long_spath) *spath_info)
{
	khash_t(long_int) *L_pre = kh_init(long_int);
	for (int i = 0; i < g->n_e; ++i){
		int v = i;
		if (g->edges[v].seq_len > MAX_RADIUS)
			continue;
		kh_long_int_set(L_pre, GET_CODE(v, v), g->edges[v].seq_len);

		struct shortest_path_info_t *wrapper = calloc(1,
				sizeof(struct shortest_path_info_t));
		wrapper->len = g->edges[v].seq_len;
		wrapper->trace = -1;

		kh_long_spath_set(spath_info, GET_CODE(v, v), wrapper);
	}
	for (int i = 2; i <= MAX_PATH_LEN; ++i){
		khash_t(long_int) *L_cur = kh_init(long_int);
		for (int j = 0; j < g->n_e; ++j){
			if ((j + 1) % 10000 == 0)
				log_debug("%d-th iteration: %d edges processed",
						i, j + 1);
			int v = j;
			if (g->edges[v].seq_len > MAX_RADIUS)
				continue;
			int *edges;
			int n_e;
			bfs_nearby(g, v, i, &edges, &n_e);
			for (int k = 0; k < n_e; ++k){
				int u = edges[k];
				if (g->edges[u].seq_len > MAX_RADIUS)
					continue;
				uint64_t code = GET_CODE(v, u);
				int tg = g->edges[v].target;
				int min_len = 1e9;
				int next = -1;
				for (int h = 0; h < g->nodes[tg].deg; ++h){
					int w = g->nodes[tg].adj[h];
					uint64_t code = GET_CODE(w, u);
					if (kh_long_int_exist(L_pre, code)
						== 0)
						continue;
					int new_len = kh_long_int_get(L_pre, code)
						+ g->edges[v].seq_len - g->ksize;
					if (min_len > new_len){
						min_len = new_len;
						next = w;
					}
				}
				if (min_len > MAX_RADIUS)
					continue;
				kh_long_int_set(L_cur, code, min_len);

				struct shortest_path_info_t *wrapper;
				if (kh_long_spath_exist(spath_info, code)){
					wrapper = kh_long_spath_get(spath_info, code);
				} else {
					wrapper = calloc(1, sizeof(struct shortest_path_info_t));
					wrapper->len = 1e9;
					wrapper->trace = 0;
					kh_long_spath_set(spath_info, code, wrapper);
				}
				if (wrapper->len > min_len){
					wrapper->len = min_len;
					wrapper->trace = next;
				}
			}
			free(edges);
		}
		kh_destroy(long_int, L_pre);
		L_pre = L_cur;
	}
	kh_destroy(long_int, L_pre);
}

void bfs_nearby(struct asm_graph_t *g, int source, int radius, int **edges, int *n_e)
{
	khash_t(int_int) *L = kh_init(int_int);
	struct queue_t q;
	init_queue(&q, 1024);
	push_queue(&q, pointerize(&source, sizeof(int)));
	kh_int_int_set(L, source, 1);
	while (!is_queue_empty(&q)){
		int v = *(int *) get_queue(&q);
		free(get_queue(&q));
		pop_queue(&q);

		int len = kh_int_int_get(L, v);
		if (len > radius)
			break;
		int tg = g->edges[v].target;
		for (int i = 0; i < g->nodes[tg].deg; ++i){
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
	for (khiter_t it = kh_begin(L); it != kh_end(L); ++it){
		if (!kh_exist(L, it))
			continue;
		(*edges)[(*n_e)++] = kh_key(L, it);
	}
	kh_destroy(int_int, L);
}

