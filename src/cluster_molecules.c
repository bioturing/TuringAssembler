#include "attribute.h"
#include "cluster_molecules.h"
#include "helper.h"
#include "verbose.h"
#include "assembly_graph.h"
#include "minimizers/count_barcodes.h"
#include "minimizers/smart_load.h"
#include "minimizers/minimizers.h"

#define MAX_RADIUS 10000
#define MAX_PATH_LEN 30
#define MIN_BC_READ_COUNT 10
#define MAX_BC_READ_COUNT 88
#define MIN_BARCODE_EDGE_COUNT 100

int cmp_dijkstra(void *node1, void *node2)
{
	struct dijkstra_node_t *n1 = node1;
	struct dijkstra_node_t *n2 = node2;
	if (n1->len < n2->len)
		return -1;
	if (n1->len > n2->len)
		return 1;
	if (n1->n_nodes < n2->n_nodes)
		return -1;
	if (n1->n_nodes > n2->n_nodes)
		return 1;
	return 0;
}

void dijkstra(struct asm_graph_t *g, int source, khash_t(int_int) *distance)
{
	struct heap_t *heap = calloc(1, sizeof(struct heap_t));
	init_heap(heap, &cmp_dijkstra);
	struct dijkstra_node_t wrapper = {
		.vertex = source,
		.len = g->edges[source].seq_len,
		.n_nodes = 0
	};
	push_heap(heap, pointerize(&wrapper, sizeof(struct dijkstra_node_t)));

	put_in_map(distance, source, wrapper.len);
	khash_t(int_int) *n_nodes = kh_init(int_int);
	put_in_map(n_nodes, source, wrapper.n_nodes);

	int ignore_rc = !is_repeat(g, source);
	while (!is_heap_empty(heap)){
		struct dijkstra_node_t *node = get_heap(heap);
		pop_heap(heap);
		int v = node->vertex;
		int len = node->len;
		int path_len = node->n_nodes;
		free(node);
		if (get_in_map(distance, v) != len
			|| get_in_map(n_nodes, v) != path_len)
			continue;
		//printf("%d %d %d %d\n", source, v, len, path_len);
		int tg = g->edges[v].target;
		for (int i = 0; i < g->nodes[tg].deg; ++i){
			int u = g->nodes[tg].adj[i];
			if (u == g->edges[source].rc_id && ignore_rc)
				continue;
			int new_len = len + g->edges[u].seq_len - g->ksize;
			int new_path_len = path_len + 1;
			if (new_len > MAX_RADIUS || new_path_len > MAX_PATH_LEN)
				continue;
			if (check_in_map(distance, u) == 0){
				put_in_map(distance, u, 2e9);
				put_in_map(n_nodes, u, 2e9);
			}
			int cur_len = get_in_map(distance, u);
			int cur_path_len = get_in_map(n_nodes, u);
			if (cur_len > new_len || (cur_len == new_len
					&& cur_path_len > new_path_len)){
				khiter_t it = kh_get(int_int, distance, u);
				kh_val(distance, it) = new_len;
				it = kh_get(int_int, n_nodes, u);
				kh_val(n_nodes, it) = new_path_len;
				wrapper.vertex = u;
				wrapper.len = new_len;
				wrapper.n_nodes = new_path_len;
				push_heap(heap, pointerize(&wrapper,
					sizeof(struct dijkstra_node_t)));
			}
		}
	}
	kh_destroy(int_int, n_nodes);
	heap_destroy(heap);
	free(heap);
}

int get_shortest_path(struct asm_graph_t *g, int source, int target)
{
	khash_t(int_int) *L = kh_init(int_int);
	dijkstra(g, source, L);
	int res = -1;
	if (check_in_map(L, target))
		res = get_in_map(L, target);
	kh_destroy(int_int, L);
	return res;
}

void get_all_shortest_paths(struct asm_graph_t *g, khash_t(long_int) *distance)
{
	for (int i = 0; i < g->n_e; ++i){
		if ((i + 1) % 1000 == 0 || i + 1 == g->n_e)
			log_debug("%d/%d edges processed", i + 1, g->n_e);
		if (g->edges[i].seq_len > MAX_RADIUS)
			continue;
		khash_t(int_int) *D = kh_init(int_int);
		dijkstra(g, i, D);
		for (khiter_t it = kh_begin(D); it != kh_end(D); ++it){
			if (!kh_exist(D, it))
				continue;
			int v = i;
			int u = kh_key(D, it);
			int val = kh_val(D, it);
			uint64_t code = (((uint64_t) v) << 32) | u;
			int ret;
			khiter_t it = kh_put(long_int, distance, code, &ret);
			kh_val(distance, it) = val;
		}
		kh_destroy(int_int, D);
	}
	log_debug("Done finding all shortest paths");
}

int get_pair_distance(int v, int u, khash_t(long_int) *distance)
{
	uint64_t code = (((uint64_t) v) << 32) | u;
	khiter_t it = kh_get(long_int, distance, code);
	if (it == kh_end(distance))
		return -1;
	return kh_val(distance, it);
}

int check_connected(struct asm_graph_t *g, int v, int u,
		khash_t(long_int) *distance)
{
	int tg = g->edges[v].target;
	int sr = g->nodes[g->edges[u].source].rc_id;
	for (int i = 0; i < g->nodes[tg].deg; ++i){
		for (int j = 0; j < g->nodes[sr].deg; ++j){
			int w = g->nodes[tg].adj[i];
			int t = g->edges[g->nodes[sr].adj[j]].rc_id;
			if (w == u || t == v)
				return 1;
			int d = get_pair_distance(w, t, distance);
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
		khash_t(long_int) *distance, khash_t(long_int) *is_connected,
		khash_t(long_int) *count_link)
{
	for (int i = 0; i < n_e; ++i){
		for (int j = 0; j < n_e; ++j){
			if (i == j)
				continue;
			int v = edges[i];
			int u = edges[j];
			uint64_t code = (((uint64_t) v) << 32) | u;
			int ok;
			khiter_t it = kh_get(long_int, is_connected, code);
			if (it == kh_end(is_connected)){
				ok = check_connected(g, v, u, distance);
				int ret;
				it = kh_put(long_int, is_connected, code, &ret);
				kh_val(is_connected, it) = ok;
			} else {
				ok = kh_val(is_connected, it);
			}
			if (!ok)
				continue;
			it = kh_get(long_int, count_link, code);
			if (it == kh_end(count_link)){
				int ret;
				it = kh_put(long_int, count_link, code, &ret);
				kh_val(count_link, it) = 0;
			}
			++kh_val(count_link, it);
		}
	}
}

void count_edge_links_bc(struct opt_proc_t *opt)
{
	khash_t(long_int) *link_count = kh_init(long_int);

	struct bc_hit_bundle_t bc_hit_bundle;
	get_bc_hit_bundle(opt, &bc_hit_bundle);
	struct asm_graph_t *g = bc_hit_bundle.g;
	khash_t(bcpos) *bc_pos_dict = bc_hit_bundle.bc_pos_dict;
	struct read_path_t *read_sorted_path = bc_hit_bundle.read_sorted_path;

	khash_t(long_int) *distance = kh_init(long_int);
	khash_t(long_int) *is_connected = kh_init(long_int);
	get_all_shortest_paths(bc_hit_bundle.g, distance);

	struct barcode_list_t blist;
	get_barcode_list(opt->bx_str, &blist);

	for (int i = 0; i < blist.n_bc; ++i){
		if ((i + 1) % 10000 == 0)
			log_debug("%d/%d barcodes processed", i + 1, blist.n_bc);
		if (blist.read_count[i] < MIN_BC_READ_COUNT
			|| blist.read_count[i] > MAX_BC_READ_COUNT)
			continue;
		uint64_t bx_encoded = barcode_hash_mini(blist.bc_list[i]);
		uint64_t bx[1] = {bx_encoded};

		khint_t k = kh_get(bcpos, bc_pos_dict, bx_encoded);
		if (k == kh_end(bc_pos_dict)) {
			log_error("Barcode does not exist");
		}

		struct mm_hits_t *hits = get_hits_from_barcode(blist.bc_list[i],
				&bc_hit_bundle);

		int *edges = calloc(kh_size(hits->edges), sizeof(int));
		int n_e = 0;
		for (khiter_t it = kh_begin(hits->edges); it != kh_end(hits->edges);
				++it){
			if (!kh_exist(hits->edges, it))
				continue;
			edges[n_e++] = kh_key(hits->edges, it);
		}
		mm_hits_destroy(hits);

		get_edge_links_by_distance(bc_hit_bundle.g, edges, n_e, distance,
				is_connected, link_count);
		free(edges);
	}

	barcode_list_destroy(&blist);
	bc_hit_bundle_destroy(&bc_hit_bundle);
	kh_destroy(long_int, distance);
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

void get_sub_graph(struct asm_graph_t *g, struct mm_hits_t *hits,
		khash_t(long_int) *pair_count)
{
	khash_t(set_int) *edges = kh_init(set_int);
	for (khiter_t it = kh_begin(hits->edges); it != kh_end(hits->edges); ++it){
		if (!kh_exist(hits->edges, it))
			continue;
		int e = kh_key(hits->edges, it);
		put_in_set(edges, e);
		put_in_set(edges, g->edges[e].rc_id);
	}

	int *tmp = calloc(kh_size(edges), sizeof(int));
	int n = 0;
	for (khiter_t it = kh_begin(edges); it != kh_end(edges); ++it){
		if (!kh_exist(edges, it))
			continue;
		tmp[n++] = kh_key(edges, it);
	}
	kh_destroy(set_int, edges);

	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j){
			if (i == j)
				continue;
			int v = tmp[i];
			int u = tmp[j];

			int tg = g->edges[v].target;
			int sr = g->nodes[g->edges[u].source].rc_id;
			int ok = 0;
			for (int h = 0; !ok && h < g->nodes[tg].deg; ++h){
				for (int k = 0; !ok && k < g->nodes[sr].deg; ++k){
					int w = g->nodes[tg].adj[h];
					int t = g->edges[g->nodes[sr].adj[k]].rc_id;
					if (w == u || t == v){
						ok = 1;
						break;
					}
					int len = get_shortest_path(g, w, t);
					if (len == -1)
						continue;
					if(len > MAX_RADIUS)
						log_error("Something went wrong, probably Dijkstra is incorrect");
					ok = 1;
				}
			}

			if (!ok)
				continue;
			uint64_t code = (((uint64_t) v) << 32) | u;
			int ret;
			kh_put(long_int, pair_count, code, &ret);
		}
	}
	free(tmp);
}

void print_barcode_graph(struct opt_proc_t *opt)
{
	khash_t(long_int) *h_all = kh_init(long_int);
	FILE *f = fopen(opt->in_fasta, "r");
	int u, v, c;
	while (fscanf(f, "%d %d %d\n", &u, &v, &c) == 3){
		uint64_t code = (((uint64_t) u) << 32) | v;
		khiter_t it = kh_get(long_int, h_all, code);
		if (it != kh_end(h_all))
			log_error("Something went wrong");
		int ret;
		it = kh_put(long_int, h_all, code, &ret);
		kh_val(h_all, it) = c;
	}
	fclose(f);

	khash_t(long_int) *h_1 = kh_init(long_int);
	struct asm_graph_t g;
	load_asm_graph(&g, opt->in_file);
	struct bc_hit_bundle_t bc_hit_bundle;
	get_bc_hit_bundle(opt, &bc_hit_bundle);
	struct mm_hits_t *hits = get_hits_from_barcode(opt->bx_str, &bc_hit_bundle);
	get_sub_graph(&g, hits, h_1);

	f = fopen(opt->lc, "w");
	fprintf(f, "digraph %s{\n", opt->bx_str);
	float unit_cov = get_genome_coverage(&g);
	khash_t(set_int) *nodes = kh_init(set_int);
	for (khiter_t it = kh_begin(h_1); it != kh_end(h_1); ++it){
		if (!kh_exist(h_1, it))
			continue;
		uint64_t code = kh_key(h_1, it);
		khiter_t it2 = kh_get(long_int, h_all, code);
		if (it2 != kh_end(h_all)){
			int val = kh_val(h_all, it2);
			if (val < 300)
				continue;
			int u = code >> 32;
			int v = code & ((((uint64_t) 1) << 32) - 1);
			put_in_set(nodes, u);
			put_in_set(nodes, v);
			fprintf(f, "\t%d -> %d [label=\"%d\"]\n", u, v, val);
		}
	}

	for (khiter_t it = kh_begin(nodes); it != kh_end(nodes); ++it){
		if (!kh_exist(nodes, it))
			continue;
		int u = kh_key(nodes, it);
		float cov = __get_edge_cov(g.edges + u, g.ksize);
		float ratio = cov / unit_cov;
		if (ratio >= 0.8 && ratio <= 1.2)
			fprintf(f, "%d [style=\"filled\",fillcolor=green]\n",
					u);
		else
			fprintf(f, "%d [style=\"filled\",fillcolor=violet]\n",
					u);
	}
	kh_destroy(set_int, nodes);
	fprintf(f, "}");
	fclose(f);
}

void get_barcode_edges_path(struct opt_proc_t *opt)
{
	struct bc_hit_bundle_t *bc_hit_bundle = calloc(1,
			sizeof(struct bc_hit_bundle_t));
	get_bc_hit_bundle(opt, bc_hit_bundle);
	struct barcode_list_t blist;
	get_barcode_list(opt->bx_str, &blist);

	khash_t(long_int) *all_pairs = kh_init(long_int);
	get_all_pair_edge_count(opt->in_fasta, all_pairs); // Option -f

	FILE *f = fopen(opt->lc, "w");
	for (int i = 0; i < blist.n_bc; ++i){
		if ((i + 1) % 10000 == 0)
			log_debug("Processing %d-th barcode", i + 1);
		struct mm_hits_t *hits = get_hits_from_barcode(blist.bc_list[i],
				bc_hit_bundle);
		khash_t(long_int) *pair_count = kh_init(long_int);
		get_sub_graph(bc_hit_bundle->g, hits, pair_count);

		struct simple_graph_t sg;
		init_simple_graph(&sg);
		build_simple_graph(pair_count, all_pairs, &sg);
		find_DAG(&sg, bc_hit_bundle->g);
		get_longest_path(&sg);

		khash_t(set_int) *not_source = kh_init(set_int);
		for (khiter_t it = kh_begin(sg.nodes); it != kh_end(sg.nodes);
				++it){
			if (!kh_exist(sg.nodes, it))
				continue;
			int u = kh_key(sg.nodes, it);
			struct simple_node_t *snode = kh_val(sg.nodes, it);
			for (int j = 0; j < snode->deg; ++j)
				put_in_set(not_source, snode->adj[j]);
		}

		for (khiter_t it = kh_begin(sg.nodes); it != kh_end(sg.nodes);
				++it){
			if (!kh_exist(sg.nodes, it))
				continue;
			int u = kh_key(sg.nodes, it);
			if (check_in_set(not_source, u) || check_in_set(sg.is_loop, u))
				continue;
			fprintf(f, "%s: ", blist.bc_list[i]);
			fprintf(f, "%d", u);
			for (int v = get_in_map(sg.next, u); v != -1;
					v = get_in_map(sg.next, v))
				fprintf(f, " --> %d", v);
			fprintf(f, "\n");
		}
		fflush(f);
		kh_destroy(set_int, not_source);


		/*for (khiter_t it = kh_begin(sg.nodes); it != kh_end(sg.nodes);
				++it){
			if (!kh_exist(sg.nodes, it))
				continue;
			int u = kh_key(sg.nodes, it);
			struct simple_node_t *snode = kh_val(sg.nodes, it);
			for (int j = 0; j < snode->deg; ++j)
				fprintf(f, "\t%d -> %d\n", u, snode->adj[j]);
		}*/

		simple_graph_destroy(&sg);
		kh_destroy(long_int, pair_count);

	}
	fclose(f);

	kh_destroy(long_int, all_pairs);
	bc_hit_bundle_destroy(bc_hit_bundle);
	free(bc_hit_bundle);
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

void get_all_pair_edge_count(char *file_path, khash_t(long_int) *pair_count)
{
	FILE *f = fopen(file_path, "r");
	int u, v, count;
	while (fscanf(f, "%d %d %d\n", &u, &v, &count) == 3){
		uint64_t code = (((uint64_t) u) << 32) | v;
		khiter_t it = kh_get(long_int, pair_count, code);
		if (it != kh_end(pair_count))
			log_error("Key is already in set, something went wrong");
		int ret;
		it = kh_put(long_int, pair_count, code, &ret);
		kh_val(pair_count, it) = count;
	}
	fclose(f);
}

void add_simple_node(struct simple_graph_t *sg, int u)
{
	khiter_t it = kh_get(int_node, sg->nodes, u);
	if (it == kh_end(sg->nodes)){
		int ret;
		struct simple_node_t *snode = calloc(1,
				sizeof(struct simple_node_t));
		it = kh_put(int_node, sg->nodes, u, &ret);
		kh_val(sg->nodes, it) = snode;
	}
}

void add_simple_edge(struct simple_graph_t *sg, int u, int v)
{
	khiter_t it = kh_get(int_node, sg->nodes, u);
	if (it == kh_end(sg->nodes))
		log_error("Key not in hash");
	struct simple_node_t *snode = kh_val(sg->nodes, it);
	snode->adj = realloc(snode->adj, sizeof(int) * (snode->deg + 1));
	snode->adj[snode->deg] = v;
	++snode->deg;
}

void init_simple_graph(struct simple_graph_t *sg)
{
	sg->nodes = kh_init(int_node);
	sg->is_loop = kh_init(set_int);
	sg->path_len = kh_init(int_int);
	sg->next = kh_init(int_int);
}

void build_simple_graph(khash_t(long_int) *one_bc, khash_t(long_int) *all_bc,
		struct simple_graph_t *sg)
{
	for (khiter_t it = kh_begin(one_bc); it != kh_end(one_bc); ++it){
		if (!kh_exist(one_bc, it))
			continue;
		if (kh_val(one_bc, it) == -1)
			continue;
		uint64_t code = kh_key(one_bc, it);
		khiter_t it2 = kh_get(long_int, all_bc, code);
		if (it2 == kh_end(all_bc))
			continue;
		int val = kh_val(all_bc, it2);
		if (val < MIN_BARCODE_EDGE_COUNT)
			continue;
		int u = code >> 32;
		int v = code & ((uint32_t) -1);
		add_simple_node(sg, u);
		add_simple_node(sg, v);
		add_simple_edge(sg, u, v);
	}
}

void simple_graph_destroy(struct simple_graph_t *sg)
{
	for (khiter_t it = kh_begin(sg->nodes); it != kh_end(sg->nodes); ++it){
		if (!kh_exist(sg->nodes, it))
			continue;
		struct simple_node_t *snode = kh_val(sg->nodes, it);
		free(snode->adj);
		free(snode);
	}
	kh_destroy(int_node, sg->nodes);
}

void check_loop_dfs(struct simple_graph_t *sg, int u, khash_t(set_int) *visited,
		khash_t(set_int) *in_dfs)
{
	if (check_in_set(in_dfs, u)){
		put_in_set(sg->is_loop, u);
		return;
	}
	if (check_in_set(visited, u))
		return;
	put_in_set(visited, u);
	put_in_set(in_dfs, u);
	khiter_t it = kh_get(int_node, sg->nodes, u);
	struct simple_node_t *snode = kh_val(sg->nodes, it);
	for (int i = 0; i < snode->deg; ++i){
		int v = snode->adj[i];
		check_loop_dfs(sg, v, visited, in_dfs);
	}
	erase_from_set(in_dfs, u);
}

void find_DAG(struct simple_graph_t *sg, struct asm_graph_t *g)
{
	khash_t(int_node) *nodes = sg->nodes;
	khash_t(set_int) *visited = kh_init(set_int);
	khash_t(set_int) *in_dfs = kh_init(set_int);
	for (khiter_t it = kh_begin(nodes); it != kh_end(nodes); ++it){
		if (!kh_exist(nodes, it))
			continue;
		int u = kh_key(nodes, it);
		check_loop_dfs(sg, u, visited, in_dfs);
	}

	struct queue_t q;
	init_queue(&q, 1024);
	kh_destroy(set_int, visited);
	visited = kh_init(set_int);
	for (khiter_t it = kh_begin(sg->is_loop); it != kh_end(sg->is_loop); ++it){
		if (!kh_exist(sg->is_loop, it))
			continue;
		int u = kh_key(sg->is_loop, it);
		int u_rc = g->edges[u].rc_id;
		push_queue(&q, pointerize(&u, sizeof(int)));
		if (check_in_set(sg->is_loop, u_rc) == 0)
			push_queue(&q, pointerize(&u_rc, sizeof(int)));
		put_in_set(visited, u);
		put_in_set(visited, u_rc);
	}

	while (!is_queue_empty(&q)){
		int u = *(int *) get_queue(&q);
		pop_queue(&q);
		khiter_t it = kh_get(int_node, sg->nodes, u);
		struct simple_node_t *snode = kh_val(sg->nodes, it);
		for (int i = 0; i < snode->deg; ++i){
			int v = snode->adj[i];
			if (check_in_set(visited, v) == 0){
				put_in_set(visited, v);
				push_queue(&q, pointerize(&v, sizeof(int)));
			}
		}
	}
	for (khiter_t it = kh_begin(visited); it != kh_end(visited); ++it){
		if (!kh_exist(visited, it))
			continue;
		int u = kh_key(visited, it);
		int u_rc = g->edges[u].rc_id;
		put_in_set(sg->is_loop, u);
		put_in_set(sg->is_loop, u_rc);
	}


	/*for (khiter_t it = kh_begin(sg->is_loop); it != kh_end(sg->is_loop); ++it){
		if (!kh_exist(sg->is_loop, it))
			continue;
		__VERBOSE("%d\n", kh_key(sg->is_loop, it));
	}*/

	kh_destroy(set_int, in_dfs);
	kh_destroy(set_int, visited);
}

void get_longest_path_dfs(struct simple_graph_t *sg, int u,
		khash_t(set_int) *done_dfs)
{
	if (check_in_set(done_dfs, u))
		return;
	khiter_t it = kh_get(int_node, sg->nodes, u);
	struct simple_node_t *snode = kh_val(sg->nodes, it);
	int max_len = 0;
	int next = -1;
	for (int i = 0; i < snode->deg; ++i){
		int v = snode->adj[i];
		get_longest_path_dfs(sg, v, done_dfs);
		int next_len = get_in_map(sg->path_len, v);
		if (max_len < next_len){
			max_len = next_len;
			next = v;
		}
	}
	put_in_map(sg->path_len, u, max_len + 1);
	put_in_map(sg->next, u, next);
	put_in_set(done_dfs, u);
}

void get_longest_path(struct simple_graph_t *sg)
{
	khash_t(set_int) *done_dfs = kh_init(set_int);
	for (khiter_t it = kh_begin(sg->nodes); it != kh_end(sg->nodes); ++it){
		if (!kh_exist(sg->nodes, it))
			continue;
		int u = kh_key(sg->nodes, it);
		if (check_in_set(sg->is_loop, u))
			continue;
		get_longest_path_dfs(sg, u, done_dfs);
	}
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
