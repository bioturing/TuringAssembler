#include "attribute.h"
#include "cluster_molecules.h"
#include "helper.h"
#include "verbose.h"

#include "minimizers/count_barcodes.h"
#include "minimizers/smart_load.h"
#include "minimizers/minimizers.h"

#define MAX_RADIUS 7000
#define MAX_PATH_LEN 10
#define MIN_BC_READ_COUNT 10
#define MAX_BC_READ_COUNT 88
#define MIN_BARCODE_EDGE_COUNT 100

int get_shortest_path(struct asm_graph_t *g, int source, int target)
{
	struct queue_t *q = calloc(1, sizeof(struct queue_t));
	struct dijkstra_node_t wrapper = {
		.vertex = source,
		.len = 0
	};
	push_queue(q, pointerize(&wrapper, sizeof(struct dijkstra_node_t)));

	khash_t(int_int) *L = kh_init(int_int);
	khash_t(int_int) *n_nodes = kh_init(int_int);
	put_in_map(L, source, wrapper.len);
	put_in_map(n_nodes, source, 0);
	while (!is_queue_empty(q)){
		int p = q->front;
		for (int i = q->front + 1; i < q->back; ++i){
			if (((struct dijkstra_node_t *)q->data[p])->len >
				((struct dijkstra_node_t *)q->data[i])->len)
				p = i;
		}
		if (p != q->front){
			struct dijkstra_node_t *tmp = q->data[q->front];
			q->data[q->front] = q->data[p];
			q->data[p] = tmp;
		}
		struct dijkstra_node_t *node = get_queue(q);
		pop_queue(q);
		int v = node->vertex;
		int len = node->len;
		free(node);
		int path_len = get_in_map(n_nodes, v);
		if (get_in_map(L, v) != len)
			continue;
		if (v == target)
			break;
		if (len > MAX_RADIUS)
			break;
		if (path_len > MAX_PATH_LEN)
			continue;
		int tg = g->edges[v].target;
		for (int i = 0; i < g->nodes[tg].deg; ++i){
			int u = g->nodes[tg].adj[i];
			if (check_in_map(L, u) == 0){
				put_in_map(L, u, 2e9);
				put_in_map(n_nodes, u, path_len + 1);
			}
			if ((uint32_t) get_in_map(L, u) > len + g->edges[u].seq_len){
				khiter_t it = kh_get(int_int, L, u);
				kh_val(L, it) = len + g->edges[u].seq_len;
				wrapper.vertex = u;
				wrapper.len = len + g->edges[u].seq_len;
				push_queue(q, pointerize(&wrapper,
					sizeof(struct dijkstra_node_t)));
				it = kh_get(int_int, n_nodes, u);
				kh_val(n_nodes, it) = path_len + 1;
			}
		}
	}
	int res = -1;
	if (check_in_map(L, target) != 0)
		res = get_in_map(L, target) - g->edges[target].seq_len;
	free_queue_content(q);
	destroy_queue(q);
	free(q);

	/*if (res != -1){
		__VERBOSE("PATH FROM %d to %d:\n", source, target);
		__VERBOSE("DISTANCE %d %d\n", get_in_map(L, target),
				get_in_map(n_nodes, target));
	}*/
	kh_destroy(int_int, L);
	kh_destroy(int_int, n_nodes);
	return res;
}

void count_edge_links_bc(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct read_path_t *read_sorted_path, khash_t(bcpos) *bx_pos_dict,
		struct mm_db_edge_t *mm_edges, char **bc_list, int n_bc)
{
	khash_t(long_int) *pair_count = kh_init(long_int);
	for (int i = 0; i < n_bc; ++i){
		if (i % 10000 == 0)
			log_debug("%d barcodes processed", i + 1);
		//log_debug("Barcode: %s", bc_list[i]);
		uint64_t bx_encoded = barcode_hash_mini(bc_list[i]);
		uint64_t bx[1] = {bx_encoded}; //43 15 mock barcode pseudo hash id here

		khint_t k = kh_get(bcpos, bx_pos_dict, bx_encoded);          // query the hash table
		if (k == kh_end(bx_pos_dict)) {
			log_error("Barcode does not exist");
		}

		char *buf1, *buf2;
		uint64_t m_buf1, m_buf2;
		stream_filter_read(read_sorted_path, bx_pos_dict, bx, 1,
				&buf1, &buf2, &m_buf1, &m_buf2);

		struct read_t r1, r2;
		int pos1 = 0, pos2 = 0;
		int n_reads = 0;
		struct mm_hits_t *hits;
		hits = mm_hits_init();


		while (get_read_from_fq(&r1, buf1, &pos1) == READ_SUCCESS
			&& get_read_from_fq(&r2, buf2, &pos2) == READ_SUCCESS ) {
			n_reads++;
			struct mm_db_t *db1, *db2;
			db1 = mm_index_char_str(r1.seq, MINIMIZERS_KMER, MINIMIZERS_WINDOW, r1.len);
			db2 = mm_index_char_str(r2.seq, MINIMIZERS_KMER, MINIMIZERS_WINDOW, r2.len);

			mm_hits_cmp(db1, mm_edges, hits, g);
			mm_hits_cmp(db2, mm_edges, hits, g);
			mm_db_destroy(db1);
			mm_db_destroy(db2);
		}
		
		get_sub_graph(g, hits, pair_count);
		mm_hits_destroy(hits);
		free(buf1);
		free(buf2);
	}
	FILE *f = fopen(opt->lc, "w");
	for (khiter_t it = kh_begin(pair_count); it != kh_end(pair_count); ++it){
		if (!kh_exist(pair_count, it))
			continue;
		uint64_t code = kh_key(pair_count, it);
		int count = kh_val(pair_count, it);
		int u = code >> 32;
		int v = code & ((((uint64_t) 1) << 32) - 1);
		if (count != -1)
			fprintf(f, "%d %d %d\n", u, v, count);
	}
	fclose(f);
	kh_destroy(long_int, pair_count);
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
			uint64_t code = (((uint64_t) tmp[i])) << 32 | tmp[j];
			khiter_t it = kh_get(long_int, pair_count, code);
			if (it == kh_end(pair_count)){
				int len = get_shortest_path(g, tmp[i], tmp[j]);
				int ret;
				it = kh_put(long_int, pair_count, code, &ret);
				if (len != -1 && len <= MAX_RADIUS)
					kh_val(pair_count, it) = 1;
				else
					kh_val(pair_count, it) = -1;
			} else {
				if (kh_val(pair_count, it) != -1)
					++kh_val(pair_count, it);
			}
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
	for (khiter_t it = kh_begin(h_1); it != kh_end(h_1); ++it){
		if (!kh_exist(h_1, it))
			continue;
		uint64_t code = kh_key(h_1, it);
		khiter_t it2 = kh_get(long_int, h_all, code);
		if (it2 != kh_end(h_all)){
			int val = kh_val(h_all, it2);
			if (val < 100)
				continue;
			int u = code >> 32;
			int v = code & ((((uint64_t) 1) << 32) - 1);
			fprintf(f, "\t%d -> %d [label=\"%d\"]\n", u, v, val);
		}
	}
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
	char bc[19];
	int read_count;
	FILE *f = fopen(bc_count_path, "r");
	while (fscanf(f, "%s\t%d\n", bc, &read_count) == 2){
		if (read_count < 10 || read_count > 100)
			continue;
		if (n == m){
			m <<= 1;
			bc_list = realloc(bc_list, sizeof(char *) * m);
		}
		bc_list[n] = calloc(19, sizeof(char));
		memcpy(bc_list[n], bc, sizeof(char) * 19);
		++n;
	}
	fclose(f);
	bc_list = realloc(bc_list, sizeof(char *) * n);
	blist->n_bc = n;
	blist->bc_list = bc_list;
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

