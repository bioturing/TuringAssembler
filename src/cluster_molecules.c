#include "cluster_molecules.h"
#include "helper.h"
#include "verbose.h"
#include "minimizers/count_barcodes.h"
#include "minimizers/smart_load.h"
#include "helper.h"

#define MAX_RADIUS 7000
#define MAX_PATH_LEN 10

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
		int path_len = get_in_map(n_nodes, v);
		if (get_in_map(L, v) != len)
			goto dijkstra_node_outdated;
		if (v == target)
			break;
		if (len > MAX_RADIUS || path_len > MAX_PATH_LEN)
			break;
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
dijkstra_node_outdated:
		free(node);
	}
	int res = -1;
	if (check_in_map(L, target) != 0)
		res = get_in_map(L, target) - g->edges[target].seq_len;
	free_queue_content(q);
	destroy_queue(q);

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
		log_debug("Barcode: %s", bc_list[i]);
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
		struct mm_db_t *db1, *db2;
		struct mm_hits_t *hits;
		hits = mm_hits_init();


		while (get_read_from_fq(&r1, buf1, &pos1) == READ_SUCCESS
			&& get_read_from_fq(&r2, buf2, &pos2) == READ_SUCCESS ) {
			n_reads++;
			db1 = mm_index_char_str(r1.seq, 17, 17, r1.len);
			db2 = mm_index_char_str(r2.seq, 17, 17, r2.len);

			mm_hits_cmp(db1, mm_edges, hits, g);
			mm_hits_cmp(db2, mm_edges, hits, g);
		}
		
		get_sub_graph(opt, g, hits, pair_count);
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
		fprintf(f, "%d %d %d\n", u, v, count);
	}
	fclose(f);
	kh_destroy(long_int, pair_count);
}

void get_sub_graph(struct opt_proc_t *opt, struct asm_graph_t *g,
		struct mm_hits_t *hits, khash_t(long_int) *pair_count)
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

	/*FILE *f = fopen(opt->lc, "w");
	fprintf(f, "digraph %s{\n", opt->bx_str);*/
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j){
			if (i == j)
				continue;
			int len = get_shortest_path(g, tmp[i], tmp[j]);
			if (len != -1 && len <= MAX_RADIUS){
				/*fprintf(f, "\t%d -> %d [label=\"%d\"];\n", tmp[i], tmp[j],
					len);*/
				//__VERBOSE("%d %d %d\n", tmp[i], tmp[j], len);
				uint64_t code = (((uint64_t) tmp[i])) << 32 | tmp[j];
				khiter_t it = kh_get(long_int, pair_count, code);
				if (it == kh_end(pair_count)){
					int ret;
					it = kh_put(long_int, pair_count, code,
							&ret);
					kh_val(pair_count, it) = 0;
				}
				++kh_val(pair_count, it);
			}
		}
	}
	//fprintf(f, "}\n");
	//fclose(f);
	free(tmp);
}

