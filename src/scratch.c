
void main_process(struct opt_count_t *opt)
{
	struct kmhash_t *V;
	__VERBOSE("Counting kmer...\n");
	V = count_kmer(opt);
	khash_t(kvert) *hvert;
	__VERBOSE("Filtering vertices...\n");
	hvert = filter_kmer(V, opt);
	kmhash_destroy(V);

	uint32_t *edge_count;
	__VERBOSE("Counting edges...\n");
	edge_count = get_edges(opt, hvert);

	// __VERBOSE("Dumping raw graph...\n");
	// dump_graph(opt, hvert, edge_count);

	__VERBOSE("Dumping reduced graph...\n");
	reduce_graph(opt, hvert, edge_count);
}

khash_t(kvert) *filter_kmer(struct kmhash_t *V, struct opt_count_t *opt)
{
	// can be done in parallel
	// char dump_path[1024];
	// strcpy(dump_path, opt->out_dir);
	// strcat(dump_path, "/dump.tsv");
	// FILE *fp = xfopen(dump_path, "wb");
	kmint_t i;
	kmkey_t tombstone;
	tombstone = (kmkey_t)-1;
	khash_t(kvert) *h = kh_init(kvert);
	int n_chosen;
	n_chosen = 0;
	khiter_t k;
	int ret;
	for (i = 0; i < V->size; ++i) {
		if (V->bucks[i].idx == tombstone)
			continue;
		if (V->bucks[i].cnt > (uint32_t)opt->filter_thres) {
			k = kh_put(kvert, h, V->bucks[i].idx, &ret);
			kh_value(h, k).cnt = V->bucks[i].cnt;
			kh_value(h, k).idx = n_chosen++;
			// fprintf(fp, "%llu\t%u\n", (unsigned long long)V->bucks[i].idx, (unsigned int)V->bucks[i].cnt);
		}
	}
	__VERBOSE_LOG("Result", "Number of filtered vertices        : %20d\n", (int)V->n_items - n_chosen);
	// fclose(fp);
	return h;
}

static void dump_seq(uint64_t num, char *seq, int len)
{
	seq[len] = '\0';
	while (len) {
		seq[--len] = nt4_char[num & 3];
		num >>= 2;
	}
}

#define __get_revc_num(y, x, l, mask)					       \
(	(x) = (y) << (64 - ((l) << 1)),					       \
	(x) = (((x) & 0xffffffff00000000ull) >> 32) | (((x) & 0x00000000ffffffffull) << 32), \
	(x) = (((x) & 0xffff0000ffff0000ull) >> 16) | (((x) & 0x0000ffff0000ffffull) << 16), \
	(x) = (((x) & 0xff00ff00ff00ff00ull) >>  8) | (((x) & 0x00ff00ff00ff00ffull) <<  8), \
	(x) = (((x) & 0xf0f0f0f0f0f0f0f0ull) >>  4) | (((x) & 0x0f0f0f0f0f0f0f0full) <<  4), \
	(x) = (((x) & 0xccccccccccccccccull) >>  2) | (((x) & 0x3333333333333333ull) <<  2), \
	(x) ^= 0xffffffffffffffffull, (x) &= (mask))

void dump_graph(struct opt_count_t *opt, khash_t(kvert) *h, int16_t *edges)
{
	char dump_path[1024];
	strcpy(dump_path, opt->out_dir);
	strcat(dump_path, "/graph_raw.gfa");
	FILE *fp = xfopen(dump_path, "w");
	int ksize, node_id, c;
	uint64_t kmask, idx, rev_idx, adj_idx, rev_adj_idx;
	int16_t *adj;
	ksize = opt->kmer_size;
	kmask = (1ull << (ksize << 1)) - 1;
	char *seq = malloc(ksize + 1);

	khint_t i, ik;
	for (i = kh_begin(h); i != kh_end(h); ++i) {
		if (!kh_exist(h, i))
			continue;
		dump_seq(kh_key(h, i), seq, ksize);
		fprintf(fp, "S\t%d\t%s\tKC:i:%d\n", kh_value(h, i).idx, seq, kh_value(h, i).cnt);
	}

	for (i = kh_begin(h); i != kh_end(h); ++i) {
		if (!kh_exist(h, i))
			continue;
		idx = kh_key(h, i);
		node_id = kh_value(h, i).idx;
		__get_revc_num(idx, rev_idx, ksize, kmask);

		adj = edges + (node_id * 8);
		for (c = 0; c < 4; ++c) {
			if (adj[c]) {
				adj_idx = ((idx << 2) & kmask) | c;
				rev_adj_idx = (rev_idx >> 2) | ((uint64_t)(c ^ 3) << ((ksize << 1) - 2));
				if (adj_idx < rev_adj_idx) {
					ik = kh_get(kvert, h, adj_idx);
					assert(ik != kh_end(h));
					fprintf(fp, "L\t%d\t+\t%d\t+\t%dM\n",
						node_id, kh_value(h, ik).idx, ksize - 1);
				} else {
					ik = kh_get(kvert, h, rev_adj_idx);
					assert(ik != kh_end(h));
					fprintf(fp, "L\t%d\t+\t%d\t-\t%dM\n",
						node_id, kh_value(h, ik).idx, ksize - 1);
				}
			}
			if (adj[c + 4]) {
				adj_idx = ((rev_idx << 2) & kmask) | c;
				rev_adj_idx = (idx >> 2) | ((uint64_t)(c ^ 3) << ((ksize << 1) - 2));
				if (adj_idx < rev_adj_idx) {
					ik = kh_get(kvert, h, adj_idx);
					assert(ik != kh_end(h));
					fprintf(fp, "L\t%d\t-\t%d\t+\t%dM\n",
						node_id, kh_value(h, ik).idx, ksize - 1);
				} else {
					ik = kh_get(kvert, h, rev_adj_idx);
					assert(ik != kh_end(h));
					fprintf(fp, "L\t%d\t-\t%d\t-\t%dM\n",
						node_id, kh_value(h, ik).idx, ksize - 1);
				}
			}
		}
	}

	free(seq);
	fclose(fp);
}

