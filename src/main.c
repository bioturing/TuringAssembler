#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "get_buffer.h"
#include "graph.h"
#include "io_utils.h"
#include "khash.h"
#include "kmer_count.h"
#include "kmhash.h"
#include "utils.h"
#include "verbose.h"

__KHASH_IMPL(kvert, kh_inline, kmkey_t, struct kvert_info_t, 1, kh_int64_hash_func, kh_int64_hash_equal)

struct opt_count_t *init_opt_count()
{
	struct opt_count_t *opt;
	opt = calloc(1, sizeof(struct opt_count_t));
	opt->n_threads = 1;
	opt->hash_size = (1 << 24);
	opt->kmer_size = 29;
	opt->n_files = 0;
	opt->filter_thres = 0;
	opt->files_1 = opt->files_2 = NULL;
	opt->out_dir = ".";
	return opt;
}

int opt_count_list(int argc, char *argv[])
{
	int n;
	for (n = 0; n < argc - 1; ++n) {
		if (argv[n + 1][0] == '-')
			break;
	}
	if (n == 0)
		__ERROR("Emtpy list %s", argv[0]);
	return n;
}

struct opt_count_t *parse_count_option(int argc, char *argv[])
{
	int pos = 0, n;
	struct opt_count_t *opt = init_opt_count();
	while (pos < argc) {
		if (!strcmp(argv[pos], "-t")) {
			opt->n_threads = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-s")) {
			opt->hash_size = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-k")) {
			opt->kmer_size = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-o")) {
			opt->out_dir = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-1")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__ERROR("Inconsistent number of files");
			opt->n_files = n;
			opt->files_1 = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-2")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__ERROR("Inconsistent number of files");
			opt->n_files = n;
			opt->files_2 = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "--filter-threshold")) {
			opt->filter_thres = atoi(argv[pos + 1]);
			pos += 2;
		} else if (argv[pos][0] != '-') {
			if (opt->n_files != 0)
				__ERROR("Unknown %s", argv[pos]);
			opt->files_1 = argv + pos;
			while (pos < argc && argv[pos][0] != '-') {
				++pos;
				++opt->n_files;
			}
		} else {
			__ERROR("Unknown option %s", argv[pos]);
		}
	}
	mkdir(opt->out_dir, 0755);
	return opt;
}

void print_usage(const char *prog)
{
	__VERBOSE("Usage: %s [options] -1 read_1.fq -2 read_2.fq\n", prog);
	__VERBOSE("Options: -t                     <number of threads>\n");
	__VERBOSE("         -s                     <pre-alloc size>\n");
	__VERBOSE("         -k                     <kmer size>\n");
	__VERBOSE("         -o                     <output directory>\n");
	__VERBOSE("         --filter-threshold     <kmer count cut off>\n");
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
			if (V->bucks[i].idx == 0)
				fprintf(stderr, "Count????? %lu\n", V->bucks[i].cnt);
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

void dump_seq(uint64_t num, char *seq, int len)
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
	strcat(dump_path, "/graph.gfa");
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

void main_process(struct opt_count_t *opt)
{
	struct kmhash_t *V;
	__VERBOSE("Counting kmer...\n");
	V = count_kmer(opt);
	khash_t(kvert) *hvert;
	__VERBOSE("Filtering vertices...\n");
	hvert = filter_kmer(V, opt);
	kmhash_destroy(V);

	int16_t *edge_count;
	__VERBOSE("Counting edges...\n");
	edge_count = get_edges(opt, hvert);

	__VERBOSE("Dumping graph...\n");
	dump_graph(opt, hvert, edge_count);
}

void opt_process(int argc, char *argv[])
{
	struct opt_count_t *opt;
	opt = parse_count_option(argc - 1, argv + 1);
	char tmp_dir[1024];
	strcpy(tmp_dir, opt->out_dir); strcat(tmp_dir, "/count.log");
	init_log(tmp_dir);

	int cmd_len = 0, i;
	for (i = 0; i < argc; ++i)
		cmd_len += strlen(argv[i]) + 1;
	char *cmd = malloc(cmd_len);
	cmd_len = 0;
	for (i = 0; i < argc; ++i)
		cmd_len += sprintf(cmd + cmd_len, i + 1 == argc ? "%s" : "%s ", argv[i]);
	__VERBOSE_LOG("INFO", "command: \"%s\"\n", cmd);
	free(cmd);

	__VERBOSE_LOG("INFO", "kmer size: %d\n", opt->kmer_size);
	__VERBOSE_LOG("INFO", "pre-allocated hash table size: %d\n", opt->hash_size);
	__VERBOSE_LOG("INFO", "number of threads: %d\n", opt->n_threads);
	__VERBOSE_LOG("INFO", "cut off with kmer count less or equal: %d\n", opt->filter_thres);
	if (opt->n_files == 0) {
		__VERBOSE_LOG("INFO", "input: { stdin }\n");
	} else {
		if (opt->files_2 == NULL) {
			int len = 10, i;
			for (i = 0; i < opt->n_files; ++i)
				len += strlen(opt->files_1[i]) + 2;
			char *list_files = malloc(len);
			len = 0;
			len += sprintf(list_files, "{ ");
			for (i = 0; i < opt->n_files; ++i)
				len += sprintf(list_files + len,
						i + 1 == opt->n_files ? "%s" : "%s, ",
						opt->files_1[i]);
			sprintf(list_files + len, " }");
			__VERBOSE_LOG("INFO", "input: %s\n", list_files);
			free(list_files);
		} else {
			int len = 10, i;
			for (i = 0; i < opt->n_files; ++i)
				len += strlen(opt->files_1[i]) + strlen(opt->files_2[i]) + 6;
			char *list_files = malloc(len);
			len = 0;
			len += sprintf(list_files, "{ ");
			for (i = 0; i < opt->n_files; ++i)
				len += sprintf(list_files + len,
						i + 1 == opt->n_files ? "(%s, %s)" : "(%s, %s), ",
						opt->files_1[i], opt->files_2[i]);
			sprintf(list_files + len, " }");
			__VERBOSE_LOG("INFO", "input: %s\n", list_files);
			free(list_files);
		}
	}
	main_process(opt);
}

int main(int argc, char *argv[])
{
	if (argc < 4) {
		print_usage(argv[0]);
		return -1;
	}
	opt_process(argc, argv);
	return 0;
}

