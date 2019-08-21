#include "map_contig.h"

void init_map_contig(struct map_contig_t *mct, struct asm_edge_t global_edge,
		struct asm_graph_t local_graph)
{
	mct->global_edge = global_edge;
	mct->local_graph = local_graph;
	mct->kmers = NULL;
	mct->pos = 0;
	mct->n_candidates = local_graph.n_e;
	mct->best_match = -1;
	init_local_kmers(mct);
}

void init_local_kmers(struct map_contig_t *mct)
{
	mct->kmers = (khash_t(int32_int) **) calloc(mct->n_candidates,
			sizeof(khash_t(int32_int) *));
	for (int i = 0; i < mct->local_graph.n_e; ++i){
		struct asm_edge_t cur_edge = mct->local_graph.edges[i];
		get_all_edge_kmers(cur_edge, mct->kmers + i);
	}
}

void get_all_seq_kmers(char *seq, int len, khash_t(int32_int) **kmers)
{
	uint32_t power = 1;
	for (int i = 0; i < KSIZE - 1; ++i)
		power *= HASH_BASE;
	*kmers = kh_init(int32_int);
	khint32_t hash = 0;
	for (int i = 0; i < KSIZE - 1; ++i){
		char base = seq[i];
		hash = hash * 4 + base_to_int(base);
	}
	for (int i = 0; i + KSIZE <= len; ++i){
		char cur_base = seq[i + KSIZE - 1];
		char pre_base = i > 0 ? seq[i - 1] : 0;
		hash = (hash - base_to_int(pre_base) * power) * HASH_BASE +
			base_to_int(cur_base);
		add_kmer(hash, *kmers);
	}
}

void get_all_edge_kmers(struct asm_edge_t edge, khash_t(int32_int) **kmers)
{
	char *seq;
	decode_seq(&seq, edge.seq, edge.seq_len);
	get_all_seq_kmers(seq, edge.seq_len, kmers);
	free(seq);
}

int find_match_from_pos(struct map_contig_t *mct)
{
	khash_t(int32_int) *kmers;
	int len = min(WINDOW_SIZE, mct->global_edge.seq_len - mct->pos);
	if (len == 0)
		return -1;
	char *seq;
	decode_seq(&seq, mct->global_edge.seq, mct->global_edge.seq_len);
	get_all_seq_kmers(seq + mct->pos, len, &kmers);

	int res = -1;
	int max_point = 0;
	int match_id = -1;
	int *points = (int *) calloc(mct->n_candidates, sizeof(int));
	for (int i = 0; i < mct->n_candidates; ++i){
		if (mct->local_graph.edges[i].seq_len < WINDOW_SIZE)
			continue;
		points[i] = count_match_kmer(kmers, mct->kmers[i]);
		max_point = max(max_point, points[i]);

		if (check_good_match(points[i], POINT_MEDIUM_TRHESH)){
			match_id = i;
			res = i;
			break;
		}
	}

	free(seq);
	free(points);
	kh_destroy(int32_int, kmers);
	return res;
}

int count_match_kmer(khash_t(int32_int) *first, khash_t(int32_int) *second)
{
	int res = 0;
	for (khiter_t it = kh_begin(first); it != kh_end(first); ++it){
		if (!kh_exist(first, it))
			continue;
		khiter_t it2 = kh_get(int32_int, second, kh_key(first, it));
		if (it2 != kh_end(second))
			res += min(kh_val(first, it), kh_val(second, it2));
	}
	return res;
}

int find_match(struct map_contig_t *mct)
{
	while (mct->pos < (int) mct->global_edge.seq_len){
		int res = find_match_from_pos(mct);
		if (res != -1){
			mct->best_match = res;
			return res;
		}
		advance_pos(mct);
	}
	return -1;
}

void map_contig_destroy(struct map_contig_t *mct)
{
	for (int i = 0; i < mct->n_candidates; ++i)
		kh_destroy(int32_int, mct->kmers[i]);
}

void advance_pos(struct map_contig_t *mct)
{
	int len = min(WINDOW_SIZE, mct->global_edge.seq_len - mct->pos);
	mct->pos += len;
}

void get_global_match_pos(struct map_contig_t *mct, struct subseq_pos_t *pos)
{
	pos->start = -1;
	pos->end = -1;
	if (mct->pos == (int) mct->global_edge.seq_len)
		return;
	pos->start = mct->pos;
	int cur_best_match = mct->best_match;
	while (check_stop(mct) == 0){
		int tmp = find_match_from_pos(mct);
		if (tmp != cur_best_match)
			break;
		int next_len = min(WINDOW_SIZE, mct->global_edge.seq_len
				- mct->pos);
		if (next_len == WINDOW_SIZE)
			pos->end = mct->pos;
		advance_pos(mct);
	}
	mct->best_match = cur_best_match;

	int len = mct->local_graph.edges[cur_best_match].seq_len;
	while (pos->end - pos->start + WINDOW_SIZE > len)
		pos->end -= WINDOW_SIZE;
	pos->end = max(pos->end, pos->start);
}

void get_match_pos(struct map_contig_t *mct, struct subseq_pos_t *global,
		struct subseq_pos_t *local)
{
	get_global_match_pos(mct, global);
	get_local_match_pos(mct, global, local);
}

void get_local_match_pos(struct map_contig_t *mct, struct subseq_pos_t *global,
		struct subseq_pos_t *local)
{
	char *global_seq;
	decode_seq(&global_seq, mct->global_edge.seq, mct->global_edge.seq_len);

	khash_t(int32_int) *global_start, *global_end;
	get_all_seq_kmers(global_seq + global->start,
			get_next_len_global(mct, global->start), &global_start);
	get_all_seq_kmers(global_seq + global->end,
			get_next_len_global(mct, global->end), &global_end);
	free(global_seq);

	char *local_seq;
	struct asm_edge_t best_match = mct->local_graph.edges[mct->best_match];
	if (best_match.seq_len < WINDOW_SIZE)
		__ERROR("Local contig is too short!\n");
	decode_seq(&local_seq, best_match.seq, best_match.seq_len);

	khash_t(int32_int) *local_kmers;
	get_all_seq_kmers(local_seq, WINDOW_SIZE, &local_kmers);

	int *start_point = (int *) calloc(best_match.seq_len, sizeof(int));
	int *end_point = (int *) calloc(best_match.seq_len, sizeof(int));
	int cur_start_point = count_match_kmer(global_start, local_kmers);
	int cur_end_point = count_match_kmer(global_end, local_kmers);
	
	for (khiter_t it = kh_begin(local_kmers); it != kh_end(local_kmers);
			++it){
		if (!kh_exist(local_kmers, it))
			continue;
		khint32_t key = kh_key(local_kmers, it);
		int val = kh_val(local_kmers, it);
		khiter_t it2;
		it2 = kh_get(int32_int, global_start, key);
		if (it2 != kh_end(global_start)){
			int min_val = min(kh_val(global_start, it2), val);
			kh_val(global_start, it2) -= min_val;
		}

		it2 = kh_get(int32_int, global_end, key);
		if (it2 != kh_end(global_end)){
			int min_val = min(kh_val(global_end, it2), val);
			kh_val(global_end, it2) -= min_val;
		}
	}
	int found = 0;
	for (int i = 0; i + WINDOW_SIZE <= (int) best_match.seq_len; ++i){
		start_point[i] = cur_start_point;
		end_point[i] = cur_end_point;
		if (i > 0){
			khint32_t pre_hash = get_one_seq_kmer_hash(local_seq
					+ i - 1);
			khiter_t it;
 			it = kh_get(int32_int, global_start, pre_hash);
			if (it != kh_end(global_start)
				&& kh_val(global_start, it) >= 0)
				--cur_start_point;

			it = kh_get(int32_int, global_end, pre_hash);
			if (it != kh_end(global_end)
				&& kh_val(global_end, it) >= 0)
				--cur_end_point;

			add_kmer(pre_hash, global_start);
			add_kmer(pre_hash, global_end);
		}
		if (i + WINDOW_SIZE < (int) best_match.seq_len){
			khint32_t next_hash = get_one_seq_kmer_hash(local_seq
					+ i + WINDOW_SIZE - KSIZE + 1);
			khiter_t it;
			it = kh_get(int32_int, global_start, next_hash);
			if (it != kh_end(global_start)
				&& kh_val(global_start, it) >= 1)
				++cur_start_point;

			it = kh_get(int32_int, global_end, next_hash);
			if (it != kh_end(global_end)
				&& kh_val(global_end, it) == 1)
				++cur_end_point;
			remove_kmer(next_hash, global_start);
			remove_kmer(next_hash, global_end);
		}
	}

	float max_ratio = 0;
	local->start = 0;
	local->end = 0;
	for (int i = 0; i < (int) best_match.seq_len; ++i){
		if (start_point[i] > start_point[local->start])
			local->start = i;
		if (end_point[i] > end_point[local->end])
			local->end = i;
	}
	free(local_seq);
	kh_destroy(int32_int, global_start);
	kh_destroy(int32_int, global_end);
	kh_destroy(int32_int, local_kmers);
	free(start_point);
	free(end_point);
}

khint32_t get_one_seq_kmer_hash(char *seq)
{
	khint32_t res = 0;
	for (int i = 0; i < KSIZE; ++i)
		res = res * HASH_BASE + base_to_int(seq[i]);
	return res;
}

int check_good_match(int point, float thresh)
{
	return point >= thresh * (WINDOW_SIZE - KSIZE + 1);
}


int get_next_len_global(struct map_contig_t *mct, int pos)
{
	return min(WINDOW_SIZE, mct->global_edge.seq_len - pos);
}

int get_next_len_local(struct map_contig_t *mct, int pos)
{
	return min(WINDOW_SIZE, mct->local_graph.edges[mct->best_match].seq_len
			- pos);
}

void add_kmer(khint32_t hash, khash_t(int32_int) *kmers)
{
	khiter_t it = kh_get(int32_int, kmers, hash);
	if (it == kh_end(kmers)){
		int ret;
		it = kh_put(int32_int, kmers, hash, &ret);
		kh_val(kmers, it) = 0;
	}
	++kh_val(kmers, it);
}

void remove_kmer(khint32_t hash, khash_t(int32_int) *kmers)
{
	khiter_t it = kh_get(int32_int, kmers, hash);
	if (it == kh_end(kmers)){
		int ret;
		it = kh_put(int32_int, kmers, hash, &ret);
		kh_val(kmers, it) = 0;
	}
	--kh_val(kmers, it);
}

int check_stop(struct map_contig_t *mct)
{
	int len = min(WINDOW_SIZE, mct->global_edge.seq_len - mct->pos);
	if (len < WINDOW_SIZE)
		return 1;
	return 0;
}
