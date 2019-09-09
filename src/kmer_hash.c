#include "kmer_hash.h"
#include "verbose.h"

khash_t(kmer_int) *get_kmer_hash(char *r1_path, char *r2_path, int ksize)
{
	struct read_list_t rlist1, rlist2;
	load_reads_from_fastq(r1_path, &rlist1);
	load_reads_from_fastq(r2_path, &rlist2);

	khash_t(kmer_int) *res = kh_init(kmer_int);

	for (int i = 0; i < rlist1.n; ++i){
		count_hash_from_read(rlist1.reads[i].seq, ksize, res);
		char *rv = (char *) calloc(rlist1.reads[i].len + 1, sizeof(char));
		strcpy(rv, rlist1.reads[i].seq);
		flip_reverse(rv);
		count_hash_from_read(rv, ksize, res);
		free(rv);
	}
	for (int i = 0; i < rlist2.n; ++i){
		count_hash_from_read(rlist2.reads[i].seq, ksize, res);
		char *rv = (char *) calloc(rlist2.reads[i].len + 1, sizeof(char));
		strcpy(rv, rlist2.reads[i].seq);
		flip_reverse(rv);
		count_hash_from_read(rv, ksize, res);
		free(rv);
	}
	return res;
}

void count_hash_from_read(char *seq, int ksize, khash_t(kmer_int) *h)
{
	uint64_t power = 1;
	for (int i = 0; i < ksize - 1; ++i)
		power *= 4;

	uint64_t hash = 0;
	for (int i = 0; i < ksize - 1; ++i)
		hash = hash * 4 + base_to_int(seq[i]);
	int len = strlen(seq);
	for (int i = 0; i + ksize <= len; ++i){
		if (i > 0)
			hash -= power * base_to_int(seq[i - 1]);
		hash = hash * 4 + base_to_int(seq[i + ksize - 1]);

		khiter_t it = kh_get(kmer_int, h, hash);
		if (it == kh_end(h)){
			int ret;
			it = kh_put(kmer_int, h, hash, &ret);
			kh_val(h, it) = 0;
		}
		++kh_val(h, it);
	}
}

int count_kmer_on_seq(khash_t(kmer_int) *h, char *seq, int ksize)
{
	khash_t(kmer_int) *h1 = kh_init(kmer_int);
	count_hash_from_read(seq, ksize, h1);
	int res = 0;
	for (khiter_t it1 = kh_begin(h1); it1 != kh_end(h1); ++it1){
		if (!kh_exist(h1, it1))
			continue;
		uint64_t key = kh_key(h1, it1);
		int val = kh_val(h1, it1);
		khiter_t it = kh_get(kmer_int, h, key);
		if (it != kh_end(h))
			//res += min(val, kh_val(h, it));
			res += kh_val(h, it);
	}
	kh_destroy(kmer_int, h1);
	return res;
}

void print_kmer_count_on_seq(khash_t(kmer_int) *h, char *seq, int ksize)
{
	uint64_t power = 1;
	for (int i = 0; i < ksize - 1; ++i)
		power *= 4;

	uint64_t hash = 0;
	for (int i = 0; i < ksize - 1; ++i)
		hash = hash * 4 + base_to_int(seq[i]);
	int len = strlen(seq);
	for (int i = 0; i + ksize <= len; ++i){
		if (i > 0)
			hash -= power * base_to_int(seq[i - 1]);
		hash = hash * 4 + base_to_int(seq[i + ksize - 1]);
		khiter_t it = kh_get(kmer_int, h, hash);
		if (it != kh_end(h))
			printf("%d ", kh_val(h, it));
		else
			printf("0 ");
	}
	printf("\n");
}

int kmer_check(char *first, char *second, int overlap_ksize, int check_ksize,
		khash_t(kmer_int) *kmer_count)
{
	int len1 = strlen(first);
	int len2 = strlen(second);
	if (len1 + len2 - overlap_ksize < check_ksize)
		return 0;
	int m = min(check_ksize, len1);
	char *join = (char *) calloc(m + len2 - overlap_ksize + 1,
			sizeof(char));
	int pos = len1 - m;
	strcpy(join, first + pos);
	strcat(join, second + overlap_ksize);
	khash_t(kmer_int) *h = kh_init(kmer_int);
	count_hash_from_read(join, check_ksize, h);
	int res = 1;
	for (khiter_t it = kh_begin(h); it != kh_end(h); ++it){
		if (!kh_exist(h, it))
			continue;
		uint64_t key = kh_key(h, it);
		khiter_t it2 = kh_get(kmer_int, kmer_count, key);
		if (it2 == kh_end(kmer_count)){
			res = 0;
			break;
		}
	}
	kh_destroy(kmer_int, h);
	free(join);
	return res;
}

int count_zero_kmer_map(char *first, char *second, int overlap_ksize,
		int check_ksize, khash_t(kmer_int) *kmer_count)
{
	int len1 = strlen(first);
	int len2 = strlen(second);
	if (len1 + len2 - overlap_ksize < check_ksize)
		return 0;
	int m = min(check_ksize, len1);
	char *join = (char *) calloc(m + len2 - overlap_ksize + 1,
			sizeof(char));
	int pos = len1 - m;
	strcpy(join, first + pos);
	strncpy(join + m, second + overlap_ksize,
			max(0, min(len2 - overlap_ksize, m)));
	khash_t(kmer_int) *h = kh_init(kmer_int);
	count_hash_from_read(join, check_ksize, h);
	int count = 0;
	for (khiter_t it = kh_begin(h); it != kh_end(h); ++it){
		if (!kh_exist(h, it))
			continue;
		uint64_t key = kh_key(h, it);
		khiter_t it2 = kh_get(kmer_int, kmer_count, key);
		if (it2 == kh_end(kmer_count))
			++count;
	}
	/*int total_kmer = strlen(join) - check_ksize + 1;
	float res = (float) count / total_kmer / total_kmer;*/
	kh_destroy(kmer_int, h);
	free(join);
	return count;
}

int count_max_consecutive_zero_kmer(char *first, char *second, int overlap_ksize,
		int check_ksize, khash_t(kmer_int) *kmer_count)
{
	int len1 = strlen(first);
	int len2 = strlen(second);
	if (len1 + len2 - overlap_ksize < check_ksize)
		return 0;
	int m = min(check_ksize, len1);
	char *join = (char *) calloc(m + len2 - overlap_ksize + 1,
			sizeof(char));
	int pos = len1 - m;
	strcpy(join, first + pos);
	strncpy(join + m, second + overlap_ksize,
			min(len2, check_ksize) - overlap_ksize);
	khash_t(kmer_int) *h = kh_init(kmer_int);
	count_hash_from_read(join, check_ksize, h);
	int max_consecutive = 0;
	int count = 0;
	for (khiter_t it = kh_begin(h); it != kh_end(h); ++it){
		if (!kh_exist(h, it))
			continue;
		uint64_t key = kh_key(h, it);
		khiter_t it2 = kh_get(kmer_int, kmer_count, key);
		if (it2 == kh_end(kmer_count))
			++count;
		else
			count = 0;
		max_consecutive = max(max_consecutive, count);
	}
	/*int total_kmer = strlen(join) - check_ksize + 1;
	float res = (float) count / total_kmer / total_kmer;*/
	kh_destroy(kmer_int, h);
	free(join);
	return max_consecutive;
}

