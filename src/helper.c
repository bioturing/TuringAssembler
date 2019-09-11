#include "helper.h"
char __base[5] = {'A', 'C', 'G', 'T', 'N'};

char int_to_base(int x)
{
	return __base[x];
}

int base_to_int(char ch)
{
	for (int i = 0; i < 4; ++i)
		if (ch == __base[i])
			return i;
	return 0;
}

char flip(char ch)
{
	for (int i = 0; i < 4; ++i)
		if (__base[i] == ch)
			return __base[i ^ 3];
	return 'N';
}

void flip_reverse(char *seq)
{
	int len = strlen(seq);
	for (int i = 0; i < len / 2; ++i){
		char a = seq[i];
		char b = seq[len - i - 1];
		seq[i] = flip(b);
		seq[len - i - 1] = flip(a);
	}
	if (len % 2)
		seq[len / 2] = flip(seq[len / 2]);
}

void decode_seq(char **dest, uint32_t *source, int len)
{
	*dest = (char *) calloc(len + 1, sizeof(char));
	for (int i = 0; i < len; ++i)
		(*dest)[i] = int_to_base(__binseq_get(source, i));
}

void encode_seq(uint32_t **dest, char *source)
{
	int len = strlen(source);
	*dest = (uint32_t *) calloc((len + 15) / 16, sizeof(uint32_t));
	int tmp = 0;
	for (int i = 0; i < len; ++i){
		__binseq_set(*dest, tmp, base_to_int(source[i]));
		++tmp;
	}
}

float get_cov(struct asm_graph_t g, int edge)
{
	return (float) g.edges[edge].count / (g.edges[edge].seq_len - g.ksize);
}

void swap(void *a, void *b, int size)
{
	void *c = calloc(1, sizeof(size));
	memcpy(c, a, size);
	memcpy(a, b, size);
	memcpy(b, c, size);
	free(c);
}

void get_dump_N(char **N)
{
	*N = (char *) calloc(101, sizeof(char));
	for (int i = 0; i < 100; ++i)
		(*N)[i] = 'N';
}

int min(int a, int b)
{
	return a < b ? a : b;
}

int max(int a, int b)
{
	return a > b ? a : b;
}
