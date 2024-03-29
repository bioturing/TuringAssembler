#include <fts.h>
#include <errno.h>
#include "helper.h"
#include "utils.h"
#include "log.h"
#include "unistd.h"
char __base[5] = {'A', 'C', 'G', 'T', 'N'};

char int_to_base(int x)
{
	return __base[x];
}

int base_to_int(char ch)
{
	return nt4_table[(int) ch];
}

char flip(char ch)
{
	return rev_nt4_char[nt4_table[(int) ch]];
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

void swap(void *a, void *b, int size)
{
	void *c = calloc(1, sizeof(size));
	memcpy(c, a, size);
	memcpy(a, b, size);
	memcpy(b, c, size);
	free(c);
}

int min(int a, int b)
{
	return a < b ? a : b;
}

int max(int a, int b)
{
	return a > b ? a : b;
}

// https://stackoverflow.com/questions/2256945/removing-a-non-empty-directory-programmatically-in-c-or-c
int recursive_delete_files(char *dir)
{
	int ret = 0;
	FTS *ftsp = NULL;
	FTSENT *curr;

	char *files[] = {dir, NULL};

	ftsp = fts_open(files, FTS_NOCHDIR | FTS_PHYSICAL | FTS_XDEV, NULL);
	if (!ftsp){
		log_debug("%s: fts open failed: %s", dir, strerror(errno));
		ret = -1;
		goto finish;
	}
	while ((curr = fts_read(ftsp))){
		switch(curr->fts_info){
			case FTS_NS:
			case FTS_DNR:
			case FTS_ERR:
				log_debug("%s fts_read error: %s",
					curr->fts_accpath, strerror(errno));
				break;
			case FTS_DC:
			case FTS_DOT:
			case FTS_NSOK:
				break;
			case FTS_D:
				break;
			case FTS_DP:
			case FTS_F:
			case FTS_SL:
			case FTS_SLNONE:
			case FTS_DEFAULT:
				if (remove(curr->fts_accpath) < 0){
					log_debug("%s: Failed to remove %s",
						curr->fts_path, strerror(errno));
					ret = -1;
				}
				break;
		}
	}
finish:
	if (ftsp)
		fts_close(ftsp);
	return ret;
}

int check_file_exist(char *path)
{
	return access(path, F_OK) != -1;
}

int delete_file(char *file_path)
{
	int status = remove(file_path);
	if (status == 0)
		log_debug("%s deleted succesfully", file_path);
	else
		log_debug("Cannot delete %s", file_path);
}
