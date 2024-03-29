#ifndef __ASM_HELPER__
#define __ASM_HELPER__
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "assembly_graph.h"
char int_to_base(int x);
int base_to_int(char base);
int min(int a, int b);
int max(int a, int b);
void flip_reverse(char *seq);
void encode_seq(uint32_t **dest, char *source);
void decode_seq(char **dest, uint32_t *source, int len);
char flip(char base);
void swap(void *a, void *b, int size);
int recursive_delete_files(char *dir);
int check_file_exist(char *path);
int delete_file(char *file_path);
#endif
