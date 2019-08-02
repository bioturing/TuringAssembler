#ifndef _UTILS_H_
#define _UTILS_H_

/* No verbose/log/error macros using here */

#if defined(_MSC_VER)
#pragma warning(disable:4996)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include "pthread_barrier.h"
// #include <stdint.h>
// #include <stdarg.h>
// #include <math.h>
// #include <errno.h>
// #include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
// #include <pthread.h>
// #if defined(_MSC_VER)
// #include <time.h>
// #include <windows.h>
// #include <getopt.h>
// #include <BaseTsd.h>
// typedef SSIZE_T ssize_t;
// #else
// #include <semaphore.h>
// #include <unistd.h>
// #include <sys/resource.h>
// #include <sys/time.h>
// #endif

// #include "attribute.h"
// #include "khash.h"


/* color terminal */

#define MAX_INT32		2147483647
#define MIN_INT32		-2147483648

#define MASK32			4294967295ULL

#define BUFSZ			4096

#define THREAD_STACK_SIZE	16777216

#define FORWARD			0
#define REVERSE			1
#define LEFT			0
#define RIGHT			1

//#define WRITE_TEXT
//#define RUN_AFTER_ANALYSIS

/*
 * Built in macros
 */

#define __abs(x) 		((x) < 0 ? -(x) : (x))

#define __min(a, b) 		((a) < (b) ? (a) : (b))

#define __max(a, b) 		((a) > (b) ? (a) : (b))

#define __min3(a, b, c)		__min(__min((a), (b)), (c))

#define __max3(a, b, c)		__max(__max((a), (b)), (c))

#define __round_up_32(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16, ++(x))

#define __round_up_64(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16,	       \
				 (x) |= (x) >> 32, ++(x))

#define __is_sep(c)		((c) == ' ' || (c) == '\t')

#define __rotl64(x, r) (((x) << (r)) | ((x) >> (64 - (r))))

#define normalize_mapq(s)	do {					       \
	if ((s) < 0) (s) = 0;						       \
	if ((s) > 60) (s) = 60;						       \
} while (0)

/*
 * Built-in macros function
 */

#define __ALLOC(ptr, sz)	(ptr) = xmalloc(sizeof(*(ptr)) * (sz))

#define __REALLOC(ptr, sz)	(ptr) = xrealloc((ptr), sizeof(*(ptr)) * (sz))

/* push back val to ptr, ptr has sz element, realloc + 1 */
#define __PUSH_BACK(ptr, sz, val) do {					       \
	assert((sz) >= 0);						       \
	__REALLOC((ptr), (sz) + 1);					       \
	(ptr)[(sz)++] = (val);						       \
} while(0)

#define __FREE_AND_NULL(ptr) do {					       \
	free(p);							       \
	(p) = NULL;							       \
} while (0)

#define __SWAP(x, y) do {						       \
	assert(sizeof(x) == sizeof(y));					       \
	int8_t temp[sizeof(x)];						       \
	memcpy(temp, &(y), sizeof(x));					       \
	memcpy(&(y), &(x), sizeof(x));					       \
	memcpy(&(x), temp, sizeof(x));					       \
} while (0)

/*
 * Built-in function
 */

/* get time */
double realtime();

/* reverse compelemnt */
char *get_rev_complement(const char *seq, int len);

/* reverse string */
char *get_rev(const char *seq, int len);

/* return new char* concate s1 and s2 */
char *str_concate(const char *s1, const char *s2);

/* convert from [ACGTN]+ seq to number */
int64_t seq2num(const char *seq, int len);

/* convert from number to [ACGTN]+ to number */
char *num2seq(int64_t num, int len);
/*
 * Global variable
 */

extern int8_t nt4_table[256];
extern char *nt4_char, *rev_nt4_char;

static inline uint32_t unpack_int32(uint8_t *buf)
{
	uint32_t ret = 0;
	ret =	(uint32_t)buf[0]		|
		((uint32_t)buf[1] << 8)	|
		((uint32_t)buf[2] << 16)	|
		((uint32_t)buf[3] << 24);
	return ret;
}

static inline uint64_t unpack_int64(uint8_t *buf)
{
	uint64_t ret = 0;
	ret =	(uint64_t)buf[0]		|
		((uint64_t)buf[1] << 8)	|
		((uint64_t)buf[2] << 16)	|
		((uint64_t)buf[3] << 24)	|

		((uint64_t)buf[4] << 32)	|
		((uint64_t)buf[5] << 40)	|
		((uint64_t)buf[6] << 48)	|
		((uint64_t)buf[7] << 56);
	return ret;
}

static inline void pack_int32(uint8_t *buffer, uint32_t value)
{
	buffer[0] = value;
	buffer[1] = value >> 8;
	buffer[2] = value >> 16;
	buffer[3] = value >> 24;
}

static inline void pack_int64(uint8_t *buffer, uint64_t value)
{
	buffer[0] = value;
	buffer[1] = value >> 8;
	buffer[2] = value >> 16;
	buffer[3] = value >> 24;

	buffer[4] = value >> 32;
	buffer[5] = value >> 40;
	buffer[6] = value >> 48;
	buffer[7] = value >> 56;
}

#endif /* _UTILS_H_ */
