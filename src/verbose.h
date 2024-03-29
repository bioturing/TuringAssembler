#ifndef _VERBOSE_H_
#define _VERBOSE_H_

#include <stdio.h>
#include <unistd.h>

#include "log.h"

#if defined(_MSC_VER)
#define __VERBOSE_INFO(tag, fmt, ...) do {				       \
	fprintf(stderr, "[" tag "] " fmt, __VA_ARGS__);			       \
	fflush(stderr);							       \
} while (0) /* VERBOSE_INFO */

#define __VERBOSE_LOG(tag, fmt, ...) do {				       \
	fprintf(stderr, "[" tag "] " fmt, __VA_ARGS__);			       \
	log_write("[" tag "] " fmt, __VA_ARGS__);				       \
} while (0) /* VERBOSE_AND_LOG */

#define __VERBOSE(fmt, ...) do {					       \
	fprintf(stderr, fmt, __VA_ARGS__);				       \
	fflush(stderr);							       \
} while (0) /* VERBOSE */

#if defined(NDEBUG)
#define __DEBUG(fmt, ...) 0
#else
#define __DEBUG(fmt, ...) do {						       \
	fprintf(stderr, "[DEBUG] " fmt, __VA_ARGS__);			       \
	fflush(stderr);							       \
} while (0) /* __DEBUG */
#endif /* NDEBUG */

#else /* __MSC_VER */

#define __VERBOSE_INFO(tag, fmt, args...) do {				       \
	fprintf(stderr, "[" tag "] " fmt, ##args);				       \
	fflush(stderr);							       \
} while (0) /* VERBOSE_INFO */

#define __VERBOSE_LOG(tag, fmt, args...)

#define __LOG(tag, fmt, args...) do {				       \
	log_write("[" tag "] " fmt, ##args);					       \
} while (0) /* LOG */

#define __VERBOSE(fmt, args...) do {					       \
	fprintf(stderr, fmt, ##args);					       \
	fflush(stderr);							       \
} while (0) /* VERBOSE */

#if defined(NDEBUG)
#define __DEBUG(fmt, ...) 0
#else
#define __DEBUG(fmt, args...)	fprintf(stderr, "[DEBUG] " fmt, ##args)
#endif /* NDEBUG */

#endif /* __MSC_VER */

void init_log(const char *path);

void log_write(const char *fmt, ...);

void close_log();

#endif /* _VERBOSE_H_ */
