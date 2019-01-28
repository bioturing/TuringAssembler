#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include "time_utils.h"

static double time_start, time_current;

static double get_current_time()
{
	struct timeval tp;
#if defined(_MSC_VER)
	_gettimeofday(&tp,  NULL);
#else
	gettimeofday(&tp,  NULL);
#endif
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

void init_clock()
{
	time_start = time_current = get_current_time();
}

void set_time_now()
{
	time_current = get_current_time();
}

double sec_from_prev_time()
{
	return get_current_time() - time_current;
}

double sec_from_initial_time()
{
	return get_current_time() - time_start;
}
