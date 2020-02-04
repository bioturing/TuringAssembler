//
// Created by che on 04/02/2020.
//


#include "log.h"

void log_info_wrap(const char *fmt)
{
	log_info(fmt);
}

void log_warn_wrap(const char *fmt)
{
	log_warn(fmt);
}

void log_debug_wrap(const char *fmt)
{
	log_debug(fmt);
}
