
#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#include <string.h>
#include "main.h"

int _gettimeofday(struct timeval *tv, void *tzvp)
{
	uint32_t t = HAL_GetTick(); // in ms
	tv->tv_sec = t / 1000;
	tv->tv_usec = (t % 1000) * 1000;
	return 0;
}
