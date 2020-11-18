#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#include <string.h>
#include "main.h"

#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define STDERR_FILENO 2

extern UART_HandleTypeDef huart2;
// Possibly the user could be asked to just not call _write again before the previous transmit finished
// stdout has a buffer so as long as there are no simultaneous newlines this would be fine
#define WRITEBUF_SIZE 256
uint8_t _writeBuf[WRITEBUF_SIZE];

extern int _write(int file, char *data, int len)
{
	// The special files number 0, 1, 2 are for stdin, stdout, stderr (typically)
	switch (file)
	{
	case STDIN_FILENO:
	case STDOUT_FILENO:
	case STDERR_FILENO:
		// yield and wait if there is a TX ongoing. This is what Transmit_IT checks too.
		while (huart2.gState != HAL_UART_STATE_READY)
			HAL_Delay(1);

		// Buffer the data so that printf can use it again
		memcpy(_writeBuf, data, len);
		HAL_UART_Transmit_IT(&huart2, _writeBuf, len);
		return len;
	default:
		return EBADF;
	}
}

extern int _read(int file, char *ptr, int len)
{
	switch (file)
	{
	case STDIN_FILENO:
	case STDOUT_FILENO:
	case STDERR_FILENO:
		//TODO: 
		return 0;

	default:
		return EBADF;
	}
}

int _gettimeofday(struct timeval *tv, void *tzvp)
{
	uint32_t t = HAL_GetTick(); // in ms
	tv->tv_sec = t / 1000;
	tv->tv_usec = (t % 1000) * 1000;
	return 0;
}
