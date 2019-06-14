/**
 * @file osal.c
 * @author Avik De (avikde@gmail.com)
 * @brief OS abstraction stuff
 * @version 0.1
 * @date 2019-06-14
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#include <stm32h7xx_hal.h>
#include "cmsis_os.h"
#include <string.h>

// stdio --------------------------------------------

#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define STDERR_FILENO 2

extern UART_HandleTypeDef huart1;
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
		while (huart1.gState != HAL_UART_STATE_READY)
			osDelay(1);

		// Buffer the data so that printf can use it again
		memcpy(_writeBuf, data, len);
		HAL_UART_Transmit_IT(&huart1, _writeBuf, len);
		return len;
	default:
		return EBADF;
	}
}

extern int _read(int file, char *ptr, int len)
{
	// static int avail;
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

// System clock ----------------------------------------------------

extern volatile uint32_t _millis;
uint32_t millis()
{
	return _millis;
}

uint32_t micros()
{
	return (_millis * 1000) + TIM6->CNT;
}

static void delayMicroseconds(uint32_t n)
{
	// fudge for function call overhead
	if (n > 0)
	{
		n--;
		uint32_t start = micros();
		while (micros() - start < n)
			;
	}
}

int osal_usleep(uint32_t usec)
{
	if (usec < 1000)
		delayMicroseconds(usec);
	else
	{
		delayMicroseconds(usec % (uint32_t)1000);
		HAL_Delay(usec / 1000);
	}

	return 1;
}

// Other tasks -------------------------------------------------------

extern float extest, eytest, eztest;

void startPrintfTask(void const *argument)
{
	for (;;)
	{
		printf("hi %d %d\t", millis(), micros());
		printf("%.2f\t%.2f\t%.2f\t", extest, eytest, eztest);
		printf("\n");
		osDelay(20);
	}
}
