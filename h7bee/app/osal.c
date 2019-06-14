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
	static int avail;
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
uint32_t micros()
{
	return (_millis * 1000) + TIM6->CNT;
}
