/**
 * Copyright (C) Ghost Robotics - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Avik De <avik@ghostrobotics.io>
 */
#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#include "main.h"

#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define STDERR_FILENO 2

extern UART_HandleTypeDef huart1;

extern int _write(int file, char *data, int len)
{
	// The special files number 0, 1, 2 are for stdin, stdout, stderr (typically)
	switch (file)
	{
	case STDIN_FILENO:
	case STDOUT_FILENO:
	case STDERR_FILENO:
		return HAL_UART_Transmit_IT(&huart1, (uint8_t *)data, len);
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
