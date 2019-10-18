
#include "main.h"
#include "cmsis_os.h"
#include <stdio.h>
#include "bos1901.h"
#include <stdbool.h>

extern UART_HandleTypeDef huart2;
extern SPI_HandleTypeDef hspi1;

void startBlinkTask(void *argument)
{
	UART_HandleTypeDef *uart = &huart2;
	for (;;)
	{
		HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
		volatile uint32_t millis = HAL_GetTick();
		// printf("hello %d\n", millis);
		const char *_writeBuf = "hello";
		HAL_UART_Transmit_IT(uart, _writeBuf, 5);
		osDelay(100);
	}
}

void startBosTask(void *argument)
{
	BOS1901 bos;
	bos1901Init(&bos, &hspi1, SS1_GPIO_Port, SS1_Pin);
	bos1901Config(&bos, 0x1a, 0, 0);

	for (;;)
	{
		volatile uint16_t dat1a = bos1901rw(&bos, 0, 0);
		volatile bool full = dat1a & (0b1 << 14);
		osDelay(100);
	}
}
