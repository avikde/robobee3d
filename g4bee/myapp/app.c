
#include "main.h"
#include "cmsis_os.h"
#include <stdio.h>
#include "bos1901.h"

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

	for (;;)
	{
		// first 4 bits are address
		bos.txBuf[0] = 0; // set SDO to have 0x5
		bos1901rw(&bos);

		osDelay(100);
	}
}
