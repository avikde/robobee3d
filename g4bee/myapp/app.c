
#include "main.h"
#include "cmsis_os.h"
#include <stdio.h>
#include "bos1901.h"
#include <stdbool.h>
#include <math.h>

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

float rampWave(float phase)
{
	return 30 * sinf(2 * M_PI * phase);
	// return phase > 0.5 ? 0.0 : 20 * sinf(M_PI * phase);
}

void startBosTask(void *argument)
{
	BOS1901 bos;
	bos1901Init(&bos, &hspi1, SS1_GPIO_Port, SS1_Pin);
	bos1901Config(&bos, 0x0C, 0, 0); // read status

	for (;;)
	{
		bos1901AddWave(&bos, rampWave);
		// Play
		volatile uint16_t status = bos1901Config(&bos, 0x0C, 1, 0x8);
		// osDelay(10);

		// check if FIFO empty again?
		// volatile uint16_t dat1a = bos1901rw(&bos, 0, 0);
		// volatile bool full = dat1a & (0b1 << 14);
		osDelay(100);
	}
}
