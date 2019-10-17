
#include "main.h"
#include "cmsis_os.h"
#include <stdio.h>

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
	SPI_HandleTypeDef *spi = &hspi1;

	volatile uint8_t rxBuf[2] = {0, 0};
	uint8_t txBuf[2] = {0, 0};

	HAL_Delay(5000);
	HAL_GPIO_WritePin(SS1_GPIO_Port, SS1_Pin, GPIO_PIN_RESET);
	// first 4 bits are address
	txBuf[0] = (0x5) << 4 | (0x4 >> 1); // set SDO to have 0x4
	HAL_SPI_TransmitReceive(spi, txBuf, rxBuf, 2, 100);
	HAL_GPIO_WritePin(SS1_GPIO_Port, SS1_Pin, GPIO_PIN_SET);

	for (;;)
	{
		HAL_GPIO_WritePin(SS1_GPIO_Port, SS1_Pin, GPIO_PIN_RESET);
		// first 4 bits are address
		txBuf[0] = 0; // set SDO to have 0x5
		HAL_SPI_TransmitReceive(spi, txBuf, rxBuf, 2, 100);
		HAL_GPIO_WritePin(SS1_GPIO_Port, SS1_Pin, GPIO_PIN_SET);

		osDelay(100);
	}
}
