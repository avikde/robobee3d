
#include "main.h"
#include "cmsis_os.h"
#include <stdio.h>

extern UART_HandleTypeDef huart2;

void startBlinkTask(void *argument)
{
	for (;;)
	{
		HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
		volatile uint32_t millis = HAL_GetTick();
		// printf("hello %d\n", millis);
		const char *_writeBuf = "hello";
		HAL_UART_Transmit_IT(&huart2, _writeBuf, 5);
		osDelay(100);
	}
}
