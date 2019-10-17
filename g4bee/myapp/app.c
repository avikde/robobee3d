
#include "main.h"
#include "cmsis_os.h"

void startBlinkTask(void *argument)
{
	for (;;)
	{
		HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
		volatile uint32_t millis = HAL_GetTick();
		// printf("hello %d\n", millis);
		osDelay(100);
	}
}
