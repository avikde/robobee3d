
#include "stm32h7xx_hal.h"
#include "main.h"

extern TIM_HandleTypeDef htim3, htim4;

void flapUpdate(void const *argument)
{
	HAL_GPIO_TogglePin(LED1_GPIO_Port, LED1_Pin);
	// out of 1000 (depending on ARR setting)
	__HAL_TIM_SetCompare(&htim3, TIM_CHANNEL_1, 100);
	__HAL_TIM_SetCompare(&htim3, TIM_CHANNEL_2, 200);
	__HAL_TIM_SetCompare(&htim4, TIM_CHANNEL_1, 300);
	__HAL_TIM_SetCompare(&htim4, TIM_CHANNEL_2, 400);
}
