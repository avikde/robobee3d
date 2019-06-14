/**
 * @file FlapTask.c
 * @author Avik De (avikde@gmail.com)
 * @brief Task for flap to flap control
 * @version 0.1
 * @date 2019-06-14
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "stm32h7xx_hal.h"
#include "main.h"
#include <stdint.h>
#include "osal.h"

extern TIM_HandleTypeDef htim3, htim4;

static void setPWM(float h1, float l1, float h2, float l2)
{
	uint16_t arr = 1000; // out of 1000 (depending on ARR setting)
	__HAL_TIM_SetCompare(&htim3, TIM_CHANNEL_1, (uint16_t)(arr * h1));
	__HAL_TIM_SetCompare(&htim3, TIM_CHANNEL_2, (uint16_t)(arr * l1));
	__HAL_TIM_SetCompare(&htim4, TIM_CHANNEL_1, (uint16_t)(arr * h2));
	__HAL_TIM_SetCompare(&htim4, TIM_CHANNEL_2, (uint16_t)(arr * l2));
}

void flapUpdate(void const *argument)
{
	HAL_GPIO_TogglePin(LED1_GPIO_Port, LED1_Pin); // timing debugging

	float sfreq = 150;
	float t = 0.001 * millis() + 0.000001 * (micros() % 1000);
	float spwm = 0.5 * (1 + sinf(2 * PI * sfreq * t));
	setPWM(0, spwm, 0, spwm);
}
