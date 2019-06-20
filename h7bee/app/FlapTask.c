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
extern ADC_HandleTypeDef hadc1;

static void voltageControl(float vdes, float vact, TIM_HandleTypeDef *htim)
{
	const uint16_t arr = 10000; //depends on ARR setting

	// hi-side, and lo-side duty cycles. only one of them can be > 0 for each period
	static float dch, dcl;
	// TODO: on time (duty cycle) related to the magnitude of the difference?
	float mag = 0.02 * (vact - vdes);
	if (vact > vdes)
	{
		dcl = mag;
		dch = 0;
	}
	else
	{
		dch = -mag;
		dcl = 0;
	}
	__HAL_TIM_SetCompare(htim, TIM_CHANNEL_1, (uint16_t)(arr * dch));
	__HAL_TIM_SetCompare(htim, TIM_CHANNEL_2, (uint16_t)(arr * dcl));
}

// This is called from a timer update of the PWM generating timer (see *_it.c)
void flapUpdate(void const *argument)
{
	HAL_GPIO_TogglePin(LED1_GPIO_Port, LED1_Pin); // timing debugging

	float sfreq = 100;
	float t = 0.001 * millis() + 0.000001 * (micros() % 1000);
	float vdes = 0.5 * (1 + sinf(2 * PI * sfreq * t));
	// This is the "reference"
	// TODO: look at ADC for V1, V2
	volatile uint32_t a1 = ADC1->JDR1;
	volatile uint32_t a2 = ADC1->JDR2;
	volatile uint32_t a3 = ADC1->JDR3;
	volatile uint32_t a4 = ADC1->JDR4;
	HAL_ADCEx_InjectedStart(&hadc1);
	
	voltageControl(vdes, 0.5, &htim3);
	voltageControl(vdes, 0.5, &htim4);
}
