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
#include "math.h"
#include <stdbool.h>

extern TIM_HandleTypeDef htim3, htim4;
extern ADC_HandleTypeDef hadc1;
// FIXME: debugging
volatile float vmax[2] = {NAN, NAN}, vmin[2] = {NAN, NAN};
bool adcValidDataYet = false;

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

static void analogGetValues(float *vact, float *iact)
{
	static uint32_t vsens1, vsensl1, vsens2, vsensl2;
	// from hw design
	vsens1 = ADC1->JDR3;
	vsensl1 = ADC1->JDR4;
	vsens2 = ADC1->JDR1;
	vsensl2 = ADC1->JDR2;
	HAL_ADCEx_InjectedStart(&hadc1); // restart

	// if all 0, probably data not ready yet (invalid)
	if (!adcValidDataYet && vsens1 == 0 && vsensl1 == 0 && vsens2 == 0 && vsensl2 == 0)
		return;
	adcValidDataYet = true; // will short circuit the if above after this
	
	// output voltage TODO: scale
	vact[0] = vsens1;
	vact[1] = vsens2;
	// TODO: currents
}

// This is called from a timer update of the PWM generating timer (see *_it.c)
void flapUpdate(void const *argument)
{
	HAL_GPIO_TogglePin(LED1_GPIO_Port, LED1_Pin); // timing debugging

	float sfreq = 100;
	float t = 0.001 * millis() + 0.000001 * (micros() % 1000);
	float vdes = 0.5 * (1 + sinf(2 * PI * sfreq * t)); // This is the "reference"

	float vact[2], iact[2];
	analogGetValues(vact, iact);

	// store the min and max
	vmax[0] = fmaxf(vmax[0], vact[0]);
	vmax[1] = fmaxf(vmax[1], vact[1]);
	vmin[0] = fminf(vmin[0], vact[0]);
	vmin[1] = fminf(vmin[1], vact[1]);
	// TODO: look at ADC for V1, V2

	voltageControl(vdes, 0.5, &htim3);
	voltageControl(vdes, 0.5, &htim4);
}
