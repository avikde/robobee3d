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

#define constrain(amt, low, high) ((amt) < (low) ? (low) : ((amt) > (high) ? (high) : (amt)))

static void voltageControl(float vdes, float vact, TIM_HandleTypeDef *htim)
{
	const uint16_t arr = 10000; //depends on ARR setting
	// hi-side, and lo-side duty cycles. only one of them can be > 0 for each period
	static float dch, dcl;
	// TODO: on time (duty cycle) related to the magnitude of the difference?
	float mag = constrain(0.005 * (vact - vdes), -0.05, 0.05);
	if (vact > vdes)
	{
		dcl = mag;
		dch = 0;
	}
	else
	{
		dch = -0.5 * mag; // based on looking at sine output
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

	// scale from https://docs.google.com/spreadsheets/d/1NQQbD_Zaig3STnTJG6g7wlEqw0pASV7lOaqTsuFVBwo/edit?usp=sharing
	const float VADC_C0 = -4.67, VADC_C1 = 0.0191;
	const float VSMOOTH = 0.3;
	// output voltage
	vact[0] = VSMOOTH * vact[0] + (1 - VSMOOTH) * (VADC_C0 + VADC_C1 * vsens1);
	vact[1] = VSMOOTH * vact[1] + (1 - VSMOOTH) * (VADC_C0 + VADC_C1 * vsens2);
	// TODO: currents
}

// This is called from a timer update of the PWM generating timer (see *_it.c) at 10KHz
void flapUpdate(void const *argument)
{
	const float UPDATE_DT = 0.0001; // 10KHz timebase update rate
	static float phase = 0, vdes;
	HAL_GPIO_TogglePin(LED1_GPIO_Port, LED1_Pin); // timing debugging

	// PARAMETERS
	float sfreq = 100; // wave freq
	float Vpp = 20;
	// --

	phase += sfreq * UPDATE_DT;
	// // keep between 0 and 1
	// while (phase >= 1)
	// 	phase -= 1;
	float p1 = fmodf(phase, 1.0);
	
	// Sinusoid
	vdes = 0.5 * Vpp * (1 + sinf(2 * PI * phase));
	// // Square
	// vdes = (p1 > 0.5) ? Vpp : 0;
	// // Triangle
	// vdes = (p1 < 0.5) ? 2 * p1 * Vpp : 2 * (1 - p1) * Vpp;

	float vact[2], iact[2];
	analogGetValues(vact, iact);

	// store the min and max
	vmax[0] = fmaxf(vmax[0], vact[0]);
	vmax[1] = fmaxf(vmax[1], vact[1]);
	vmin[0] = fminf(vmin[0], vact[0]);
	vmin[1] = fminf(vmin[1], vact[1]);
	// TODO: look at ADC for V1, V2

	voltageControl(vdes, vact[0], &htim3);
	voltageControl(vdes, vact[1], &htim4);
}
