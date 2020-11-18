/**
 * @file loop_update.cpp
 * @author Avik De (avikde@gmail.com)
 * @brief Stupid single task update
 * @version 0.1
 * @date 2020-11-18
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include <uprightmpc2.h>
#include <main.h>
#include <stdio.h>
#include <stdint.h>

// Definitions
#define micros() (HAL_GetTick()*1000 + 1000 - SysTick->VAL*1000000/SystemCoreClock)
// Globals
UprightMPC_t up;

extern "C" void setup() {
	// Parameters
	float dt = 5;
	float g = 9.81e-3f;
	float TtoWmax = 3;
	float ws = 1e1;
	float wds = 1e3;
	float wpr = 1;
	float wvr = 1e3;
	float wpf = 5;
	float wvf = 2e3;
	float wthrust = 1e-1;
	float wmom = 1e-2;
	float Ib[3] = {3333,3333,1000};
	int maxIter = 10;
	// Call init
	umpcInit(&up, dt, g, TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom, Ib, maxIter);
}

extern "C" void loop() {
	uint32_t t1 = micros();

	// upri

	uint32_t t2 = micros();
	
	HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
	printf("hello %lu %lu\n", HAL_GetTick(), t2 - t1);
	HAL_Delay(100);
}

