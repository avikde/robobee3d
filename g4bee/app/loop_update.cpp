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

#define micros() (HAL_GetTick()*1000 + 1000 - SysTick->VAL*1000000/SystemCoreClock)

extern "C" void loopUpdate() {
	uint32_t t1 = micros();
	HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
	uint32_t t2 = micros();
	
	printf("hello %lu %lu\n", HAL_GetTick(), t2 - t1);
	HAL_Delay(100);
}

