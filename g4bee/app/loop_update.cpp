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

extern "C" void loopUpdate() {
	HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
	volatile uint32_t millis = HAL_GetTick();
	printf("hello %d\n", millis);
	HAL_Delay(1000);
}

