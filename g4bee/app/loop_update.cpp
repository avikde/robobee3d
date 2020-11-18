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
#include <ematmult.hpp>
#include <main.h>
#include <stdio.h>
#include <stdint.h>

// Definitions
extern "C" volatile uint32_t _millis;
// #define micros() (HAL_GetTick()*1000 + 1000 - SysTick->VAL*1000000/SystemCoreClock)
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
	int maxIter = 20;
	// Call init
	umpcInit(&up, dt, g, TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom, Ib, maxIter);
}

extern "C" void loop() {
	// States
	Eigen::Vector3f p0 = Eigen::Vector3f::Zero();
	Eigen::Matrix3f R0 = Eigen::Matrix3f::Identity();
	Vec6_t dq0 = Vec6_t::Zero();
	Eigen::Vector3f pdes = Eigen::Vector3f(0, 0, 10);
	Eigen::Vector3f dpdes = Eigen::Vector3f(0, 0, 0.1);
	Eigen::Vector3f sdes = Eigen::Vector3f::UnitZ();

	// Outputs
	Eigen::Vector3f uquad;
	Vec6_t accdes;

	// uint32_t t1 = micros();

	umpcUpdate(&up, uquad.data(), accdes.data(), p0.data(), R0.data(), dq0.data(), pdes.data(), dpdes.data(), sdes.data());

	// uint32_t t2 = micros();
	
	HAL_GPIO_TogglePin(LD2_GPIO_Port, LD2_Pin);
	printf("hello %lu %lu\n", HAL_GetTick(), _millis);
	HAL_Delay(100);
}

