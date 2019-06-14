/**
 * @file ImuTask.cpp
 * @author Avik De (avikde@gmail.com)
 * @brief Task for flap to flap control
 * @version 0.1
 * @date 2019-06-14
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "FreeRTOS.h"
#include "task.h"
#include "main.h"
#include "cmsis_os.h"
#include "InvensenseIMU.h"
// TODO: include Eigen

// defined in main.c
extern "C" SPI_HandleTypeDef hspi2;
extern "C" {
// TODO: RTOS message passing
float extest, eytest, eztest;
}

extern "C" void startImuTask(void const *argument)
{
	InvensenseIMU imu;
	int invres = invensenseIMUInit(&imu, &hspi2);
	// TODO: EKF with Eigen

	for (;;)
	{
		if (invres >= 0)
		{
			invensenseIMUUpdate(&imu);
			// TODO: EKF
			extest = imu.acc[0];
			eytest = imu.acc[1];
			eztest = imu.acc[2];
		}
		osDelay(1);
	}
}
