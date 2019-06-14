
#include "FreeRTOS.h"
#include "task.h"
#include "main.h"
#include "cmsis_os.h"
#include "InvensenseIMU.h"

// defined in main.c
extern SPI_HandleTypeDef hspi2;
// TODO: RTOS message passing
float extest, eytest, eztest;

void startImuTask(void const *argument)
{
	InvensenseIMU imu;
	invensenseIMUInit(&imu, &hspi2);

	for (;;)
	{
		invensenseIMUUpdate(&imu);
		extest = imu.acc[0];
		eytest = imu.acc[1];
		eztest = imu.acc[2];
		osDelay(1);
	}
}
