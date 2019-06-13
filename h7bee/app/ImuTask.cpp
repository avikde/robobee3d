
#include "FreeRTOS.h"
#include "task.h"
#include "main.h"
#include "cmsis_os.h"

// defined in main.c
extern "C" SPI_HandleTypeDef hspi2;

extern "C" {
// TODO: RTOS message passing
float extest, eytest, eztest;
}

extern "C" void startImuTask(void const *argument)
{
	// InvensenseIMU imu(hspi2);
	// imu.type = InvensenseIMU::Invensense_ICM20649;
	// // Init
	// imu.init(0);

	for (;;)
	{
			// imu.update();
			// extest = imu.acc[0];
			// eytest = imu.acc[1];
			// eztest = imu.acc[2];
		osDelay(1);
	}
}
