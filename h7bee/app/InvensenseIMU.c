/**
 * @file InvensenseIMU.c
 * @author Avik De (avikde@gmail.com)
 * @brief Driver for ICM-20649
 * @version 0.1
 * @date 2019-06-14
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "InvensenseIMU.h"
#include "main.h" // for pin/port
#include <string.h>
#include "osal.h"

static uint8_t *_invread(InvensenseIMU *imu, uint8_t reg, uint16_t length)
{
	imu->txBuf[0] = reg | 0x80;
	memset(&imu->txBuf[1], 0, length);
	HAL_GPIO_WritePin(SSIMU_GPIO_Port, SSIMU_Pin, GPIO_PIN_RESET);
	HAL_SPI_TransmitReceive(imu->_SPI, imu->txBuf, imu->rxBuf, length + 1, 100);
	HAL_GPIO_WritePin(SSIMU_GPIO_Port, SSIMU_Pin, GPIO_PIN_SET);
	return &imu->rxBuf[1];
}

static void _invwrite(InvensenseIMU *imu, uint8_t reg, uint8_t data)
{
	imu->txBuf[0] = reg & 0x7f;
	imu->txBuf[1] = data;
	HAL_GPIO_WritePin(SSIMU_GPIO_Port, SSIMU_Pin, GPIO_PIN_RESET);
	HAL_SPI_Transmit(imu->_SPI, imu->txBuf, 2, 100);
	HAL_GPIO_WritePin(SSIMU_GPIO_Port, SSIMU_Pin, GPIO_PIN_SET);
}

static void _invswapData(uint8_t *data, uint16_t datalen)
{
	datalen /= 2;
	uint8_t t;
	while (datalen--)
	{
		t = data[0];
		data[0] = data[1];
		data[1] = t;
		data += 2;
	}
}

int invensenseIMUInit(InvensenseIMU *imu, SPI_HandleTypeDef *_SPI)
{
	imu->_SPI = _SPI;
	// For ICM, use these as to set FS_SEL before calling init()
	imu->gyroFS = 0b0;   // 00 = 500dps, 01 = 1000 dps, 10 = 2000 dps, 11 = 4000 dps
	imu->accelFS = 0b10; // 00 = 4g, 01 = 8g, 10 = 16g, 11 = 30g

	// Select bank 0 for config stuff
	_invwrite(imu, ICMREG_BANK_SEL, 0 << 4);

	// Start up in reg bank 0
	int ndetectretries = 10;
	int res = -1; // not detected yet
	while (res < 0 && ndetectretries > 0)
	{
		uint8_t *magic = _invread(imu, ICMREG_WHOAMI, 1);
		if (magic[0] == 0xe1)
			res = 0;
		ndetectretries--;
		osal_usleep(10000);
	}
	if (ndetectretries == 0)
		return -1;

	// Chip reset
	_invwrite(imu, ICMREG_PWR_MGMT_1, BIT_H_RESET); // writing 1 to bit 7 will reset device
	osal_usleep(100000); // Startup time delay

	// write USER_CTRL
	_invwrite(imu, ICMREG_USER_CTRL, BIT_I2C_IF_DIS);
	// Read PWR_MGMT_1 here: says 0x41 (sleep mode; need to wake up)
	_invwrite(imu, ICMREG_PWR_MGMT_1, 0x1); // auto clock
	// volatile uint8_t pm = read(ICMREG_PWR_MGMT_1);

	osal_usleep(10000); // Writing to config registers wasn't working (wait after wake up?)
	// default sample rates are the highest for both gyro and accel (see 8.3.1, 8.3.11-12)

	// the following are in reg bank 2
	_invwrite(imu, ICMREG_BANK_SEL, 2 << 4);
	// 8.3.2: DLPF=0, FCHOICE=1 (enable DLPF)
	_invwrite(imu, ICMREG_GYRO_CONFIG_1, (0 << 3) | (imu->gyroFS << 1) | 0b1);
	// gyro config 2: use defaults
	// 8.3.15: DLPF=0, FCHOICE=1
	_invwrite(imu, ICMREG_ACCEL_CONFIG, (0 << 3) | (imu->accelFS << 1) | 0b1);

	// volatile uint8_t gc = read(ICMREG_GYRO_CONFIG_1);
	// volatile uint8_t ac = read(ICMREG_ACCEL_CONFIG);

	// Select bank 0 for reading data again
	_invwrite(imu, ICMREG_BANK_SEL, 0 << 4);

	return 0;
}

// datasheet lists LSB/dps => 2*pi/(180x) is rad/s/LSB
const static float GYRO_SENS[4] = {0.00053292496, 0.00106422515, 0.0021284503, 0.00425690061};
// LSB/g => 9.81/x = m/s^2/LSB
const static float ACCEL_SENS[4] = {0.00119750976, 0.00239501953, 0.00479003906, 0.00958007812};
const static float TEMP_SENS = 0.00299517776;

void invensenseIMUUpdate(InvensenseIMU *imu)
{
	// Added some empirical tuning to match VN100
	float accSens = ACCEL_SENS[imu->accelFS];
	float gyrSens = GYRO_SENS[imu->gyroFS];
	uint8_t dataReg = ICMREG_ACCEL_XOUT_H;
	// Read
	MPUData *data = (MPUData *)_invread(imu, dataReg, sizeof(MPUData));
	// Data is big endian; swap in place
	_invswapData(data->rawByte, sizeof(MPUData));

	imu->acc[0] = accSens * data->dataICM.accel.x;
	imu->acc[1] = accSens * data->dataICM.accel.y;
	imu->acc[2] = accSens * data->dataICM.accel.z;
	imu->gyr[0] = gyrSens * data->dataICM.gyro.x;
	imu->gyr[1] = gyrSens * data->dataICM.gyro.y;
	imu->gyr[2] = gyrSens * data->dataICM.gyro.z;
	imu->temperature = 21 + TEMP_SENS * data->dataICM.temperature;
}
