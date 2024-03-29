/**
 * @file InvensenseIMU.h
 * @author Avik De (avikde@gmail.com)
 * @brief Driver for ICM-20649
 * @version 0.1
 * @date 2019-06-14
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once

#include <stdint.h>
#include <stm32h7xx_hal.h>

#ifdef __cplusplus
extern "C" {
#endif

// MPU 6000 registers
#define MPUREG_WHOAMI 0x75
#define MPUREG_SMPLRT_DIV 0x19
#define MPUREG_CONFIG 0x1A
#define MPUREG_GYRO_CONFIG 0x1B
#define MPUREG_ACCEL_CONFIG 0x1C
#define MPUREG_FIFO_EN 0x23
#define MPUREG_INT_PIN_CFG 0x37
#define MPUREG_INT_ENABLE 0x38
#define MPUREG_INT_STATUS 0x3A
#define MPUREG_ACCEL_XOUT_H 0x3B
#define MPUREG_ACCEL_XOUT_L 0x3C
#define MPUREG_ACCEL_YOUT_H 0x3D
#define MPUREG_ACCEL_YOUT_L 0x3E
#define MPUREG_ACCEL_ZOUT_H 0x3F
#define MPUREG_ACCEL_ZOUT_L 0x40
#define MPUREG_TEMP_OUT_H 0x41
#define MPUREG_TEMP_OUT_L 0x42
#define MPUREG_GYRO_XOUT_H 0x43
#define MPUREG_GYRO_XOUT_L 0x44
#define MPUREG_GYRO_YOUT_H 0x45
#define MPUREG_GYRO_YOUT_L 0x46
#define MPUREG_GYRO_ZOUT_H 0x47
#define MPUREG_GYRO_ZOUT_L 0x48
#define MPUREG_USER_CTRL 0x6A
#define MPUREG_PWR_MGMT_1 0x6B
#define MPUREG_PWR_MGMT_2 0x6C
#define MPUREG_FIFO_COUNTH 0x72
#define MPUREG_FIFO_COUNTL 0x73
#define MPUREG_FIFO_R_W 0x74

// Configuration bits
#define BIT_SLEEP 0x40
#define BIT_H_RESET 0x80
#define BITS_CLKSEL 0x07
#define MPU_CLK_SEL_PLLGYROX 0x01
#define MPU_CLK_SEL_PLLGYROZ 0x03
#define MPU_EXT_SYNC_GYROX 0x02
#define BITS_FS_250DPS 0x00
#define BITS_FS_500DPS 0x08
#define BITS_FS_1000DPS 0x10
#define BITS_FS_2000DPS 0x18
#define BITS_FS_MASK 0x18
#define BITS_DLPF_CFG_256HZ_NOLPF2 0x00
#define BITS_DLPF_CFG_188HZ 0x01
#define BITS_DLPF_CFG_98HZ 0x02
#define BITS_DLPF_CFG_42HZ 0x03
#define BITS_DLPF_CFG_20HZ 0x04
#define BITS_DLPF_CFG_10HZ 0x05
#define BITS_DLPF_CFG_5HZ 0x06
#define BITS_DLPF_CFG_2100HZ_NOLPF 0x07
#define BITS_DLPF_CFG_MASK 0x07
#define BIT_INT_ANYRD_2CLEAR 0x10
#define BIT_RAW_RDY_EN 0x01
#define BIT_I2C_IF_DIS 0x10
#define BIT_INT_STATUS_DATA 0x01

// ICM 20649
#define ICMREG_BANK_SEL 0x7f // same address in banks 0--3
#define ICMREG_WHOAMI 0x80
#define ICMREG_USER_CTRL 0x3
#define ICMREG_PWR_MGMT_1 0x6
#define ICMREG_PWR_MGMT_2 0x7
#define ICMREG_GYRO_CONFIG_1 0x1
#define ICMREG_GYRO_CONFIG_2 0x2
#define ICMREG_ACCEL_CONFIG 0x14
#define ICMREG_ACCEL_CONFIG_2 0x15
#define ICMREG_ACCEL_XOUT_H 0x2d

#pragma pack(push, 1)

typedef struct
{
	int16_t x, y, z;
} tAxis;

typedef union _MPUData {
	uint8_t rawByte[14];
	struct
	{
		tAxis accel;
		tAxis gyro;
		int16_t temperature;
	} dataICM;
} MPUData;

#pragma pack(pop)

#define INVENSENSE_BUF_SIZE 16

typedef struct _InvensenseIMU {
	SPI_HandleTypeDef *_SPI;
	uint8_t accelFS, gyroFS;
	uint8_t txBuf[INVENSENSE_BUF_SIZE], rxBuf[INVENSENSE_BUF_SIZE];
	float acc[3], gyr[3];
	float temperature;
} InvensenseIMU;

int invensenseIMUInit(InvensenseIMU *imu, SPI_HandleTypeDef *_SPI);

void invensenseIMUUpdate(InvensenseIMU *imu);

#ifdef __cplusplus
}
#endif
