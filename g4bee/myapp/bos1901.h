
#pragma once

#include "main.h"
#include "cmsis_os.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	SPI_HandleTypeDef *spi;
	// CS pin
	GPIO_TypeDef *GPIOx;
	uint16_t GPIO_Pin;
	// Data
	volatile uint8_t rxBuf[2];
	uint8_t txBuf[2];
} BOS1901;

void bos1901Init(BOS1901 *bos, SPI_HandleTypeDef *spi, GPIO_TypeDef *GPIOx, uint16_t GPIO_Pin);

void bos1901SetReadReg(BOS1901 *bos, uint8_t reg);

#ifdef __cplusplus
}
#endif

