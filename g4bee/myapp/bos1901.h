
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
	uint8_t rxBuf[2];
	uint8_t txBuf[2];
	// which reg is being broadcast
	uint8_t sdoReg;
} BOS1901;

void bos1901Init(BOS1901 *bos, SPI_HandleTypeDef *spi, GPIO_TypeDef *GPIOx, uint16_t GPIO_Pin);

/**
 * @brief See datasheet 6.4.2
 * 
 * @param bos 
 * @param reg 
 */
void bos1901SetSDOBroadcast(BOS1901 *bos, uint8_t reg);

uint16_t bos1901rw(BOS1901 *bos, uint8_t addr, uint16_t data);

#ifdef __cplusplus
}
#endif

