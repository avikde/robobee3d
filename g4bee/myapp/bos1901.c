
#include "bos1901.h"
#include <string.h>

// This is a helper and the only HAL-dependent part
static void _bos1901rw(BOS1901 *bos)
{
	HAL_GPIO_WritePin(bos->GPIOx, bos->GPIO_Pin, GPIO_PIN_RESET);
	HAL_SPI_TransmitReceive(bos->spi, bos->txBuf, bos->rxBuf, 2, 100);
	HAL_GPIO_WritePin(bos->GPIOx, bos->GPIO_Pin, GPIO_PIN_SET);
}

void bos1901Init(BOS1901 *bos, SPI_HandleTypeDef *spi, GPIO_TypeDef *GPIOx, uint16_t GPIO_Pin)
{
	bos->spi = spi;
	bos->GPIOx = GPIOx;
	bos->GPIO_Pin = GPIO_Pin;
	memset(bos->rxBuf, 0, sizeof(bos->rxBuf));
	memset(bos->txBuf, 0, sizeof(bos->txBuf));

	// HAL_Delay(100);
	bos1901SetReadReg(bos, 0x4);
}

void bos1901SetReadReg(BOS1901 *bos, uint8_t reg)
{
	// first 4 bits are address
	bos->txBuf[0] = (0x5) << 4 | (reg >> 1); // set SDO to have 0x4
	bos->txBuf[1] = (reg & 0b1) << 7;
	_bos1901rw(bos);
}

