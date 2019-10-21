
#include "bos1901.h"
#include <string.h>

// This is a helper and the only HAL-dependent part
static void _bos1901transfer(BOS1901 *bos)
{
	HAL_GPIO_WritePin(bos->GPIOx, bos->GPIO_Pin, GPIO_PIN_RESET);
	HAL_SPI_TransmitReceive(bos->spi, bos->txBuf, bos->rxBuf, 2, 100);
	HAL_GPIO_WritePin(bos->GPIOx, bos->GPIO_Pin, GPIO_PIN_SET);
}


static uint16_t _bos1901ParseSDO(BOS1901 *bos)
{
	// Interpret the SDO data according to sec 6.4.2
	if (bos->BC <= 0x0e)
	{
		bos->rxBuf[0] &= 0x0f;
	}
	else if (bos->BC == 0x10)
	{
		memset(bos->rxBuf, 0, sizeof(bos->rxBuf));
	}
	// Swap for endianness
	uint16_t ret = (bos->rxBuf[0] << 8) | bos->rxBuf[1];
	return ret;
}

void bos1901Init(BOS1901 *bos, SPI_HandleTypeDef *spi, GPIO_TypeDef *GPIOx, uint16_t GPIO_Pin)
{
	bos->spi = spi;
	bos->GPIOx = GPIOx;
	bos->GPIO_Pin = GPIO_Pin;
	memset(bos->rxBuf, 0, sizeof(bos->rxBuf));
	memset(bos->txBuf, 0, sizeof(bos->txBuf));

	// HAL_Delay(100);
}

uint16_t bos1901Config(BOS1901 *bos, uint8_t BC, uint8_t OE, uint8_t PLAY)
{
	// first 4 bits are address
	bos->txBuf[0] = (0x5) << 4;
	// See datasheet table 16
	bos->BC = BC;
	bos->txBuf[0] |= (BC >> 1);
	bos->txBuf[1] = ((BC & 0b1) << 7) | ((OE & 0b1) << 4) | (PLAY & 0b111);
	_bos1901transfer(bos);
	return _bos1901ParseSDO(bos);
}

uint16_t bos1901rw(BOS1901 *bos, uint8_t addr, uint16_t data)
{
	// first 4 bits are address
	bos->txBuf[0] = (addr & 0xf) << 4 | ((data >> 8) & 0xf);
	bos->txBuf[1] = data & 0xff;
	_bos1901transfer(bos);
	return _bos1901ParseSDO(bos);
}

#define V_TO_CT(x) (100*x)
#define NFIFO (64)
#define DPHASE (1/((float)NFIFO))
void bos1901AddWave(BOS1901 *bos, WaveFunc f)
{
	// float ct = V_TO_CT(V);
	for (int i = 0; i < NFIFO; ++i)
	{
		float ct = V_TO_CT(f(DPHASE * i)); // get the wave function at this phase
		if (ct < 0)
		{
			// 2's complement 12-bit TODO:
		}
		bos1901rw(bos, 0, (uint16_t)ct);
	}
}