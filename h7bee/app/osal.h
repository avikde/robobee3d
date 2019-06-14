/**
 * @file osal.h
 * @author Avik De (avikde@gmail.com)
 * @brief OS abstraction stuff
 * @version 0.1
 * @date 2019-06-14
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

uint32_t millis();
uint32_t micros();

#ifdef __cplusplus
}
#endif
