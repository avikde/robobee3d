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

#ifndef PI
#define PI 3.1415926535897932384626433832795f
#endif
#ifndef HALF_PI
#define HALF_PI 1.5707963267948966192313216916398f
#endif
#ifndef TWO_PI
#define TWO_PI 6.283185307179586476925286766559f
#endif
#define DEG_TO_RAD 0.017453292519943295769236907684886f
#define RAD_TO_DEG 57.295779513082320876798154814105f

uint32_t millis();
uint32_t micros();
int osal_usleep(uint32_t usec);

#ifdef __cplusplus
}
#endif
