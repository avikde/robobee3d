/**
 * @file wlcontroller.h
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once

#include "eigenc.h"

#ifdef __cplusplus
extern "C" {
#endif

void wlControllerInit();

void wlControllerUpdate(float *u, const float *u0init, const float *h0, const float *pdotdes);

#ifdef __cplusplus
}
#endif

