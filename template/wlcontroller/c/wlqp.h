/**
 * @file wlqp.h
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

// Need to define these
void wrenchMap(float *w, const float *u);
void wrenchJacMap(float *dw_du, const float *u);

typedef struct {
	float umin[NU], umax[NU], U0[NU];
	float w0[NW];
	float Q[NW * NW];
} WLQP_t;

void wlqpInit(WLQP_t *wlqp);

void wlqpUpdate(WLQP_t *wlqp, float *u, const float *u0, const float *h0, const float *pdotdes);

#ifdef __cplusplus
}
#endif
