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

typedef struct {
	float umin[NU], umax[NU], U0[NU];
	float w0[NW], Qdiag[NW];
} WLQP_t;

void wlqpInit(WLQP_t *wlqp);

void wlqpUpdate(WLQP_t *wlqp, const float *u0, const float *p0, const float *h0, const float *pdes, const float *kpmom);


#ifdef __cplusplus
}
#endif
