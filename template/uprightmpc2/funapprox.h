/**
 * @file funapprox.h
 * @author Avik De
 * @brief Scalar function approx
 * @version 0.1
 * @date 2020-08-12
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#define NDELU 4

typedef struct {
	int k;
	float a0;
	float a1[NDELU];
	float A2[NDELU * NDELU];
} FunApprox_t;

/**
 * @brief 
 * 
 * @param fa 
 * @param popts A2 components are row major!
 */
void funApproxInit(FunApprox_t *fa, const float popts[/* 1 + k + k * (k + 1) / 2 */]);

float funApproxF(const FunApprox_t *fa, const float *xi);

void funApproxDf(float *Df, const FunApprox_t *fa, const float *xi);

// 
typedef struct {
	float u0[NDELU], umin[NDELU], umax[NDELU], dumax[NDELU];
	float Qw[6*6];
  FunApprox_t fa[6];
} WLCon_t;

void wlConInit(WLCon_t *wl, const float u0[/* 4 */], const float umin[/* 4 */], const float umax[/* 4 */], const float dumax[/* 4 */], const float Qw[/* 6 */], float controlRate, const float popts[/* 90 */]);

void wlConUpdate(WLCon_t *wl, float u1[/* 4 */], float w0[/* 6 */], const float h0[/* 6 */], const float pdotdes[/* 6 */]);

// S function
void wlconS(float u1_y1[/* 4 */], float w0_y2[/* 6 */], const float u0init_u1[/* 4 */], const float umin_u2[/* 4 */], const float umax_u3[/* 4 */], const float dumax_u4[/* 4 */], const float Qw_u5[/* 6 */], float controlRate_u6, const float popts_u7[/* 90 */], const float h0_u8[/* 6 */], const float pdotdes_u9[/* 6 */]);

#ifdef __cplusplus
}
#endif
