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

typedef struct {
	int k;
	float a0;
	float a1[4];
	float A2[4 * 4];
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

#ifdef __cplusplus
}
#endif
