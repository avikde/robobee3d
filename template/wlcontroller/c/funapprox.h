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
#include <matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	float a0;
	mxArray *a1;
	mxArray *A2;
} FunApprox_t;

void funApproxInit(FunApprox_t *fa, const float popts[/* 1 + k + k * (k + 1) / 2 */]);

float funApproxF(const FunApprox_t *fa, const mxArray *xi);

void funApproxDf(mxArray *Df, const FunApprox_t *fa, const mxArray *xi);

void funApproxClear(FunApprox_t *fa);

#ifdef __cplusplus
}
#endif
