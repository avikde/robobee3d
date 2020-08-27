/**
 * @file eigenc.h
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Dimensions
#define NU (4)
#define NW (6)

// index for col-major access
#define Cind(n, i, j) ((i) + (j)*n)

void matMult(float *C, const float *A, const float *B, const int m, const int n, const int k, const float alpha, int AT, int BT);

#ifdef __cplusplus
}
#endif

