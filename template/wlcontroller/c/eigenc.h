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

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Dimensions
#define NU (4)
#define NW (6)

// index for col-major access
#define Cind(n, i, j) ((i) + (j)*n)

void matMult(float *C, const float *A, const float *B, int64_t m, int64_t n, int64_t k, float alpha, bool AT, bool BT);

#ifdef __cplusplus
}
#endif

