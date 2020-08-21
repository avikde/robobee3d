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

#ifdef __cplusplus
extern "C" {
#endif

// index for col-major access
#define Cind(n, i, j) ((i) + (j)*n)

void matMult(float *C, const float *A, const float *B, size_t m, size_t n, size_t k, float alpha, bool AT, bool BT);

#ifdef __cplusplus
}
#endif

