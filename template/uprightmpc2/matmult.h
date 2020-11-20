/**
 * @file matmult.h
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

// index for col-major access
#define Cind(n, i, j) ((i) + (j)*n)

/**
 * @brief Multiply matrices
 * 
 * @param C Output C = alpha * op(A) * op(B)
 * @param A Column major (can use AT to switch)
 * @param B Column major (can use BT to switch)
 * @param m Number of rows of op(A), C
 * @param n Number of columns of op(B), C
 * @param k Number of columns of op(A), rows of op(B)
 * @param alpha
 * @param AT transpose A?
 * @param BT transpose B?
 */
void matMult(float *C, const float *A, const float *B, const int m, const int n, const int k, const float alpha, int AT, int BT);

/**
 * @brief Solve A*x=b
 * 
 * @param x 
 * @param A 
 * @param b 
 * @param m rows of A
 * @param n cols of A
 */
void lsSolve(float *x, const float *A, const float *b, const int m, const int n);

#ifdef __cplusplus
}
#endif

