/**
 * @file eigenc.c
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "eigenc.h"
// Use MATLAB -> BLAS https://www.mathworks.com/help/matlab/matlab_external/calling-lapack-and-blas-functions-from-mex-files.html
#include <matrix.h>
#include <blas.h>

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
void matMult(float *C, const float *A, const float *B, int64_t m, int64_t n, int64_t k, float alpha, bool AT, bool BT)
{
	const char *chn = "N";
	const char *cht = "T";
	// scalar values to use in sgemm
	float zero = 0.0f;

	// Source for sgemm http://www.netlib.org/clapack/cblas/sgemm.c
	const char *transa = chn;
	const char *transb = chn;
	if (AT)
		transa = cht;
	if (BT)
		transb = cht;

	sgemm(transa, transb, &m, &n, &k, &alpha, A, AT ? &k : &m, B, BT ? &n : &k, &zero, C, &m);
}
