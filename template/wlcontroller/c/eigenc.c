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
#ifdef BASIC_BLAS
#include <blas.h> // to test at home
#else
#include <cblas.h> // target computer has this
#endif

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
void matMult(float *C, const float *A, const float *B, const int m, const int n, const int k, const float alpha, int AT, int BT)
{
	// scalar values to use in sgemm
	const float zero = 0.0f;

#ifdef BASIC_BLAS
	// Source for sgemm http://www.netlib.org/clapack/cblas/sgemm.c
	char transa = 'N';
	char transb = 'N';
	if (AT == 1)
		transa = 'T';
	if (BT == 1)
		transb = 'T';
	
	ptrdiff_t um = m, un = n, uk = k;

	sgemm(&transa, &transb, 
		&um, &un, &uk, &alpha, 
		A, AT ? &uk : &um, 
		B, BT ? &un : &uk, 
		&zero, C, &um);
#else
	// Source for sgemm http://www.netlib.org/blas/cblas.h
	enum CBLAS_TRANSPOSE transa = CblasNoTrans;
	enum CBLAS_TRANSPOSE transb = CblasNoTrans;
	if (AT == 1)
		transa = CblasTrans;
	if (BT == 1)
		transb = CblasTrans;

	cblas_sgemm(CblasColMajor, transa, transb, 
		m, n, k, alpha, 
		A, AT ? k : m,
		B, BT ? n : k, 
		zero, C, m);
#endif
}
