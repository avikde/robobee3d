/**
 * @file funapprox.c
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "funapprox.h"
#include <string.h>
#include <stdio.h>

void funApproxInit(FunApprox_t *fa, const float popts[/* 1 + k + k * (k + 1) / 2 */]) {
	const int k = 4;
	double *ptr;

	// unpack
	fa->a0 = popts[0];

	fa->a1 = mxCreateNumericMatrix(k, 1, mxSINGLE_CLASS, mxREAL);
	fa->A2 = mxCreateNumericMatrix(k, k, mxSINGLE_CLASS, mxREAL);
	
	ptr = mxGetPr(fa->a1);
	for (int i = 0; i < k; ++i) {
		printf("%d\n", ptr[i]);
	}
	// memcpy(ptr, &popts[1], k * sizeof(double));
	
	// // eigenUpperTriangularSet<4>(A2, p.segment</*  k *(k + 1) / 2  */10>(/* k + 1 */5));
	// // eigenUpperTriangularSet<4>(A2, p.segment<10>(5));
		
	// int kk = 0, N = 4;
	
	// // Upper triangular must be filled in row major to match python:
	// // https://github.com/avikde/robobee3d/pull/173#issuecomment-674082069
	// auto tvals = p.segment<k * (k + 1) / 2>(k + 1);
	// for (int i = 0; i < N; ++i) {
	// 	for (int j = i; j < N; ++j) {
	// 		A2(i, j) = tvals[kk];
	// 		kk++;
	// 	}
	// }
	// // Use A2.selfadjointView<Eigen::Upper>() after this
}
float funApproxF(const FunApprox_t *fa, const mxArray *xi) {
		// return a0 + xi.dot(a1) + 0.5 * xi.transpose() * A2.selfadjointView<Eigen::Upper>() * xi;
}

void funApproxDf(mxArray *Df, const FunApprox_t *fa, const mxArray *xi) {
		// return a1 + A2.selfadjointView<Eigen::Upper>() * xi;
}

void funApproxClear(FunApprox_t *fa) {
	mxDestroyArray(fa->a1);
	mxDestroyArray(fa->A2);
}
