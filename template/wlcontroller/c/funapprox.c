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
	const int k = FA_K;

	// unpack
	fa->a0 = popts[0];
	memcpy(fa->a1, &popts[1], k * sizeof(float));
	
  int kk = 0;
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i <= j; ++i) {
      fa->A2[Cind(k, i, j)] = fa->A2[Cind(k, j, i)] = popts[(k + 1) + kk];
      kk++;
    }
  }
}

float funApproxF(const FunApprox_t *fa, const float *xi) {
	static float vout[FA_K], fout;
	float res = fa->a0;
	// xi.dot(a1)
	matMult(vout, xi, fa->a1, 1, 1, FA_K, 1.0f, false, false);
	res += vout[0];
	// 0.5 * xi.transpose() * A2.selfadjointView<Eigen::Upper>() * xi
	// First A2 * xi
	matMult(vout, fa->A2, xi, FA_K, 1, FA_K, 1.0f, false, false);
	// Then 0.5 * xi.dot(vout)
	matMult(&fout, xi, vout, 1, 1, FA_K, 0.5f, false, false);
	return res + fout;
}

void funApproxDf(float *Df, const FunApprox_t *fa, const float *xi) {
	static float vout[FA_K];
	memcpy(Df, fa->a1, FA_K * sizeof(float));
	
	// First A2 * xi
	matMult(vout, fa->A2, xi, FA_K, 1, FA_K, 1.0f, false, false);
	for (int i = 0; i < FA_K; ++i)
		Df[i] += vout[i];
}
