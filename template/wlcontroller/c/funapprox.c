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
	int i, j;
    fa->k = NU;

	// unpack
	fa->a0 = popts[0];
	memcpy(fa->a1, &popts[1], fa->k * sizeof(float));
	
	// Unpack A2 components row major, fill out upper triangular A2
  int kk = 0;
  for (i = 0; i < fa->k; ++i) {
    for (j = i; j < fa->k; ++j) {
      fa->A2[Cind(fa->k, i, j)] = fa->A2[Cind(fa->k, j, i)] = popts[(fa->k + 1) + kk];
      kk++;
    }
  }
}

float funApproxF(const FunApprox_t *fa, const float *xi) {
	static float vout[NU], fout;
	float res = fa->a0;
	// xi.dot(a1)
	matMult(vout, xi, fa->a1, 1, 1, fa->k, 1.0f, 0, 0);
	res += vout[0];
	// 0.5 * xi.transpose() * A2.selfadjointView<Eigen::Upper>() * xi
	// First A2 * xi
	matMult(vout, fa->A2, xi, fa->k, 1, fa->k, 1.0f, 0, 0);
	// Then 0.5 * xi.dot(vout)
	matMult(&fout, xi, vout, 1, 1, fa->k, 0.5f, 0, 0);
	return res + fout;
}

void funApproxDf(float *Df, const FunApprox_t *fa, const float *xi) {
    int i;
	static float vout[NU];
	memcpy(Df, fa->a1, fa->k * sizeof(float));
	
	// First A2 * xi
	matMult(vout, fa->A2, xi, fa->k, 1, fa->k, 1.0f, 0, 0);
	for (i = 0; i < fa->k; ++i)
		Df[i] += vout[i];
}
