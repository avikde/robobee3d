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
#include "matmult.h"
#include <string.h>
#include <stdio.h>

#define NW 6

void funApproxInit(FunApprox_t *fa, const float popts[/* 1 + k + k * (k + 1) / 2 */]) {
	int i, j;
	fa->k = NDELU;

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
	static float vout[NDELU], fout;
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
	static float vout[NDELU];
	memcpy(Df, fa->a1, fa->k * sizeof(float));
	
	// First A2 * xi
	matMult(vout, fa->A2, xi, fa->k, 1, fa->k, 1.0f, 0, 0);
	for (i = 0; i < fa->k; ++i)
		Df[i] += vout[i];
}

// Wrench map
static void wrenchMap(const WLCon_t *wl, float *w, const float *u) {
	int i;
	for (i = 0; i < 6; ++i) {
		w[i] = funApproxF(&wl->fa[i], u);
	}
}

static void wrenchJacMap(const WLCon_t *wl, float *dw_du, const float *u) {
	int i, j;
	static float dwi_du[NDELU];
	for (i = 0; i < 6; ++i) {
		funApproxDf(dwi_du, &wl->fa[i], u);
		for (j = 0; j < NDELU; ++j) {
			// Copy into the col-major matrix
			dw_du[Cind(6, i, j)] = dwi_du[j]; // i >= 2 && i - 2 == j ? 1 : 0;  //
		}
	}
}

void wlConInit(WLCon_t *wl, const float u0[/* 4 */], const float umin[/* 4 */], const float umax[/* 4 */], const float dumax[/* 4 */], const float Qw[/* 6 */], float controlRate, const float popts[/* 90 */]) {
	int i;
	
	memset(wl->Qw, 0, 36 * sizeof(float));

	// Limits
	for (i = 0; i < NDELU; ++i) {
		wl->u0[i] = u0[i];
		wl->umin[i] = umin[i];
		wl->umax[i] = umax[i];
		wl->dumax[i] = dumax[i] / controlRate;
	}

	for (i = 0; i < 6; ++i) {
		funApproxInit(&wl->fa[i], &popts[15 * i]);
		wl->Qw[Cind(6, i, i)] = Qw[i];
	}
	// 
}

void wlConUpdate(WLCon_t *wl, float u1[/* 4 */], const float h0[/* 6 */], const float pdotdes[/* 6 */]) {
	int i;
	static float A1[NW * NDELU];
	static float a0[NW], delu[NDELU];
	static float dum[NW * NDELU];
	static float P[NDELU * NDELU];
	static float q[NDELU], L[NDELU], U[NDELU];

	// Sample numerical maps
	wrenchMap(wl, a0, wl->u0); // a0 = w(u0)
	wrenchJacMap(wl, A1, wl->u0);
	
	// auto a0 = w0 - h0 - pdotdes;
	for (i = 0; i < NW; ++i) {
		a0[i] += (-h0[i] - pdotdes[i]);
	}

	// u_t L = -this->U0;
	// u_t U = this->U0;
	for (i = 0; i < NDELU; ++i) {
		L[i] = -wl->dumax[i];
		U[i] = wl->dumax[i];
	}
	// Input limits (not just rate limits)
	for (i = 0; i < NDELU; ++i) {
		if (wl->u0[i] < wl->umin[i])
			L[i] = 0; // do not reduce further
		else if (wl->u0[i] > wl->umax[i])
			U[i] = 0; // do not increase further
	}

	// P = A1.transpose() * Qdiag.asDiagonal() * A1;
	matMult(dum, wl->Qw, A1, NW, NDELU, NW, 1.0f, 0, 0);
	matMult(P, A1, dum, NDELU, NW, NW, 1.0f, 1, 0);
	// q = -A1.transpose() * Qdiag.cwiseProduct(a0);
	matMult(dum, wl->Qw, a0, NW, 1, NW, 1.0f, 0, 0); // only using NW elements of dum
	matMult(q, A1, dum, NDELU, 1, NW, -1.0f, 1, 0);
	// Solve
	lsSolve(delu, P, q, NDELU, NDELU);

	// Limit
	for (i = 0; i < NDELU; ++i) {
		if (delu[i] < L[i])
			delu[i] = L[i];
		else if (delu[i] > U[i])
			delu[i] = U[i];
	}
	
	// Update result
	for (i = 0; i < NDELU; ++i) {
		wl->u0[i] += delu[i];
		u1[i] = wl->u0[i];
	}
}
