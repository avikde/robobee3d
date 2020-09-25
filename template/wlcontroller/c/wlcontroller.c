/**
 * @file wlcontroller.c
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "wlcontroller.h"
#include "funapprox.h"
#include "wlqp.h"
// #include <mex.h>

// Global vars needed
FunApprox_t fa[NW];
WLQP_t wlqp;
float g_u0[NU] = {-1, 0, 0, 0}; // state

static void wlControllerInit(const float *u0init, const float *popts, float controlRate) {
	// mexPrintf("HI INITING\n\n");
	int i;

	wlqpInit(&wlqp, controlRate);

	for (i = 0; i < NW; ++i) {
		funApproxInit(&fa[i], &popts[15 * i]);
	}
	for (i = 0; i < NU; ++i) {
		g_u0[i] = u0init[i];
	}
}

void wlControllerUpdate(float *u, const float *u0init, const float *h0, const float *pdotdes, const float *popts, float controlRate) {
	int i;

	if (g_u0[0] < 0) {
		wlControllerInit(u0init, popts, controlRate);
	}
	wlqpUpdate(&wlqp, u, g_u0, h0, pdotdes);
	for (i = 0; i < NU; ++i) {
		g_u0[i] = u[i];
	}
}

// Wrench map (private) ---------

void wrenchMap(float *w, const float *u) {
	int i;
	for (i = 0; i < NW; ++i) {
		w[i] = funApproxF(&fa[i], u);
	}
}

void wrenchJacMap(float *dw_du, const float *u) {
	int i, j;
	static float dwi_du[NU];
	for (i = 0; i < NW; ++i) {
		funApproxDf(dwi_du, &fa[i], u);
		for (j = 0; j < NU; ++j) {
			// Copy into the col-major matrix
			dw_du[Cind(NW, i, j)] = dwi_du[j];
		}
	}
}
