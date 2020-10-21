/**
 * @file uprightmpc2.c
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-10-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "uprightmpc2.h"
#include "matmult.h"
#include "osqp.h"
#include "workspace.h"
#include <stdio.h>
#include <string.h>

void umpcInit(UprightMPC_t *up, float dt, float g, const float smin[/* 3 */], const float smax[/* 3 */], float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, int maxIter) {
	up->dt = dt;
	up->g = dt;
	up->Tmax = TtoWmax * g;
	// Weights
	for (int i = 0; i < 3; ++i) {
		up->Qyr[i] = wpr;
		up->Qyf[i] = wpf;
		up->Qyr[3 + i] = up->Qyf[3 + i] = ws;
		up->Qdyr[i] = wvr;
		up->Qdyf[i] = wvf;
		up->Qdyr[3 + i] = up->Qdyf[3 + i] = wds;
	}
	up->R[0] = wthrust;
	up->R[1] = up->R[2] = wmom;
	// Limits
	for (int i = 0; i < 3; ++i) {
		up->smin[i] = smin[i];
		up->smax[i] = smax[i];
	}

	// OSQP
	osqp_update_max_iter(&workspace, maxIter);
	osqp_update_check_termination(&workspace, 0);
}

static void umpcUpdate1(UprightMPC_t *up, float T0sp, const float s0s[/* 3*N */], const float Btaus[/*  */], const float y0[/* 6 */], const float dy0[/* 6 */]) {

}

static void umpcUpdateConstraint(UprightMPC_t *up, float T0, const float s0s[/*  */], const float Btaus[/*  */], const float y0[/*  */], const float dy0[/*  */]) {
	static float y1[UMPC_NY];
	// Update vector
	memset(up->l, 0, UMPC_NX * sizeof(float));
	// y1 = y0 + dt * dy0
	for (int i = 0; i < UMPC_NY; ++i) {
		y1[i] = y0[i] + up->dt * dy0[i];
		// l[:ny] = -y1
		up->l[i] = -y1[i];
	}
	
	for (int k = 0; k < UMPC_N; ++k) {
		if (k == 0) {

		}
	}
	
	// copy for dynamics
	memcpy(up->u, up->l, UMPC_NC * sizeof(float));
}

void umpcUpdate(UprightMPC_t *up, float uquad[/* 3 */], float accdes[/* 6 */], const float p0[/* 6 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */]) {
	// TODO:
	// for (int i = 0; i < 9; ++i) {
	// 	printf("%.2f ", R0[i]);
	// }
	// printf("\n");
	matMult(uquad, R0, p0, 3, 1, 3, 1.0f, 1, 0);
}
