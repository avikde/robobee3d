/**
 * @file wlqp.c
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "wlqp.h"
#include <math.h>
#include <osqp.h>
#include <workspace.h>
#include <string.h>

void wlqpInit(WLQP_t *wlqp) {
	int i, j;
	// populate default values
	for (i = 0; i < NW; ++i) {
		for (j = 0; j < NW; ++j) {
			wlqp->Q[Cind(NW, i, j)] = (i == j) ? (i < 3 ? 1.0f : 0.1f) : 0.0f;
		}
	}
	// Limits
	// u = Vmean, uoffs, udiff, h2
	wlqp->umin[0] = 90.0f;
	wlqp->umin[1] = -0.5f;
	wlqp->umin[2] = -0.2f;
	wlqp->umin[3] = -0.1f;

	wlqp->umax[0] = 160.0f;
	wlqp->umax[1] = 0.5f;
	wlqp->umax[2] = 0.2f;
	wlqp->umax[3] = 0.1f;
	// This is dumax
	wlqp->U0[0] = 5.0f;
	wlqp->U0[1] = 1e-2f;
	wlqp->U0[2] = 1e-2f;
	wlqp->U0[3] = 1e-2f;
	
	// Osqp init
	OSQPWorkspace *work = &workspace;
	osqp_set_default_settings(work->settings);
	osqp_update_max_iter(work, 40);
	// osqp_update_check_termination(work, 0); // don't check at all
}

static void wlqpSolve(float *du, const float *P, const float *q, const float *L, const float *U) {
	int i, j;

	static float Px_data[NU * (NU + 1) / 2];
	OSQPWorkspace *work = &workspace;

	// Get upper triangular
	int kk = 0;
	for (j = 0; j < NU; ++j) {
		for (i = 0; i <= j; ++i) {
			Px_data[kk] = P[Cind(NU, i, j)];
			kk++;
		}
	}

	// Update
	osqp_update_P(work, Px_data, OSQP_NULL, NU * (NU + 1) / 2);
	// osqp_update_A(work, A.data(), OSQP_NULL, work->data->m * work->data->n);
	osqp_update_lin_cost(work, q);
	osqp_update_bounds(work, L, U);

	/* int res = */ osqp_solve(work);

	memcpy(du, work->solution->x, NU * sizeof(float));

	// mexPrintf("du = %.3f %.3f %.3f, %s\n", du[0], du[1], du[2], work->info->status);
}

void wlqpUpdate(WLQP_t *wlqp, float *u, const float *u0, const float *h0, const float *pdotdes) {
	int i;
	static float A1[NW * NU];
	static float a0[NW];
	static float Q[NW * NW];
	static float dum[NW * NU];
	static float P[NU * NU];
	static float q[NU], L[NU], U[NU];

	// Sample numerical maps
	wrenchMap(wlqp->w0, u0);
	wrenchJacMap(A1, u0);
	
	// auto a0 = w0 - h0 - pdotdes;
	for (i = 0; i < NW; ++i) {
		a0[i] = wlqp->w0[i] - h0[i] - pdotdes[i];
	}

	// auto P = A1.transpose() * Qdiag.asDiagonal() * A1;
	matMult(dum, wlqp->Q, A1, NW, NU, NW, 1.0f, 0, 0);
	matMult(P, A1, dum, NU, NW, NW, 1.0f, 1, 0);

	// u_t q = A1.transpose() * Qdiag.cwiseProduct(a0);
	matMult(dum, wlqp->Q, a0, NW, 1, NW, 1.0f, 0, 0); // only using NW elements of dum
	matMult(q, A1, dum, NU, 1, NW, 1.0f, 1, 0);

	// u_t L = -this->U0;
	// u_t U = this->U0;
	for (i = 0; i < NU; ++i) {
		L[i] = -wlqp->U0[i];
		U[i] = wlqp->U0[i];
	}

	// Input limits (not just rate limits)
	for (i = 0; i < NU; ++i) {
		if (u0[i] < wlqp->umin[i])
			L[i] = 0; // do not reduce further
		else if (u0[i] > wlqp->umax[i])
			U[i] = 0; // do not increase further
	}

	// Solve
	wlqpSolve(u, P, q, L, U);
	for (i = 0; i < NU; ++i) {
		u[i] += u0[i];
	}
}

