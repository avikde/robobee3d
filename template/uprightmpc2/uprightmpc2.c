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
#include "funapprox.h"
#include "matmult.h"
#include "osqp.h"
#include "workspace.h"
#include <stdio.h>
#include <string.h>

void umpcInit(UprightMPC_t *up, float dt, float g, float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, const float Ib[/* 3 */], int maxIter) {
	static float Ibi[9];
	int i, k, n1, n2, offs;

	up->dt = dt;
	up->g = g;
	up->Tmax = TtoWmax * g;
	// Weights
	for (i = 0; i < 3; ++i) {
		up->Qyr[i] = wpr;
		up->Qyf[i] = wpf;
		up->Qyr[3 + i] = up->Qyf[3 + i] = ws;
		up->Qdyr[i] = wvr;
		up->Qdyf[i] = wvf;
		up->Qdyr[3 + i] = up->Qdyf[3 + i] = wds;
	}
	up->R[0] = wthrust;
	up->R[1] = up->R[2] = wmom;

	// Other updates
	for (i = 0; i < UMPC_NY; ++i) {
		up->c0[i] = i == 2 ? -up->g : 0;
	}
	up->T0 = 0;

	// skew(e3) entered in col major order
	memset(up->e3h, 0, 9 * sizeof(float));
	up->e3h[1] = 1;
	up->e3h[3] = -1;

	// e3h*Ibi calculation
	for (i = 0; i < 9; ++i) {
		Ibi[i] = 0;
	}
	Ibi[0] = 1.0f / Ib[0];
	Ibi[4] = 1.0f / Ib[1];
	Ibi[8] = 1.0f / Ib[2];
	matMult(up->e3hIbi, up->e3h, Ibi, 3, 3, 3, 1.0f, 0, 0); // e3h*Ibi

	for (i = 0; i < UMPC_NX; ++i) {
		up->q[i] = 0;
	}
	for (i = 0; i < UMPC_NC; ++i) {
		up->l[i] = up->u[i] = 0;
	}

	// Ax idx ---
	// Left third, T0*dt
	offs = 0;
	n2 = 2*UMPC_NY + 3; // nnz in each block col on the left
	for (k = 0; k < UMPC_N-2; ++k) {
		up->Ax_idx[offs+0] = n2*k + 8;
		up->Ax_idx[offs+1] = n2*k + 11;
		up->Ax_idx[offs+2] = n2*k + 14;
		offs += 3;
	}
	up->nAxT0dt = offs;

	// Middle third, dt
	n1 = (2*UMPC_N-1)*UMPC_NY + (UMPC_N-2)*3; // All the nnz in the left third
	n2 = 3*UMPC_NY; // nnz in each of the first N-1 block cols in the middle third
	for (k = 0; k < UMPC_N; ++k) {
		up->Ax_idx[offs+0] = n1 + n2*k + 0;
		if (k < UMPC_N-1) {
			up->Ax_idx[offs+1] = n1 + n2*k + 3;
			up->Ax_idx[offs+2] = n1 + n2*k + 6;
			up->Ax_idx[offs+3] = n1 + n2*k + 9;
			up->Ax_idx[offs+4] = n1 + n2*k + 12;
			up->Ax_idx[offs+5] = n1 + n2*k + 15;
		} else {
			up->Ax_idx[offs+1] = n1 + n2*k + 2;
			up->Ax_idx[offs+2] = n1 + n2*k + 4;
			up->Ax_idx[offs+3] = n1 + n2*k + 6;
			up->Ax_idx[offs+4] = n1 + n2*k + 8;
			up->Ax_idx[offs+5] = n1 + n2*k + 10;
		}
		offs += 6;
	}
	up->nAxdt = offs;

	// Right third
	n1 += 3*UMPC_NY*(UMPC_N-1) + 2*UMPC_NY; // All the nnz in the left, middle third
	n2 = 10; // nnz in each B0 + 1 for thrust lim
	// s0
	for (k = 0; k < UMPC_N; ++k) {
		for (i = 0; i < 3; ++i)
			up->Ax_idx[offs+i] = n1 + n2*k + i;
		offs += 3;
	}
	// Btau
	for (k = 0; k < UMPC_N; ++k) {
		for (i = 0; i < 6; ++i)
			up->Ax_idx[offs+i] = n1 + n2*k + 4 + i;
		offs += 6;
	}

	// OSQP ---
	osqp_update_max_iter(&workspace, maxIter);
	osqp_update_check_termination(&workspace, 0);
}

// Return the ith element of (A0*y), where A0 = (dt*N)
static float A0_times_i(const UprightMPC_t *up, const float y[/* ny */], int i) {
	return (i < 3) ? up->T0 * y[i + 3] : 0;
}

// Same s0, Btau for each k (that's what ended up happening in python anyway)
static void umpcUpdateConstraint(UprightMPC_t *up, const float s0[/*  */], const float Btau[/*  */], const float y0[/*  */], const float dy0[/*  */]) {
	static float y1[UMPC_NY];
	int i, k, offs;
	
	// Update vector ---
	for (i = 0; i < UMPC_N * UMPC_NY; ++i) {
		up->l[i] = 0;
	}
	// y1 = y0 + dt * dy0
	for (i = 0; i < UMPC_NY; ++i) {
		y1[i] = y0[i] + up->dt * dy0[i];
		// l[:ny] = -y1
		up->l[i] = -y1[i];
	}
	// Dynamics constraints
	for (k = 0; k < UMPC_N; ++k) {
		float *pk = &up->l[UMPC_NY * (UMPC_N + k)];
		for (i = 0; i < UMPC_NY; ++i) {
			if (k == 0) {
				pk[i] = -dy0[i] - up->dt * A0_times_i(up, y0, i) - up->dt * up->c0[i];
			} else if (k == 1) {
				pk[i] = -up->dt * A0_times_i(up, y1, i) - up->dt * up->c0[i];
			} else {
				pk[i] = -up->dt * up->c0[i];
			}
		}
	}
	// copy for dynamics
	for (i = 0; i < 2 * UMPC_N * UMPC_NY; ++i) {
		up->u[i] = up->l[i];
	}
	// thrust lims
	for (k = 0; k < UMPC_N; ++k) {
		up->l[2*UMPC_N*UMPC_NY + k] = -up->T0;
		up->u[2*UMPC_N*UMPC_NY + k] = up->Tmax - up->T0;
	}

	// Update matrix ---
	for (i = 0; i < up->nAxdt; ++i) {
		up->Ax_data[i] = i < up->nAxT0dt ? up->dt * up->T0 : up->dt;
	}
	offs = up->nAxdt;
	for (k = 0; k < UMPC_N; ++k) {
		for (i = 0; i < 3; ++i) {
			up->Ax_data[offs + i] = up->dt * s0[i];
		}
		offs += 3;
	}
	for (k = 0; k < UMPC_N; ++k) {
		for (i = 0; i < 6; ++i) {
			up->Ax_data[offs + i] = up->dt * Btau[i];
		}
		offs += 6;
	}
}

static void updateObjective(UprightMPC_t *up, const float ydes[/* 6 */], const float dydes[/* 6 */]) {
	int offs, k, i;
	// q, P diag ---
	offs = 0;
	for (k = 0; k < UMPC_N; ++k) {
		for (i = 0; i < UMPC_NY; ++i) {
			up->Px_data[offs + i] = k == UMPC_N-1 ? up->Qyf[i] : up->Qyr[i];
			up->q[offs + i] = -up->Px_data[offs + i] * ydes[i];
		}
		offs += UMPC_NY;
	}
	for (k = 0; k < UMPC_N; ++k) {
		for (i = 0; i < UMPC_NY; ++i) {
			up->Px_data[offs + i] = k == UMPC_N-1 ? up->Qdyf[i] : up->Qdyr[i];
			up->q[offs + i] = -up->Px_data[offs + i] * dydes[i];
		}
		offs += UMPC_NY;
	}
	for (k = 0; k < UMPC_N; ++k) {
		for (i = 0; i < UMPC_NU; ++i) {
			up->Px_data[offs + i] = up->R[i];
			// Last rows of q remain 0
		}
		offs += UMPC_NU;
	}
}

int umpcUpdate(UprightMPC_t *up, float uquad[/* 3 */], float accdes[/* 6 */], const float p0[/* 3 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */], const float sdes[/* 3 */], float actualT0) {
	static float s0[3], ds0[3], y0[UMPC_NY], dy0[UMPC_NY], ydes[UMPC_NY], dydes[UMPC_NY], dummy[9], Btau[9];
	static float dy1des[UMPC_NY], dq1des[UMPC_NY], e3hR0T[3*3];
	int i, ret;

	// If actual thrust is known, use that
	if (actualT0 >= 0)
		up->T0 = actualT0;

	// Compute some states
	memcpy(s0, &R0[6], 3 * sizeof(float)); // column major R0, and want third col
	// ds0 = -R0 @ e3h @ dq0[3:6] # omegaB
	matMult(dummy, up->e3h, &dq0[3], 3, 1, 3, 1.0f, 0, 0); // e3h*omegaB
	matMult(ds0, R0, dummy, 3, 1, 3, -1.0f, 0, 0); // -R0*e3h*omegaB
	// Btau = (-R0 @ e3h @ self.Ibi)[:,:2]
	matMult(Btau, R0, up->e3hIbi, 3, 3, 3, -1.0f, 0, 0); // -R0*e3h*Ibi; col major => can take first 6 elems
	
	for (i = 0; i < UMPC_NY; ++i) {
		// y0 = np.hstack((p0, s0))
		y0[i] = i < 3 ? p0[i] : s0[i - 3];
		// dy0 = np.hstack((dq0[:3], ds0))
		dy0[i] = i < 3 ? dq0[i] : ds0[i - 3];
		// position parts of ydes
		if (i < 3) {
			ydes[i] = pdes[i];
			dydes[i] = dpdes[i];
		}
	}
	// s parts of ydes
	ydes[3] = sdes[0];
	ydes[4] = sdes[1];
	ydes[5] = sdes[2];
	dydes[3] = dydes[4] = dydes[5] = 0;

	umpcUpdateConstraint(up, s0, Btau, y0, dy0);
	updateObjective(up, ydes, dydes);

	// Update
	osqp_update_bounds(&workspace, up->l, up->u);
	osqp_update_lin_cost(&workspace, up->q);
	osqp_update_P_A(&workspace, up->Px_data, OSQP_NULL, UMPC_NX, up->Ax_data, up->Ax_idx, UMPC_nAdata);
	ret = osqp_solve(&workspace);

	// copy solution
	for (i = 0; i < UMPC_NU; ++i) {
		uquad[i] = workspace.solution->x[2*UMPC_NY*UMPC_N + i];
	}
	up->T0 += uquad[0];
	uquad[0] = up->T0;

	for (i = 0; i < UMPC_NY; ++i) {
		dy1des[i] = workspace.solution->x[UMPC_NY*UMPC_N + i];
		if (i < 3)
			dq1des[i] = dy1des[i];
	}
	matMult(e3hR0T, up->e3h, R0, 3, 3, 3, 1.0f, 0, 1); // e3h*R0^T
	matMult(&dq1des[3], e3hR0T, &dy1des[3], 3, 1, 3, 1.0f, 0, 0); // e3h*R0^T*s1des

	for (i = 0; i < UMPC_NY; ++i) {
		accdes[i] = (dq1des[i] - dq0[i]) / up->dt;
	}

	return ret;
}

// S function ---
UprightMPC_t _up;
int _inited = 0;

void umpcS(float uquad[/* 3 */], float accdes[/* 6 */], const float p0[/* 3 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */], const float sdes[/* 3 */], float dt, float g, float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, const float Ib[/* 3 */], int maxIter, float actualT0) {
	if (_inited == 0) {
		umpcInit(&_up, dt, g, TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom, Ib, maxIter);
		_inited = 1;
	}
	umpcUpdate(&_up, uquad, accdes, p0, R0, dq0, pdes, dpdes, sdes, actualT0);
}
