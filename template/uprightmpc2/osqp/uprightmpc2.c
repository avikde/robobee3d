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

static void PRINTVEC(const float *y, int sz) {
	for (int i = 0; i < sz; ++i) {
		printf("%.2f,", y[i]);
	}
	printf("\n");
}

void umpcInit(UprightMPC_t *up, float dt, float g, const float smin[/* 3 */], const float smax[/* 3 */], float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, const float Ib[/* 3 */], int maxIter) {
	up->dt = dt;
	up->g = g;
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

	// Other updates
	for (int i = 0; i < UMPC_NY; ++i) {
		up->c0[i] = i == 2 ? -up->g : 0;
	}
	up->T0 = 0;
	for (int i = 0; i < 3; ++i) {
		up->Ibi[i] = Ib[i];
	}
	// skew(e3) entered in col major order
	memset(up->e3h, 0, 9 * sizeof(float));
	up->e3h[1] = 1;
	up->e3h[3] = -1;

	for (int i = 0; i < UMPC_NX; ++i) {
		up->q[i] = 0;
	}
	for (int i = 0; i < UMPC_NC; ++i) {
		up->l[i] = up->u[i] = 0;
	}

	// OSQP
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
	
	// Update vector ---
	for (int i = 0; i < UMPC_N * UMPC_NY; ++i) {
		up->l[i] = 0;
	}
	// y1 = y0 + dt * dy0
	for (int i = 0; i < UMPC_NY; ++i) {
		y1[i] = y0[i] + up->dt * dy0[i];
		// l[:ny] = -y1
		up->l[i] = -y1[i];
	}
	// Dynamics constraints
	for (int k = 0; k < UMPC_N; ++k) {
		float *pk = &up->l[UMPC_NY * (UMPC_N + k)];
		for (int i = 0; i < UMPC_NY; ++i) {
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
	for (int i = 0; i < 2 * UMPC_N * UMPC_NY; ++i) {
		up->u[i] = up->l[i];
	}
	// s lims
	for (int k = 0; k < UMPC_N; ++k) {
		float *pl = &up->l[2*UMPC_N*UMPC_NY + 3*k];
		float *pu = &up->u[2*UMPC_N*UMPC_NY + 3*k];
		for (int i = 0; i < 3; ++i) {
			pl[i] = up->smin[i];
			pu[i] = up->smax[i];
		}
	}
	// thrust lims
	for (int k = 0; k < UMPC_N; ++k) {
		up->l[2*UMPC_N*UMPC_NY + 3*UMPC_N + k] = -up->T0;
		up->u[2*UMPC_N*UMPC_NY + 3*UMPC_N + k] = up->Tmax - up->T0;
	}

	// Update matrix ---TODO:
	
}

static void updateObjective(UprightMPC_t *up, const float ydes[/* 6 */], const float dydes[/* 6 */]) {
	// q ---
	int offs = 0;
	for (int k = 0; k < UMPC_N; ++k) {
		for (int i = 0; i < UMPC_NY; ++i) {
			up->q[offs + i] = -(k == UMPC_N-1 ? up->Qyf[i] : up->Qyr[i]) * ydes[i];
		}
		offs += UMPC_NY;
	}
	for (int k = 0; k < UMPC_N; ++k) {
		for (int i = 0; i < UMPC_NY; ++i) {
			up->q[offs + i] = -(k == UMPC_N-1 ? up->Qdyf[i] : up->Qdyr[i]) * dydes[i];
		}
		offs += UMPC_NY;
	}
	// Last rows remain 0
}

int umpcUpdate(UprightMPC_t *up, float uquad[/* 3 */], float accdes[/* 6 */], const float p0[/* 3 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */]) {
	static float s0[3], ds0[3], y0[UMPC_NY], dy0[UMPC_NY], ydes[UMPC_NY], dydes[UMPC_NY], dummy3[3];
	static float dy1des[UMPC_NY], dq1des[UMPC_NY];

	// Compute some states
	memcpy(s0, &R0[6], 3 * sizeof(float)); // column major R0, and want third col
	// ds0 = -R0 @ e3h @ dq0[3:6] # omegaB
	matMult(dummy3, up->e3h, &dq0[3], 3, 1, 3, 1.0f, 0, 0); // d3h*omegaB
	matMult(ds0, R0, dummy3, 3, 1, 3, -1.0f, 0, 0); // -R0*d3h*omegaB
	
	for (int i = 0; i < UMPC_NY; ++i) {
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
	ydes[3] = ydes[4] = 0;
	ydes[5] = 1;
	dydes[3] = dydes[4] = dydes[5] = 0;

	umpcUpdateConstraint(up, NULL, NULL, y0, dy0);
	updateObjective(up, ydes, dydes);

	// Update
	osqp_update_bounds(&workspace, up->l, up->u);
	osqp_update_lin_cost(&workspace, up->q);
	int ret = osqp_solve(&workspace);

	// copy solution
	for (int i = 0; i < UMPC_NU; ++i) {
		uquad[i] = workspace.solution->x[2*UMPC_NY*UMPC_N + i];
	}
	up->T0 += uquad[0];
	uquad[0] = up->T0;

	for (int i = 0; i < UMPC_NY; ++i) {
		dy1des[i] = workspace.solution->x[UMPC_NY*UMPC_N + i];
		if (i < 3)
			dq1des[i] = dy1des[i];
	}
	matMult(dummy3, R0, &dy1des[3], 3, 1, 3, 1.0f, 1, 0); // R0^T*s1des
	matMult(&dq1des[3], up->e3h, dummy3, 3, 1, 3, 1.0f, 0, 0); // e3h*R0^T*s1des

	for (int i = 0; i < UMPC_NY; ++i) {
		accdes[i] = (dq1des[i] - dq0[i]) / up->dt;
	}

	return ret;
}
