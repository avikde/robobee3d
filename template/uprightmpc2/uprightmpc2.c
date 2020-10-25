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

static void PRINTVEC(const float *y, int sz) {
	int i;
	for (i = 0; i < sz; ++i) {
		printf("%.2f,", y[i]);
	}
	printf("\n");
}
static void PRINTMAT(const float *M, int sz1, int sz2) {
	int i, j;
	for (i = 0; i < sz1; ++i) {
		for (j = 0; j < sz2; ++j) {
			printf("%.2f,", M[Cind(sz1, i, j)]);
		}
		printf("\n");
	}
}

// writes n*(n+1)/2 values in out
static int getUpperTriang(float *out, const float *Mcolmaj, int n) {
	int i, j, k;
	k = 0;
	for (j = 0; j < n; ++j) {
		for (i = 0; i <= j; ++i) {
			out[k] = Mcolmaj[Cind(n, i, j)];
			k++;
		}
	}
	return k;
}

// Wrench map
static void wrenchMap(const UprightMPC_t *up, float *w, const float *u) {
	int i;
	for (i = 0; i < WLQP_NW; ++i) {
		w[i] = funApproxF(&up->fa[i], u);
	}
}

static void wrenchJacMap(const UprightMPC_t *up, float *dw_du, const float *u) {
	int i, j;
	static float dwi_du[WLQP_NU];
	for (i = 0; i < WLQP_NW; ++i) {
		funApproxDf(dwi_du, &up->fa[i], u);
		for (j = 0; j < WLQP_NU; ++j) {
			// Copy into the col-major matrix
			dw_du[Cind(WLQP_NW, i, j)] = dwi_du[j]; // i >= 2 && i - 2 == j ? 1 : 0;  //
		}
	}
}

void umpcInit(UprightMPC_t *up, float dt, float g, float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, float mb, const float Ib[/* 3 */], const float umin[/* 4 */], const float umax[/* 4 */], const float dumax[/* 4 */], const float Qw[/* 6 */], float controlRate, int maxIter, const float popts[/* 90 */]) {
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
	// Limits
	for (i = 0; i < WLQP_NU; ++i) {
		up->u0[i] = 0;
		up->umin[i] = umin[i];
		up->umax[i] = umax[i];
		up->dumax[i] = dumax[i] / controlRate;
	}

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
	// Left 1/4, T0*dt
	offs = 0;
	n2 = 2*UMPC_NY + 3; // nnz in each block col on the left
	for (k = 0; k < UMPC_N-2; ++k) {
		up->Ax_idx[offs+0] = n2*k + 8;
		up->Ax_idx[offs+1] = n2*k + 11;
		up->Ax_idx[offs+2] = n2*k + 14;
		offs += 3;
	}
	up->nAxT0dt = offs;

	// Middle 2/4, dt
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

	// Right 3/4
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

	// WLQP ---
	for (i = 0; i < 6; ++i) {
		funApproxInit(&up->fa[i], &popts[15 * i]);
	}
	memset(up->M0, 0, 36 * sizeof(float));
	memset(up->Qw, 0, 36 * sizeof(float));
	memset(up->Tform0, 0, 36 * sizeof(float));
	for (i = 0; i < 6; ++i) {
		up->M0[Cind(6, i, i)] = i < 3 ? mb : Ib[i - 3];
		up->Qw[Cind(6, i, i)] = Qw[i];
		up->Tform0[Cind(6, i, i)] = 1; // Identity
	}
	// Gravity vector world frame
	up->h0W[0] = up->h0W[1] = 0;
	up->h0W[2] = up->M0[0] * up->g;
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
	// Delta-u input (rate) limits
	float *delUL = &up->l[2*UMPC_N*UMPC_NY + UMPC_N];
	float *delUU = &up->u[2*UMPC_N*UMPC_NY + UMPC_N];
	for (i = 0; i < WLQP_NU; ++i) {
		delUL[i] = -up->dumax[i];
		delUU[i] = up->dumax[i];
		// Cannot decrease if at limit and vice versa
		if (up->u0[i] < up->umin[i])
			delUL[i] = 0;
		else if (up->u0[i] > up->umax[i])
			delUU[i] = 0;
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

static void updateObjective(UprightMPC_t *up, const float ydes[/* 6 */], const float dydes[/* 6 */], const float dwdu0[/* 6*4 */], const float w0t[/* 6 */], const float M0t[/* 6*6 */]) {
	int offsq, offsP, k, i, ii, jj;
	static float dummy66[6*6], dummy44[4*4], Qwdwdu[6*4], dy1delu[6*4], deludelu[4*4], M0TQwM0[6*6], Qww0t[6], mM0tQww0t[6], dwduQww0t[4];

	// // Last block col
	// matMult(Qwdwdu, up->Qw, dwdu0, 6, 4, 6, 1.0f, 0, 0); // Qw*dwdu
	// matMult(dy1delu, M0t, Qwdwdu, 6, 4, 6, -1.0f, 1, 0); // -M0^T*Qw*dwdu
	// matMult(deludelu, dwdu0, Qwdwdu, 4, 4, 6, 1.0f, 1, 0); // dwdu^T*Qw*dwdu
	// matMult(dummy66, up->Qw, M0t, 6, 6, 6, 1.0f, 0, 0); // Qw*M0
	// matMult(M0TQwM0, M0t, dummy66, 6, 6, 6, 1.0f, 1, 0); // M0^T*Qw*M0
	// matMult(Qww0t, up->Qw, w0t, 6, 1, 6, 1.0f, 0, 0); // - M0t.T @ Qw @ w0t
	// matMult(mM0tQww0t, M0t, Qww0t, 6, 1, 6, -1.0f, 1, 0); // - M0t.T @ Qw @ w0t
	// matMult(dwduQww0t, dwdu0, Qww0t, 4, 1, 6, 1.0f, 1, 0); // dwdu.T @ Qw @ w0t

	// q, P diag ---
	offsq = offsP = 0;

	// First N*ny
	for (k = 0; k < UMPC_N; ++k) {
		for (i = 0; i < UMPC_NY; ++i) {
			up->Px_data[offsP + i] = k == UMPC_N-1 ? up->Qyf[i] : up->Qyr[i];
			up->q[offsq + i] = -up->Px_data[offsP + i] * ydes[i];
		}
		offsP += UMPC_NY;
		offsq += UMPC_NY;
	}
	
	// Second N*ny
	for (k = 0; k < UMPC_N; ++k) {
		if (k == 0) {
			// dy1,dy1 block upper triang
			// create diag matrix
			for (ii = 0; ii < 6; ++ii) {
				for (jj = 0; jj < 6; ++jj) {
					dummy66[Cind(6, ii, jj)] = ((ii == jj) ? up->Qdyr[ii] : 0);// + M0TQwM0[Cind(6, ii, jj)]; // add M0T0dt.T @ np.diag(Qw) @ M0T0dt
				}
			}
			offsP += getUpperTriang(&up->Px_data[offsP], dummy66, UMPC_NY);
		} else {
			for (i = 0; i < UMPC_NY; ++i) {
				up->Px_data[offsP + i] = k == UMPC_N-1 ? up->Qdyf[i] : up->Qdyr[i];
			}
			offsP += UMPC_NY;
		}
		for (i = 0; i < UMPC_NY; ++i) {
			if (k == 0) {
				up->q[offsq + i] = -up->Qdyr[i] * dydes[i];// + mM0tQww0t[i];
			} else {
				up->q[offsq + i] = -(k == UMPC_N-1 ? up->Qdyf[i] : up->Qdyr[i]) * dydes[i];
			}
		}
		offsq += UMPC_NY;
	}

	// Third N*nu
	for (k = 0; k < UMPC_N; ++k) {
		for (i = 0; i < UMPC_NU; ++i) {
			up->Px_data[offsP + i] = up->R[i];
			// N*nu of q remain 0
		}
		offsP += UMPC_NU;
		offsq += UMPC_NU;
	}

	// Last wlqp part for q
	for (i = 0; i < 4; ++i) {
		up->q[offsq + i] = 0;//dwduQww0t[i];
	}

	// populate lastcol
	for (jj = 0; jj < 4; ++jj) {
		for (ii = 0; ii < 6; ++ii) {
			up->Px_data[offsP + ii] = 0;//dy1delu[Cind(6, ii, jj)];
		}
		offsP += 6;
		for (ii = 0; ii < jj+1; ++ii) {
			up->Px_data[offsP + ii] = 0;//deludelu[Cind(6, ii, jj)];
		}
		offsP += (jj+1);
	}
}

int umpcUpdate(UprightMPC_t *up, float uquad[/* 3 */], float accdes[/* 6 */], float uwlqp[/* 4 */], const float p0[/* 3 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */]) {
	static float s0[3], ds0[3], y0[UMPC_NY], dy0[UMPC_NY], ydes[UMPC_NY], dydes[UMPC_NY], dummy[9], Btau[9];
	static float dy1des[UMPC_NY], dq1des[UMPC_NY], dwdu0[6*6], w0t[6], M0t[6*6], T0[6*6], e3hR0T[3*3];
	int i, j, ret;

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
	ydes[3] = ydes[4] = 0;
	ydes[5] = 1;
	dydes[3] = dydes[4] = dydes[5] = 0;

	// WLQP stuff. sample wrench map ---
	// wrenchMap(up, w0t, up->u0);
	// wrenchJacMap(up, dwdu0, up->u0);
	// Create Tform0
	matMult(e3hR0T, up->e3h, R0, 3, 3, 3, 1.0f, 0, 1); // e3h*R0^T
	// for (i = 0; i < 3; ++i) {
	// 	for (j = 0; j < 3; ++j) {
	// 		up->Tform0[Cind(6, 3+i, 3+j)] = e3hR0T[Cind(3, i, j)]; // T0[3:,3:] = e3h @ R0.T
	// 	}
	// }
	// // M0t = M0*T0/dt
	// matMult(M0t, up->M0, up->Tform0, 6, 6, 6, 1 / up->dt, 0, 0);
	// // w0t = w0 - h0 + M0*dq0/dt
	// matMult(dummy, up->M0, dq0, 6, 1, 6, 1 / up->dt, 0, 0); // M0*dq0/dt
	// matMult(&dummy[6], R0, up->h0W, 3, 1, 3, 1.0f, 1, 0); // h0B = R0^T*h0W
	// for (i = 0; i < 6; ++i) {
	// 	w0t[i] = w0t[i] - (i < 3 ? dummy[6+i] : 0) + dummy[i];
	// }

	umpcUpdateConstraint(up, s0, Btau, y0, dy0);
	updateObjective(up, ydes, dydes, dwdu0, w0t, M0t);

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
	matMult(&dq1des[3], e3hR0T, &dy1des[3], 3, 1, 3, 1.0f, 0, 0); // e3h*R0^T*s1des

	for (i = 0; i < UMPC_NY; ++i) {
		accdes[i] = (dq1des[i] - dq0[i]) / up->dt;
	}

	// WLQP
	for (i = 0; i < 4; ++i) {
		up->u0[i] += workspace.solution->x[(2*UMPC_NY + UMPC_NU)*UMPC_N + i]; // delta u
		uwlqp[i] = up->u0[i];
	}

	return ret;
}

// S function ---
UprightMPC_t _up;
int _inited = 0;

void umpcS(float uquad_y1[/* 3 */], float accdes_y2[/* 6 */], float uwlqp_y3[/* 4 */], const float p0_u1[/* 3 */], const float R0_u2[/* 9 */], const float dq0_u3[/* 6 */], const float pdes_u4[/* 3 */], const float dpdes_u5[/* 3 */], float dt_u6, float g_u7, float TtoWmax_u8, float ws_u9, float wds_u10, float wpr_u11, float wpf_u12, float wvr_u13, float wvf_u14, float wthrust_u15, float wmom_u16, float mb_u17, const float Ib_u18[/* 3 */], const float umin_u19[/* 4 */], const float umax_u20[/* 4 */], const float dumax_u21[/* 4 */], const float Qw_u22[/* 6 */], float controlRate_u23, int maxIter_u24, const float *popts_u25) {
	if (_inited == 0) {
		umpcInit(&_up, dt_u6, g_u7, TtoWmax_u8, ws_u9, wds_u10, wpr_u11, wpf_u12, wvr_u13, wvf_u14, wthrust_u15, wmom_u16, mb_u17, Ib_u18, umin_u19, umax_u20, dumax_u21, Qw_u22, controlRate_u23, maxIter_u24, popts_u25);
		_inited = 1;
	}
	umpcUpdate(&_up, uquad_y1, accdes_y2, uwlqp_y3, p0_u1, R0_u2, dq0_u3, pdes_u4, dpdes_u5);
}
