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

void wlqpInit(WLQP_t *wlqp) {
	// populate default values
	wlqp->Qdiag[0] = wlqp->Qdiag[1] = wlqp->Qdiag[2] = 1.0f;
	wlqp->Qdiag[3] = wlqp->Qdiag[4] = wlqp->Qdiag[5] = 0.1f;
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
	wlqp->U0[1] = 0.01f;
	wlqp->U0[2] = 0.01f;
	wlqp->U0[3] = 0.01f;
}

void wlqpUpdate(WLQP_t *wlqp, const float *u0, const float *h0, const float *pdotdes) {
	// // Input rate limit TODO: sparsity
	// Eigen::Matrix4f A = Eigen::Matrix4f::Identity();

	// auto w0 = wrenchMap(u0);
	// auto a0 = w0 - h0 - pdotdes;
	// auto A1 = wrenchJacMap(u0);

	// auto P = A1.transpose() * Qdiag.asDiagonal() * A1;
	// u_t q = A1.transpose() * Qdiag.cwiseProduct(a0);
	// u_t L = -this->U0;
	// u_t U = this->U0;

	// // Input limits (not just rate limits)
	// if (!std::isnan(umin[0])) {
	// 	for (int i = 0; i < umin.size(); ++i) {
	// 		if (u0[i] < umin[i])
	// 			L[i] = 0; // do not reduce further
	// 		else if (u0[i] > umax[i])
	// 			U[i] = 0; // do not increase further
	// 	}
	// }

	// // std::cout << w0 << A1;
	// // std::cout << P << q << L << U << A;

	// auto du = solve(P, A, q, L, U);
	// return u0 + du;
}

