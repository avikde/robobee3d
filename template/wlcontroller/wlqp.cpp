/**
 * @file wlqp.cpp
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-11
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "wlqp.hpp"
#include <osqp.h>
#include <cmath>
extern "C" {
#include <workspace.h>
}

u_t WrenchLinQP::update(const u_t &u0, const w_t &p0, const w_t &h0, const w_t &pdes, const w_t &Qdiag, const w_t &kpmom) {
	if (!wrenchMap || !wrenchJacMap)
		return u_t::Zero();
	
	// Momentum reference dynamics https://github.com/avikde/robobee3d/pull/166 TODO: incorporate as MPC
	w_t pdotdes = kpmom.cwiseProduct(pdes - p0);
	return update2(u0, h0, pdotdes, Qdiag);
}

u_t WrenchLinQP::update2(const u_t &u0, const w_t &h0, const w_t &pdotdes, const w_t &Qdiag) {
	// Input rate limit TODO: sparsity
	Eigen::Matrix4f A = Eigen::Matrix4f::Identity();

	auto w0 = wrenchMap(u0);
	auto a0 = w0 - h0 - pdotdes;
	auto A1 = wrenchJacMap(u0);

	auto P = A1.transpose() * Eigen::DiagonalMatrix<float, 6>(Qdiag) * A1;
	u_t q = A1.transpose() * Qdiag.cwiseProduct(a0);
	u_t L = -this->U0;
	u_t U = this->U0;

	// Input limits (not just rate limits)
	if (!std::isnan(umin[0])) {
		for (int i = 0; i < umin.size(); ++i) {
			if (u0[i] < umin[i])
				L[i] = 0; // do not reduce further
			else if (u0[i] > umax[i])
				U[i] = 0; // do not increase further
		}
	}

	auto du = solve(P, A, q, L, U);
	return u0 + du;
}

void WrenchLinQP::setLimits(const u_t &umin, const u_t &umax, const u_t &dumax) {
	this->umin = umin;
	this->umax = umax;
	this->U0 = dumax;
}

u_t WrenchLinQP::solve(const Eigen::Matrix4f &P, const Eigen::Matrix4f &A, const u_t &q, const u_t &L, const u_t &U) {
	OSQPWorkspace *work = &workspace;
	static Eigen::VectorXf Px_data;

	// Osqp init
	osqp_update_max_iter(work, maxIter);
	osqp_update_eps_rel(work, eps);
	osqp_update_eps_abs(work, eps);
	osqp_update_check_termination(work, 0); // don't check at all

	// Update
	eigenUpperTriangularVals<4>(P, Px_data);
	osqp_update_P(work, Px_data.data(), OSQP_NULL, Px_data.size());
	osqp_update_A(work, A.data(), OSQP_NULL, work->data->m * work->data->n);
	osqp_update_lin_cost(work, q.data());
	osqp_update_bounds(work, L.data(), U.data());

	/* int res = */ osqp_solve(work);

	return mc_u((float *)work->solution->x);
}
