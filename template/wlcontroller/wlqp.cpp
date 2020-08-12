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
extern "C" {
#include <workspace.h>
}

u_t wlqp(const u_t &u, const dw_du_t &dw_du, const w_t &w0, int maxIter) {
	OSQPWorkspace *work = &workspace;
	static Eigen::Matrix<float, 4 * (4+1)/2, 1> Px_data;

	// Osqp init
	osqp_update_max_iter(work, maxIter);
	osqp_update_check_termination(work, 0); // don't check at all

	// Input rate limit TODO: sparsity
	Eigen::Matrix4f A = Eigen::Matrix4f::Identity();

	// auto a0 = w0 - h0 - pdotdes;
	const auto &A1 = dw_du;

	w_t Qdiag = q_t::Zero();

	auto P = A1.transpose() * Eigen::DiagonalMatrix<float, 6>(Qdiag) * A1;
	auto q = A1.transpose() * Qdiag.cwiseProduct(a0);
	auto l = -u;

	// q = A1.T @ np.diag(Qd) @ a0
	// L = -self.U
	// U = np.copy(self.U)

	// # Input limits (not just rate limits)
	// if self.umin is not None:
	//     for i in range(len(self.umin)):
	//         if self.u0[i] < self.umin[i]:
	//             L[i] = 0 # do not reduce further
	//         elif self.u0[i] > self.umax[i]:
	//             U[i] = 0 # do not increase further

	// # update OSQP
	// Px = P[np.tril_indices(P.shape[0])] # need in col order
	// self.model.update(Px=Px, q=q, l=L, u=U, Ax=np.ravel(A))
	// res = self.model.solve()
	// self.u0 = res.x + self.u0
	// # self.u0[0] = 110 + 3 * (pdes[2] - p0[2])

	// if self.n >= 4:
	//     self.w0 = self.wrenchMap(self.u0)
	//     return self.u0

	// Update
	eigenUpperTriangularVals<4>(P, Px_data);
	Px_idx.resize(Pu0_inds.size() + P_yaw_inds.size());
	for (Eigen::Index i = 0; i < Px_idx.size(); ++i) {
		Px_idx[i] = (c_int)(i < Pnnz_u0 ? Pu0_inds[i] : P_yaw_inds[i - Pnnz_u0]);
	}
	osqp_update_P(work, Px_data);
	osqp_update_A(work, A.data(), OSQP_NULL, work->data->m * work->data->n);

	/* int res = */ osqp_solve(work);

	return mc_u((float *)work->solution->x);
}
