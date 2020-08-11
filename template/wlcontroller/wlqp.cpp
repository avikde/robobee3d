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

u_t wlqp(const u_t &u, const dw_du_t &dw_du, int maxIter) {
  OSQPWorkspace *work = &workspace;
	osqp_update_max_iter(work, maxIter);
	osqp_update_check_termination(work, 0); // don't check at all

  // Input rate limit TODO: sparsity
  Eigen::Matrix4f A = Eigen::Matrix4f::Identity();


  // Update
  osqp_update_A(work, (const float *)A.data(), OSQP_NULL, work->data->m * work->data->n);

  /* int res = */ osqp_solve(work);

  return mc_u((float *)work->solution->x);
}
