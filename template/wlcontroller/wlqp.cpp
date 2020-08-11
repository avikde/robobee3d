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

u_t wlqp(int maxIter) {
	osqp_update_max_iter(&workspace, maxIter);
	osqp_update_check_termination(&workspace, 0); // don't check at all

  /* int res = */ osqp_solve(&workspace);
	// memcpy(f, (&workspace)->solution->x, 12 * sizeof(float));

  return mc_u((float *)(&workspace)->solution->x);
}
