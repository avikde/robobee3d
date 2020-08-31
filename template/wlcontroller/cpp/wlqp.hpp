/**
 * @file wlqp.hpp
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-11
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once

#include "eigenutil.hpp"
#include <limits>

/**
 * @brief Abstract base class for wrench lin QP. Assumes wrench in R^6, nu = 4;
 */
class WrenchLinQP {
public:
  virtual ~WrenchLinQP() {}
  const float nan = std::numeric_limits<float>::quiet_NaN();
  // Define these to instantiate
  virtual w_t wrenchMap(const u_t &) = 0;
  virtual dw_du_t wrenchJacMap(const u_t &) = 0;

  // Functions to call, variables to read and set
  void setLimits(const u_t &umin, const u_t &umax, const u_t &dumax);
  void setWeight(const w_t &Qdiag) { this->Qdiag = Qdiag; }

  u_t update(const u_t &u0, const w_t &h0, const w_t &pdotdes);

  // Solve settings
  int maxIter = 10;
  float eps = 1e-2f;

protected:
  // u_t u0 = u_t::Zero(); 
  u_t umin = u_t(90.0f, -0.5f, -0.2f, -0.1f), umax = u_t(160.0f, 0.5f, 0.2f, 0.1f);
  // QP bound *not* same as limit on u
  u_t U0 = u_t(0.5f, 1e-3f, 1e-3f, 1e-3f);
  // Keep track of this
  w_t w0 = w_t::Zero();
  w_t Qdiag = w_t(1.0f, 1.0f, 1.0f, 0.1f, 0.1f, 0.1f);
  
  // Call osqp. NOTE: assumes P is symmetric.
  u_t solve(const Eigen::Matrix4f &P, /* const Eigen::Matrix4f &A,  */const u_t &q, const u_t &L, const u_t &U);
};
