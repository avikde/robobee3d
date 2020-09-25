/**
 * @file wlcontroller.cpp
 * @author Avik De (avikde@gmail.com)
 * @brief Use the wlqp with numerical quadratic fit wrench map
 * @version 0.1
 * @date 2020-08-12
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "wlcontroller.hpp"
// #include <iostream>

WLController::WLController(const Eigen::VectorXf &popts, float controlRate) : WrenchLinQP(controlRate) {
  // load popts
  auto poptsm = popts.reshaped<Eigen::RowMajor>(6, 15);
  for (int i = 0; i < 6; ++i)
    fa[i].init(poptsm.row(i));
}

w_t WLController::wrenchMap(const u_t &u) {
  static w_t w;
  for (int i = 0; i < w.size(); ++i) {
    w[i] = fa[i].f(u);
  }
  return w;
}

dw_du_t WLController::wrenchJacMap(const u_t &u) {
  static dw_du_t dw_du;
  for (int i = 0; i < dw_du.rows(); ++i) {
    dw_du.row(i) = fa[i].Df(u);
  }
  return dw_du;
}
