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

WLController::WLController(const float *popts) : fa(popts) {
}

w_t WLController::wrenchMap(const u_t &) {
  return w_t::Zero();
}

dw_du_t WLController::wrenchJacMap(const u_t &) {
  return dw_du_t::Zero();
}
