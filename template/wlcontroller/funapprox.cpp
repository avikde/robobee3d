/**
 * @file funapprox.cpp
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-12
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "funapprox.hpp"

FunApprox::FunApprox(const float *) {
  // TODO:
}

w_t FunApprox::wrenchMap(const u_t &u) {
  return w_t::Zero();
}

dw_du_t FunApprox::wrenchJacMap(const u_t &u) {
  return dw_du_t::Zero();
}
