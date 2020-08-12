/**
 * @file funapprox.hpp
 * @author Avik De (avikde@gmail.com)
 * @brief Simple quadratic function approximation
 * @version 0.1
 * @date 2020-08-12
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once
#include "eigenutil.hpp"

class FunApprox {
public:
  FunApprox(const float *popts);

  w_t wrenchMap(const u_t &u);
  dw_du_t wrenchJacMap(const u_t &u);
};
