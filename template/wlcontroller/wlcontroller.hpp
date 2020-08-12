/**
 * @file wlcontroller.hpp
 * @author Avik De (avikde@gmail.com)
 * @brief Use the wlqp with numerical quadratic fit wrench map
 * @version 0.1
 * @date 2020-08-12
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once
#include "wlqp.hpp"
#include "funapprox.hpp"

class WLController : public WrenchLinQP {
public:
  WLController(const float *popts);
  w_t wrenchMap(const u_t &);
  dw_du_t wrenchJacMap(const u_t &);
protected:
  FunApprox fa;
};
