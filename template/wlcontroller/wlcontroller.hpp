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
  WLController();
  w_t wrenchMap(const u_t &);
  dw_du_t wrenchJacMap(const u_t &);
protected:
  FunApprox fa[6];
};

void wlcWrapper(double uout[/* 4 */], const double u0[/* 4 */], const double p0[/* 6 */], const double h0[/* 6 */], const double pdes[/* 6 */], const double kpmom[/* 6 */], const double Qdiag[/* 6 */]);
