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

u_t wlqp(const u_t &u, const dw_du_t &dw_du, int maxIter);
