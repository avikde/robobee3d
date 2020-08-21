/**
 * @file wlmex.h
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void wlcWrapper(double uout[/* 4 */], const double u0[/* 4 */], const double p0[/* 6 */], const double h0[/* 6 */], const double pdes[/* 6 */], const double kpmom[/* 6 */], const double Qdiag[/* 6 */]);

#ifdef __cplusplus
}
#endif
