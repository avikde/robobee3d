/**
 * @file uprightmpc2.h
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-10-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once


#ifdef __cplusplus
extern "C" {
#endif

// The matMult function must be defined somewhere - dependent on C++ or C

typedef struct {
	float dt, g, TtoWmax, ws, wds, wpr, wpf, wvr, wthrust, wmom;
  float smin[3], smax[3];
} UprightMPC_t;

void umpcInit(UprightMPC_t *up, float dt, float g, const float smin[/* 3 */], const float smax[/* 3 */], float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom);

void umpcUpdate(UprightMPC_t *up, float uquad[/* 3 */], float accdes[/* 6 */], const float p0[/* 6 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */]);

#ifdef __cplusplus
}
#endif

