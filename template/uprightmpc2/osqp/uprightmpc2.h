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

// Must be changed along with the autogen code
#define UMPC_N 3
#define UMPC_NY 6
#define UMPC_NU 3
#define UMPC_NX (UMPC_N*(2*UMPC_NY + UMPC_NU))
#define UMPC_NC (2*UMPC_N*UMPC_NY + 4*UMPC_N)

typedef struct {
	float dt, g, Tmax;
  // Weights
  float Qyr[6], Qyf[6], Qdyr[6], Qdyf[6], R[3];
  // Limits
  float smin[3], smax[3];
  // Other params
  float Ibi[3]; // inertia inv
  float e3h[3*3]; // constant matrix skew(e3)
  // Workspace
  float l[UMPC_NC], u[UMPC_NC], q[UMPC_NX];
  float c0[UMPC_NY];
  float T0;
} UprightMPC_t;

void umpcInit(UprightMPC_t *up, float dt, float g, const float smin[/* 3 */], const float smax[/* 3 */], float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, const float Ib[/* 3 */], int maxIter);

void umpcUpdate(UprightMPC_t *up, float uquad[/* 3 */], float accdes[/* 6 */], const float p0[/* 6 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */]);

#ifdef __cplusplus
}
#endif

