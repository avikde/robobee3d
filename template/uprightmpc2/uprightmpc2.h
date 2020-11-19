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
#define UMPC_NC (2*UMPC_N*UMPC_NY + UMPC_N)
#define UMPC_nAdata 48 // depends on N, printed out in python script

typedef struct {
	float dt, g, Tmax;
  // Weights
  float Qyr[6], Qyf[6], Qdyr[6], Qdyf[6], R[3];
  // Limits
  float smin[3], smax[3];
  // Other params
  float e3h[3*3]; // constant matrix skew(e3)
  float e3hIbi[3*3]; // e3h*Ibi
  // Workspace
  float l[UMPC_NC], u[UMPC_NC], q[UMPC_NX];
  float Px_data[UMPC_NX]; // all diagonal elements (checked size of Pdata_i=NX)
  float Ax_data[UMPC_nAdata];
  int Ax_idx[UMPC_nAdata], nAxT0dt, nAxdt;
  float c0[UMPC_NY];
  float T0;
} UprightMPC_t;

void umpcInit(UprightMPC_t *up, float dt, float g, float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, const float Ib[/* 3 */], int maxIter);

int umpcUpdate(UprightMPC_t *up, float uquad[/* 3 */], float accdes[/* 6 */], const float p0[/* 3 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */], const float sdes[/* 3 */], float actualT0);

void umpcS(float uquad_y1[/* 3 */], float accdes_y2[/* 6 */], const float p0_u1[/* 3 */], const float R0_u2[/* 9 */], const float dq0_u3[/* 6 */], const float pdes_u4[/* 3 */], const float dpdes_u5[/* 3 */], const float sdes_u6[/* 3 */], float dt_u7, float g_u8, float TtoWmax_u9, float ws_u10, float wds_u11, float wpr_u12, float wpf_u13, float wvr_u14, float wvf_u15, float wthrust_u16, float wmom_u17, const float Ib_u18[/* 3 */], int maxIter_u19, float actualT0_u20);

#ifdef __cplusplus
}
#endif

