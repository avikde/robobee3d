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
#include "funapprox.h"

#ifdef __cplusplus
extern "C" {
#endif

// The matMult function must be defined somewhere - dependent on C++ or C

// Must be changed along with the autogen code
#define UMPC_N 3
#define UMPC_NY 6
#define UMPC_NU 3
#define WLQP_NU 4
#define WLQP_NW 6
#define UMPC_NX (UMPC_N*(2*UMPC_NY + UMPC_NU) + WLQP_NU)
#define UMPC_NC (2*UMPC_N*UMPC_NY + UMPC_N + WLQP_NU)
#define UMPC_nAdata 48 // depends on N, printed out in python script
#define UMPC_nPdata 94 // depends on N, printed out in python script

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
  float Px_data[UMPC_nPdata];
  float Ax_data[UMPC_nAdata];
  int Ax_idx[UMPC_nAdata], nAxT0dt, nAxdt;
  float c0[UMPC_NY];
  float T0;
  // WLQP stuff
  float u0[WLQP_NU], umin[WLQP_NU], umax[WLQP_NU], dumax[WLQP_NU];
  float M0[6*6], Qw[6*6], Tform0[6*6], h0W[3];
  FunApprox_t fa[6];
} UprightMPC_t;

void umpcInit(UprightMPC_t *up, float dt, float g, float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, float mb, const float Ib[/* 3 */], const float umin[/* 4 */], const float umax[/* 4 */], const float dumax[/* 4 */], const float Qw[/* 6 */], float controlRate, int maxIter, const float popts[/* 90 */]);

int umpcUpdate(UprightMPC_t *up, float uquad[/* 3 */], float accdes[/* 6 */], float uwlqp[/* 4 */], const float p0[/* 6 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */]);

void umpcS(float uquad_y1[/* 3 */], float accdes_y2[/* 6 */], float uwlqp_y3[/* 4 */],const float p0_u1[/* 3 */], const float R0_u2[/* 9 */], const float dq0_u3[/* 6 */], const float pdes_u4[/* 3 */], const float dpdes_u5[/* 3 */], float dt_u6, float g_u7, float TtoWmax_u8, float ws_u9, float wds_u10, float wpr_u11, float wpf_u12, float wvr_u13, float wvf_u14, float wthrust_u15, float wmom_u16, float mb_u17, const float Ib_u18[/* 3 */], const float umin_u19[/* 4 */], const float umax_u20[/* 4 */], const float dumax_u21[/* 4 */], const float Qw_u22[/* 6 */], float controlRate_u23, int maxIter_u24, const float *popts_u25);

#ifdef __cplusplus
}
#endif

