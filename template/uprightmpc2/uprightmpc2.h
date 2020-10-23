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
#define WLQP_NU 4
#define UMPC_NX (UMPC_N*(2*UMPC_NY + UMPC_NU) + WLQP_NU)
#define UMPC_NC (2*UMPC_N*UMPC_NY + UMPC_N + WLQP_NU)
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
  // WLQP stuff
  float u0[WLQP_NU], umin[WLQP_NU], umax[WLQP_NU], dumax[WLQP_NU];
  float M0[6*6], Qw[6*6];
} UprightMPC_t;

void umpcInit(UprightMPC_t *up, float dt, float g, float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, float mb, const float Ib[/* 3 */], const float umin[/* 4 */], const float umax[/* 4 */], const float dumax[/* 4 */], const float Qw[/* 6 */], float controlRate, int maxIter);

int umpcUpdate(UprightMPC_t *up, float uquad[/* 3 */], float accdes[/* 6 */], const float p0[/* 6 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */]);

/**
 * @brief For simulink RT, a single function with a bunch of arguments
 * 
 * @param uquad y1
 * @param accdes y2
 * @param p0 u1
 * @param R0 u2
 * @param dq0 u3
 * @param pdes u4
 * @param dpdes u5
 * @param dt u6
 * @param g u7
 * @param smin u8 
 * @param smax u9
 * @param TtoWmax u10
 * @param ws u11
 * @param wds u12
 * @param wpr u13
 * @param wpf u14
 * @param wvr u15
 * @param wvf u16
 * @param wthrust u17 
 * @param wmom u18
 * @param Ib u19
 * @param maxIter u20
 */
void umpcS(float uquad[/* 3 */], float accdes[/* 6 */], const float p0[/* 6 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */], float dt, float g, const float smin[/* 3 */], const float smax[/* 3 */], float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, const float Ib[/* 3 */], int maxIter);

#ifdef __cplusplus
}
#endif

