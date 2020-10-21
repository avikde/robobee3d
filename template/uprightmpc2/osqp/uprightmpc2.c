/**
 * @file uprightmpc2.c
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-10-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "uprightmpc2.h"
#include "matmult.h"
#include "osqp.h"
#include "workspace.h"
#include <stdio.h>

static void uprightMPCUpdate(UprightMPC_t *up, float T0sp, const float s0s[/* 3*N */], const float Btaus[/*  */], const float y0[/* 6 */], const float dy0[/* 6 */]);


void umpcInit(UprightMPC_t *up, float dt, float g, const float smin[/* 3 */], const float smax[/* 3 */], float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom) {

}

void umpcUpdate(UprightMPC_t *up, float uquad[/* 3 */], float accdes[/* 6 */], const float p0[/* 6 */], const float R0[/* 9 */], const float dq0[/* 6 */], const float pdes[/* 3 */], const float dpdes[/* 3 */]) {
	// TODO:
	// for (int i = 0; i < 9; ++i) {
	// 	printf("%.2f ", R0[i]);
	// }
	// printf("\n");
	matMult(uquad, R0, p0, 3, 1, 3, 1.0f, 1, 0);
}
