/**
 * @file ematmult.cpp
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-11-18
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "ematmult.hpp"

// Wrapper for matrix multiply - use Eigen for it
extern "C" void matMult(float *C, const float *A, const float *B, const int m, const int n, const int k, const float alpha, int AT, int BT) {
	auto opC = MMatX(C, m, n);
	if (AT == 0 && BT == 0) {
		opC = alpha * MCMatX(A, m, k) * MCMatX(B, k, n);
	} else if (AT == 0 && BT == 1) {
		opC = alpha * MCMatX(A, m, k) * MCMatX(B, n, k).transpose();
	} else if (AT == 1 && BT == 0) {
		opC = alpha * MCMatX(A, k, m).transpose() * MCMatX(B, k, n);
	} else {
		opC = alpha * MCMatX(A, k, m).transpose() * MCMatX(B, n, k).transpose();
	}
}
