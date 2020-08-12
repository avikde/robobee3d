/**
 * @file funapprox.hpp
 * @author Avik De
 * @brief Scalar function approx
 * @version 0.1
 * @date 2020-08-12
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once
#include "eigenutil.hpp"

/* template <int k> */
class FunApprox {
public:
	typedef Eigen::Matrix<float, /* 1 + k + k *(k + 1) / 2 */15, 1> popts_t;

	void init(const popts_t &p) {
		// unpack
		a0 = p[0];
		a1 = p.segment</*  k  */4>(1);
		// eigenUpperTriangularSet<4>(A2, p.segment</*  k *(k + 1) / 2  */10>(/* k + 1 */5));
		// eigenUpperTriangularSet<4>(A2, p.segment<10>(5));
			
		int k = 0, N = 4;
		auto tvals = p.segment<10>(5);
		for (int j = 0; j < N; ++j) {
			for (int i = 0; i <= j; ++i) {
				A2(i, j) = tvals[k];
				k++;
			}
		}
		// Use A2.selfadjointView<Eigen::Upper>() after this
	}

protected:
	float a0;
	Eigen::Matrix<float, /* k */4, 1> a1;
	Eigen::Matrix<float, /* k, k */4, 4> A2;
};
