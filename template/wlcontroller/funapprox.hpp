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
	const static int k = 4; // TODO: template

	typedef Eigen::Matrix<float, 1 + k + k * (k + 1) / 2, 1> popts_t;
	typedef Eigen::Matrix<float, k, 1> x_t;
	typedef Eigen::Matrix<float, k, k> A_t;

	void init(const popts_t &p) {
		// unpack
		a0 = p[0];
		a1 = p.segment<k>(1);
		// eigenUpperTriangularSet<4>(A2, p.segment</*  k *(k + 1) / 2  */10>(/* k + 1 */5));
		// eigenUpperTriangularSet<4>(A2, p.segment<10>(5));
			
		int kk = 0, N = 4;
		A2.setZero();
		// Upper triangular must be filled in row major to match python:
		// https://github.com/avikde/robobee3d/pull/173#issuecomment-674082069
		auto tvals = p.segment<k * (k + 1) / 2>(k + 1);
		for (int i = 0; i < N; ++i) {
			for (int j = i; j < N; ++j) {
				A2(i, j) = tvals[kk];
				kk++;
			}
		}
		// Use A2.selfadjointView<Eigen::Upper>() after this
	}

	inline float f(const x_t &xi) const {
		return a0 + xi.dot(a1) + 0.5 * xi.transpose() * A2.selfadjointView<Eigen::Upper>() * xi;
	}

	inline x_t Df(const x_t &xi) const {
		return a1 + A2.selfadjointView<Eigen::Upper>() * xi;
	}

protected:
	float a0;
	x_t a1;
	A_t A2;
};
