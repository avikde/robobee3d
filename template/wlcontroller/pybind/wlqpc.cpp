/**
 * @file wlqpc.cpp
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-29
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <wlcontroller.h>

typedef Eigen::Matrix<float, 4, 1> u_t;
typedef Eigen::Matrix<float, 6, 1> w_t;

u_t wlControllerWrap(const u_t &u0init, const w_t &h0, const w_t &pdotdes) {
	static u_t u;
	wlControllerUpdate(u.data(), u0init.data(), h0.data(), pdotdes.data());
	return u;
}

PYBIND11_MODULE(wlqpc, m) {
	m.doc() = "Python binding of C version";
	m.def("update", &wlControllerWrap, "wlControllerUpdate");
}
