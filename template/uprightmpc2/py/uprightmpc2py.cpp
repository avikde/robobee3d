/**
 * @file wlqpc.cpp
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-10-21
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <tuple>
#include <uprightmpc2.h>
#include <funapprox.h>
#include "ematmult.hpp"

namespace py = pybind11;

typedef std::tuple<Eigen::Vector3f, Vec6_t> ret_t;
typedef Eigen::Matrix<float, UMPC_NX, 1> Vecx_t;
typedef Eigen::Matrix<float, UMPC_NC, 1> Vecc_t;
typedef Eigen::Matrix<float, UMPC_nAdata, 1> Adatax_t;
typedef Eigen::Array<int, UMPC_nAdata, 1> Adatai_t;
typedef std::tuple<Vecc_t, Vecc_t, Vecx_t> uvecs_t;
typedef std::tuple<Vecx_t, Adatax_t, Adatai_t> umats_t;
typedef std::tuple<Eigen::Vector4f, Vec6_t> wlret_t;

// Wrapper for the C implementation
class UprightMPC2 {
public:
	UprightMPC_t umpc;

	UprightMPC2(float dt, float g, float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, const Eigen::Vector3f &Ib, int maxIter) {
		umpcInit(&umpc, dt, g, TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom, Ib.data(), maxIter);
	}

	ret_t update(const Eigen::Vector3f &p0, const Eigen::Matrix3f &R0, const Vec6_t &dq0, const Eigen::Vector3f &pdes, const Eigen::Vector3f &dpdes, const Eigen::Vector3f &sdes, float actualT0) {
		static Eigen::Vector3f uquad;
		static Vec6_t accdes;
		umpcUpdate(&umpc, uquad.data(), accdes.data(),p0.data(), R0.data(), dq0.data(), pdes.data(), dpdes.data(), sdes.data(), actualT0);
		return std::make_tuple(uquad, accdes);
	}

	// for debugging
	uvecs_t vectors() {
		return std::make_tuple(Eigen::Map<Vecc_t>(umpc.l), Eigen::Map<Vecc_t>(umpc.u), Eigen::Map<Vecx_t>(umpc.q));
	}
	umats_t matrices() {
		return std::make_tuple(Eigen::Map<Vecx_t>(umpc.Px_data), Eigen::Map<Adatax_t>(umpc.Ax_data), Eigen::Map<Adatai_t>(umpc.Ax_idx));
	}
};

class WLCon {
public:
	WLCon_t wl;
	
	WLCon(const Eigen::Vector4f &u0, const Eigen::Vector4f &umin, const Eigen::Vector4f &umax, const Eigen::Vector4f &dumax, const Vec6_t &Qw, float controlRate, const Eigen::Matrix<float, 90, 1> &popts) {
		wlConInit(&wl, u0.data(), umin.data(), umax.data(), dumax.data(), Qw.data(), controlRate, popts.data());
	}

	wlret_t update(const Vec6_t &h0, const Vec6_t &pdotdes) {
		static Eigen::Vector4f u1;
		static Vec6_t w0;
		wlConUpdate(&wl, u1.data(), w0.data(), h0.data(), pdotdes.data());
		return std::make_tuple(u1, w0);
	}
};

PYBIND11_MODULE(uprightmpc2py, m) {
	py::class_<UprightMPC2>(m, "UprightMPC2C")
	.def(py::init<float /* dt */, float /* g */, float /* TtoWmax */, float /* ws */, float /* wds */, float /* wpr */, float /* wpf */, float /* wvr */, float /* wvf */, float /* wthrust */, float /* wmom */, const Eigen::Vector3f &/* Ib */, int /* maxIter */>())
	.def("update", &UprightMPC2::update)
	.def("vectors", &UprightMPC2::vectors)
	.def("matrices", &UprightMPC2::matrices);
	
	py::class_<WLCon>(m, "WLCon")
	.def(py::init<const Eigen::Vector4f &/* u0 */, const Eigen::Vector4f &/* umin */, const Eigen::Vector4f &/* umax */, const Eigen::Vector4f &/* dumax */, const Vec6_t &/* Qw */, float /* controlRate */, const Eigen::Matrix<float, 90, 1> &/* popts */>())
	.def("update", &WLCon::update);
}
