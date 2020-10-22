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
#define EIGEN_NO_DEBUG
#include <Eigen/Core>

namespace py = pybind11;

typedef Eigen::Matrix<float, 6, 1> Vec6_t;
typedef std::tuple<Eigen::Vector3f, Vec6_t> uacc_t;
typedef Eigen::Map<const Eigen::MatrixXf> MCMatX;
typedef Eigen::Map<Eigen::MatrixXf> MMatX;
typedef Eigen::Matrix<float, UMPC_NX, 1> Vecx_t;
typedef Eigen::Matrix<float, UMPC_NC, 1> Vecc_t;
typedef std::tuple<Vecc_t, Vecc_t, Vecx_t> uvecs_t;

// Wrapper for the C implementation
class UprightMPC2 {
public:
	UprightMPC_t umpc;

	UprightMPC2(float dt, float g, const Eigen::Vector3f &smin, const Eigen::Vector3f &smax, float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, const Eigen::Vector3f &Ib, int maxIter) {
		umpcInit(&umpc, dt, g, smin.data(), smax.data(), TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom, Ib.data(), maxIter);
	}

	uacc_t update(const Eigen::Vector3f &p0, const Eigen::Matrix3f &R0, const Vec6_t &dq0, const Eigen::Vector3f &pdes, const Eigen::Vector3f &dpdes) {
		static Eigen::Vector3f uquad;
		static Vec6_t accdes;
		umpcUpdate(&umpc, uquad.data(), accdes.data(), p0.data(), R0.data(), dq0.data(), pdes.data(), dpdes.data());
		return std::make_tuple(uquad, accdes);
	}

	// for debugging
	uvecs_t vectors() {
		return std::make_tuple(Eigen::Map<Vecc_t>(umpc.l), Eigen::Map<Vecc_t>(umpc.u), Eigen::Map<Vecx_t>(umpc.q));
	}
};

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

PYBIND11_MODULE(uprightmpc2py, m) {
	py::class_<UprightMPC2>(m, "UprightMPC2C")
	.def(py::init<float /* dt */, float /* g */, const Eigen::Vector3f &/* smin */, const Eigen::Vector3f &/* smax */, float /* TtoWmax */, float /* ws */, float /* wds */, float /* wpr */, float /* wpf */, float /* wvr */, float /* wvf */, float /* wthrust */, float /* wmom */, const Eigen::Vector3f &, int>())
	.def("update", &UprightMPC2::update)
	.def("vectors", &UprightMPC2::vectors);
}