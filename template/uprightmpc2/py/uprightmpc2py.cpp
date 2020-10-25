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
typedef std::tuple<Eigen::Vector3f, Vec6_t, Eigen::Vector4f> ret_t;
typedef Eigen::Map<const Eigen::MatrixXf> MCMatX;
typedef Eigen::Map<Eigen::MatrixXf> MMatX;
typedef Eigen::Matrix<float, UMPC_NX, 1> Vecx_t;
typedef Eigen::Matrix<float, UMPC_NC, 1> Vecc_t;
typedef Eigen::Matrix<float, UMPC_nPdata, 1> Pdatax_t;
typedef Eigen::Matrix<float, UMPC_nAdata, 1> Adatax_t;
typedef Eigen::Array<int, UMPC_nAdata, 1> Adatai_t;
typedef std::tuple<Vecc_t, Vecc_t, Vecx_t> uvecs_t;
typedef std::tuple<Pdatax_t, Adatax_t, Adatai_t> umats_t;

// Wrapper for the C implementation
class UprightMPC2 {
public:
	UprightMPC_t umpc;

	UprightMPC2(float dt, float g, float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom, float mb, const Eigen::Vector3f &Ib, const Eigen::Vector4f &umin, const Eigen::Vector4f &umax, const Eigen::Vector4f &dumax, const Vec6_t &Qw, float controlRate, int maxIter, const Eigen::Matrix<float, 90, 1> &popts) {
		umpcInit(&umpc, dt, g, TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom, mb, Ib.data(), umin.data(), umax.data(), dumax.data(), Qw.data(), controlRate, maxIter, popts.data());
	}

	ret_t update(const Eigen::Vector3f &p0, const Eigen::Matrix3f &R0, const Vec6_t &dq0, const Eigen::Vector3f &pdes, const Eigen::Vector3f &dpdes) {
		static Eigen::Vector3f uquad;
		static Eigen::Vector4f uwlqp;
		static Vec6_t accdes;
		umpcUpdate(&umpc, uquad.data(), accdes.data(), uwlqp.data(), p0.data(), R0.data(), dq0.data(), pdes.data(), dpdes.data());
		return std::make_tuple(uquad, accdes, uwlqp);
	}

	// for debugging
	uvecs_t vectors() {
		return std::make_tuple(Eigen::Map<Vecc_t>(umpc.l), Eigen::Map<Vecc_t>(umpc.u), Eigen::Map<Vecx_t>(umpc.q));
	}
	umats_t matrices() {
		return std::make_tuple(Eigen::Map<Pdatax_t>(umpc.Px_data), Eigen::Map<Adatax_t>(umpc.Ax_data), Eigen::Map<Adatai_t>(umpc.Ax_idx));
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
	.def(py::init<float /* dt */, float /* g */, float /* TtoWmax */, float /* ws */, float /* wds */, float /* wpr */, float /* wpf */, float /* wvr */, float /* wvf */, float /* wthrust */, float /* wmom */, float /* mb */, const Eigen::Vector3f &/* Ib */, const Eigen::Vector4f &/* umin */, const Eigen::Vector4f &/* umax */, const Eigen::Vector4f &/* dumax */, const Vec6_t &/* Qw */, float, int /* maxIter */, const Eigen::Matrix<float, 90, 1> &/* popts */>())
	.def("update", &UprightMPC2::update)
	.def("vectors", &UprightMPC2::vectors)
	.def("matrices", &UprightMPC2::matrices);
}
