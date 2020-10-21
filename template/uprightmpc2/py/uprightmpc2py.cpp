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

class UprightMPC2 {
public:
	UprightMPC_t umpc;

	UprightMPC2(float dt, float g, const Eigen::Vector3f &smin, const Eigen::Vector3f &smax, float TtoWmax, float ws, float wds, float wpr, float wpf, float wvr, float wvf, float wthrust, float wmom) {
		umpcInit(&umpc, dt, g, smin.data(), smax.data(), TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom);
	}

	uacc_t update(const Eigen::Vector3f &p0, const Eigen::Matrix3f &R0, const Vec6_t &dq0, const Eigen::Vector3f &pdes, const Eigen::Vector3f &dpdes) {
		static Eigen::Vector3f uquad;
		static Vec6_t accdes;
		umpcUpdate(&umpc, uquad.data(), accdes.data(), p0.data(), dq0.data(), R0.data(), pdes.data(), dpdes.data());
		return std::make_tuple(uquad, accdes);
	}
};

PYBIND11_MODULE(uprightmpc2py, m) {
	py::class_<UprightMPC2>(m, "UprightMPC2C")
	.def(py::init<float /* dt */, float /* g */, const Eigen::Vector3f &/* smin */, const Eigen::Vector3f &/* smax */, float /* TtoWmax */, float /* ws */, float /* wds */, float /* wpr */, float /* wpf */, float /* wvr */, float /* wvf */, float /* wthrust */, float /* wmom */>())
	.def("update", &UprightMPC2::update);
}
