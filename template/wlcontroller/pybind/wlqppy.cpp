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
#include "../cpp/wlcontroller.hpp"
namespace py = pybind11;

PYBIND11_MODULE(wlqppy, m) {
  py::class_<WLController>(m, "WLController")
  .def(py::init<>())
  .def("update", &WLController::update, "u = update(u0, h0, pdotdes)");
}
