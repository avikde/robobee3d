/**
 * @file ematmult.hpp
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-11-18
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once
#define EIGEN_NO_DEBUG
#include <Eigen/Core>

typedef Eigen::Map<const Eigen::MatrixXf> MCMatX;
typedef Eigen::Map<Eigen::MatrixXf> MMatX;
typedef Eigen::Matrix<float, 6, 1> Vec6_t;
