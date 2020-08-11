/**
 * @file eigenutil.hpp
 * @author Avik De (avikde@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-08-11
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once

#define EIGEN_NO_DEBUG
#include <Eigen/Core>
// #include <Eigen/SparseCore>
// #include <osqp.h>
// #include <vector>

typedef Eigen::Vector4f u_t;
typedef Eigen::Map<const u_t> mc_u;

// namespace ctrl
// {

// // use c_float (which is float or double depending on cmake option DFLOAT)
// typedef Eigen::Triplet<c_float> Trip_t;
// typedef Eigen::SparseMatrix<c_float> SpMat_t;
// typedef Eigen::Matrix<c_float, Eigen::Dynamic, Eigen::Dynamic> MatX_t;
// typedef Eigen::Matrix<c_float, 2, 2> Mat2_t;
// typedef Eigen::Matrix<c_float, 3, 3> Mat3_t;
// typedef Eigen::Matrix<c_float, 4, 4> Mat4_t;
// typedef Eigen::Matrix<c_float, 6, 6> Mat6_t;
// typedef Eigen::Matrix<c_float, Eigen::Dynamic, 1> VecX_t;
// typedef Eigen::Matrix<c_float, 1, 1> Vec1_t;      // scalar for params
// typedef Eigen::Matrix<c_float, 2, 1> Vec2_t;
// typedef Eigen::Matrix<c_float, 3, 1> Vec3_t;
// typedef Eigen::Matrix<c_float, 4, 1> Vec4_t;
// typedef Eigen::Matrix<c_float, 6, 1> Vec6_t;
// typedef Eigen::Array<c_float, Eigen::Dynamic, 1> ArrX_t;
// typedef Eigen::AngleAxis<c_float> AngleAxis_t;

// inline c_float clip(c_float a, c_float l = 0.0f, c_float u = 1.0f)
// {
//   return (a > u) ? u : (a < l ? l : a);
// }

// /**
//  * @brief Convert from Eigen sparse to OSQP
//  *
//  * @param sm Eigen sparse
//  * @return csc* pass to data->P or A
//  */
// inline csc * eigenToCsc(const SpMat_t & sm)
// {
//   return csc_matrix((c_int)sm.rows(), (c_int)sm.cols(), (c_int)sm.nonZeros(), (c_float *)sm.valuePtr(), (c_int *)sm.innerIndexPtr(),
//            (c_int *)sm.outerIndexPtr());
// }

// /**
//  * @brief Given a dense matrix, return the upper-triangular elements. There must be a better way to do this within Eigen.
//  *
//  * M.triangularView<Eigen::Upper>()
//  *
//  * returns the triangular view but have not found a way to get the elements in order
//  *
//  * @param M Matrix assumed square
//  * @return ArrX_t
//  */
// inline ArrX_t eigenUpperTriangularVals(const Eigen::Ref<const MatX_t> & M)
// {
//   static ArrX_t tvals;
//   tvals.resize(M.rows() * (M.cols() + 1) / 2);
//   int k = 0;
//   for (int j = 0; j < M.rows(); ++j) {
//     for (int i = 0; i <= j; ++i) {
//       tvals[k] = M(i, j);
//       k++;
//     }
//   }
//   return tvals;
// }

// } // namespace ctrl
