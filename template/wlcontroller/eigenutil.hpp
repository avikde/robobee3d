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
typedef Eigen::Matrix<float, 6, 1> w_t;
typedef Eigen::Matrix<float, 6, 4> dw_du_t;
typedef Eigen::Map<const u_t> mc_u;
typedef Eigen::Array<float, Eigen::Dynamic, 1> ArrX_t;
typedef Eigen::Matrix<float, 6, 6> M_t;

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

/**
 * @brief 
 * 
 * @tparam N matrix number of rows or cols
 * @param M matrix
 * @param tvals Must have size N * (N + 1) / 2
 */
template< int N >
inline void eigenUpperTriangularVals(const Eigen::Ref<const Eigen::Matrix<float, N, N> > & M, Eigen::VectorXf &tvals)
{
  tvals.resize(4 * (4+1)/2);
  int k = 0;
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i <= j; ++i) {
      tvals[k] = M(i, j);
      k++;
    }
  }
}

// } // namespace ctrl
