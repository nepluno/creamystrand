/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_EIGEN_MATRIX_TRAITS_HPP
#define BOGUS_EIGEN_MATRIX_TRAITS_HPP

#include <Eigen/Core>

#include "../Block/ScalarBindings.hpp"
#include "EigenLinearSolvers.hpp"
#include "EigenSparseLinearSolvers.hpp"

namespace bogus {

template <typename _MatrixType>
struct MatrixTraits {
  typedef _MatrixType MatrixType;
  typedef typename MatrixType::Scalar Scalar;

  typedef LU<Eigen::MatrixBase<MatrixType> > LUType;
  typedef LDLT<Eigen::MatrixBase<MatrixType> > LDLTType;

  static const MatrixType& asConstMatrix(const MatrixType& src) { return src; }
};

template <typename _Scalar, int _Options, typename _Index>
struct MatrixTraits<Eigen::SparseMatrix<_Scalar, _Options, _Index> > {
  typedef _Scalar Scalar;
  typedef Eigen::SparseMatrix<Scalar, _Options, _Index> MatrixType;

  typedef LU<
      Eigen::SparseMatrixBase<Eigen::SparseMatrix<Scalar, _Options, _Index> > >
      LUType;
  typedef LDLT<
      Eigen::SparseMatrixBase<Eigen::SparseMatrix<Scalar, _Options, _Index> > >
      LDLTType;

  static const MatrixType& asConstMatrix(const MatrixType& src) { return src; }
};

#define BOGUS_PROCESS_SCALAR(Scalar)                                     \
  template <>                                                            \
  struct MatrixTraits<Scalar>                                            \
      : public MatrixTraits<Eigen::Matrix<Scalar, 1, 1> > {              \
    static Eigen::Matrix<Scalar, 1, 1> asConstMatrix(const Scalar src) { \
      Eigen::Matrix<Scalar, 1, 1> mat;                                   \
      mat(0, 0) = src;                                                   \
      return mat;                                                        \
    }                                                                    \
  };
BOGUS_BLOCK_SCALAR_TYPES
#undef BOGUS_PROCESS_SCALAR

}  // namespace bogus

#endif
