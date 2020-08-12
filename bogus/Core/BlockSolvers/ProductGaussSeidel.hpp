/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_PRODUCT_GAUSS_SEIDEL_HPP
#define BOGUS_PRODUCT_GAUSS_SEIDEL_HPP

#include <vector>

#include "GaussSeidelBase.hpp"
#include "ProductGaussSeidelUtils.hpp"

namespace bogus {

//! Matrix-free version of the GaussSeidel iterative solver.
/*!
  Assumes that the system matrix is defined as the product (M D M'),
  with D a block-diagonal matrix whose block sizes coincide with those of
  the columns of M.

  \tparam BlockMatrixType The type of the main solver matrix M (the
  constructor's first argument). \tparam DiagonalType The type of the diagonal
  matrix D. The default type (M's Scalar type) means using a constant times the
  identity matrix. \tparam PrecompupeDMt Whether to precompute and store the (D
  M') part of the product (M D M'). Defaults to true if DiagonalType is not the
  default value.

  \warning Parallelization is supported, but dangerous. If in doubt, use
  setMaxThreads(1) \note D must be block-diagonal. Works best when M is
  row-major and its columns are quite sparse \sa GaussSeidel
  */
template <typename BlockMatrixType,
          typename DiagonalType = typename BlockMatrixType::Scalar,
          bool PrecomputeDMt =
              !(IsSame<DiagonalType, typename BlockMatrixType::Scalar>::Value)>
class ProductGaussSeidel
    : public GaussSeidelBase<
          ProductGaussSeidel<BlockMatrixType, DiagonalType, PrecomputeDMt>,
          BlockMatrixType> {
 public:
  typedef GaussSeidelBase<ProductGaussSeidel, BlockMatrixType> Base;

  typedef typename Base::GlobalProblemTraits GlobalProblemTraits;
  typedef typename GlobalProblemTraits::Scalar Scalar;

  enum {
    has_trivial_diagonal =
        IsSame<DiagonalType, typename BlockMatrixType::Scalar>::Value
  };
  typedef block_solvers_impl::DiagonalMatrixWrapper<DiagonalType,
                                                    !!has_trivial_diagonal>
      DiagWrapper;
  typedef block_solvers_impl::DMtStorage<BlockMatrixType, DiagWrapper,
                                         PrecomputeDMt>
      DMtStorage;

  //! Default constructor -- you will have to call setMatrix() before using the
  //! solve() function
  ProductGaussSeidel() : Base() {}
  //! Constructor with the main system matrix M
  explicit ProductGaussSeidel(const BlockObjectBase<BlockMatrixType> &matrix)
      : Base() {
    setMatrix(matrix);
  }
  //! Constructor with both the main matrix M and the diagonal D
  ProductGaussSeidel(const BlockObjectBase<BlockMatrixType> &matrix,
                     const DiagonalType &diagonal)
      : Base(), m_diagonal(DiagWrapper(diagonal)) {
    setMatrix(matrix);
  }

  //! Sets the system matrix (M) and initializes internal structures
  ProductGaussSeidel &setMatrix(const BlockObjectBase<BlockMatrixType> &matrix);

  //! Sets the system diagonal (D) and initializes internal structures
  ProductGaussSeidel &setDiagonal(const DiagonalType &diagonal);

  //! Finds an approximate solution for a constrained linear problem
  template <typename NSLaw, typename RhsT, typename ResT>
  Scalar solve(const NSLaw &law, const RhsT &b, ResT &x,
               bool tryZeroAsWell) const;

  /*!
     Solves
     \f[
          \left\{
            \begin{array}{rclll}
                  y &=& MDM^T x   &+ p &+ b \\
                  0 &=& H^T x &+ \frac{1}{\alpha}C p &+ c \\
                  &s.t.& law (x,y)
            \end{array}
          \right.
    \f]
     where \p Cinv is such that \p Cinv * x =  \f$ C^{-1} x \f$

     \warning Requires: m_evalEvery multiple of solveEvery ;
   */
  template <typename NSLaw, typename RhsT, typename ResT, typename LSDerived,
            typename HDerived>
  Scalar solveWithLinearConstraints(const NSLaw &law,
                                    const BlockObjectBase<LSDerived> &Cinv,
                                    const BlockObjectBase<HDerived> &H,
                                    const Scalar alpha, const RhsT &b,
                                    const RhsT &c, ResT &x,
                                    bool tryZeroAsWell = true,
                                    unsigned solveEvery = 1) const;

  /*!
     Solves
     \f[
          \left\{
            \begin{array}{rcl}
                  y &=& MDM^T x + W x + b \\
                  &s.t.& law (x,y)
            \end{array}
          \right.
    \f]
     with W arbitrary linear operator ( matrix or expression )
     \warning Requires: m_evalEvery multiple of solveEvery ;
   */
  template <typename NSLaw, typename RhsT, typename ResT, typename WDerived>
  Scalar solveWithLinearConstraints(const NSLaw &law,
                                    const BlockObjectBase<WDerived> &W,
                                    const RhsT &b, ResT &x,
                                    bool tryZeroAsWell = true,
                                    unsigned solveEvery = 1) const;

  using Base::solve;

 protected:
  DiagWrapper m_diagonal;
  DMtStorage m_DMt;

  void updateLocalMatrices();

  template <typename NSLaw, typename VecT, typename ResT>
  void innerLoop(bool parallelize, const NSLaw &law, const VecT &b,
                 std::vector<unsigned char> &skip, Scalar &ndxRef, VecT &Mx,
                 ResT &x) const;

  typedef typename Base::Index Index;

  using Base::m_evalEvery;
  using Base::m_localMatrices;
  using Base::m_matrix;
  using Base::m_maxIters;
  using Base::m_maxThreads;
  using Base::m_regularization;
  using Base::m_scaling;
  using Base::m_skipIters;
  using Base::m_skipTol;
  using Base::m_tol;
};

}  // namespace bogus

#endif
