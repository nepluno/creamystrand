/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_GAUSS_SEIDEL_HPP
#define BOGUS_BLOCK_GAUSS_SEIDEL_HPP

#include <vector>

#include "Coloring.hpp"
#include "GaussSeidelBase.hpp"

namespace bogus {

//! Projected Gauss-Seidel iterative solver.
/*!
   Works by taking into account only one block-row of the system at a time, and
   iterating several times over the whole set of rows several times until
   convergence has been achieved.

   Each inner iteration of the algorithm will try to solve the local problem
          \f[
                \left\{
                  \begin{array}{rcl}
                        y_i^{k+1} &=& M_{i,i} x_i^{k+1}  + b_i^{k} \\
                        &s.t.& law (x^{k+1},y^{k+1})
                  \end{array}
                \right.
          \f]
        where \b k is the current global iteration, \b i the current row
        and \f[ b_i^{k} := b_i + \sum_{ j < i }{ M_{i,j}x_j^{k+1} } +  \sum_{ j
   > i }{ M_{i,j}x_j^{k} } \f]

   See also solve() and \cite JAJ98.
  */
template <typename BlockMatrixType>
class GaussSeidel
    : public GaussSeidelBase<GaussSeidel<BlockMatrixType>, BlockMatrixType> {
 public:
  typedef GaussSeidelBase<GaussSeidel, BlockMatrixType> Base;

  typedef typename Base::GlobalProblemTraits GlobalProblemTraits;
  typedef typename GlobalProblemTraits::Scalar Scalar;

  //! Default constructor -- you will have to call setMatrix() before using the
  //! solve() function
  GaussSeidel() : Base() {}
  //! Constructor with the system matrix
  explicit GaussSeidel(const BlockObjectBase<BlockMatrixType> &matrix)
      : Base() {
    setMatrix(matrix);
  }

  //! Sets the system matrix and initializes internal structures
  GaussSeidel &setMatrix(const BlockObjectBase<BlockMatrixType> &matrix);

  //! Finds an approximate solution for a constrained linear problem
  /*!
    Stops when the residual computed in eval() is below \ref m_tol, of the
    number of iterations exceed \ref m_maxIters

    Implements Algorithm 1. from \cite DBB11 to solve
     \f[
          \left\{
            \begin{array}{rcl}
                  y &=& M x + b \\
                  &s.t.& law (x,y)
            \end{array}
          \right.
    \f]
    \param law The (non-smooth) law that should define:
          - An error function for the local problem
          - A local solver for each row of the system ( e.g. 1 contact solver )
          \sa SOCLaw
    \param b the const part of the right hand side
    \param x the unknown. Can be warm-started
    \param tryZeroAsWell If true, the algorithm will reset r to zero if that
    would result in a lower residual
    */
  template <typename NSLaw, typename RhsT, typename ResT>
  Scalar solve(const NSLaw &law, const RhsT &b, ResT &x,
               bool tryZeroAsWell) const;

  /*!
     Solves
     \f[
          \left\{
            \begin{array}{rclll}
                  y &=& M x   &+ p &+ b \\
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
                  y &=& M x + W x + b \\
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

  //! Access to the current Coloring. Will be reset whenever the matrix is
  //! changed.
  /*! Determiniticy is achieved through the mean of contact coloring ;
    contacts that do not interact directly together can chare the same color,
    and all contacts within a given color can be solver in parallel */
  Coloring &coloring() { return m_coloring; }

  using Base::solve;

 protected:
  void updateLocalMatrices();

  template <typename NSLaw, typename RhsT, typename ResT>
  void innerLoop(bool parallelize, const NSLaw &law, const RhsT &b,
                 std::vector<unsigned char> &skip, Scalar &ndxRef,
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

  //! \sa coloring()
  Coloring m_coloring;
};

}  // namespace bogus

#endif
