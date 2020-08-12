/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_GAUSS_SEIDEL_BASE_HPP
#define BOGUS_BLOCK_GAUSS_SEIDEL_BASE_HPP

#include <vector>

#include "Coloring.hpp"
#include "ConstrainedSolverBase.hpp"

namespace bogus {

//! Abstract Gauss-Seidel interface . \sa GaussSeidel
template <typename GaussSeidelImpl, typename BlockMatrixType>
class GaussSeidelBase
    : public ConstrainedSolverBase<GaussSeidelImpl, BlockMatrixType> {
 public:
  typedef ConstrainedSolverBase<GaussSeidelImpl, BlockMatrixType> Base;

  typedef typename Base::GlobalProblemTraits GlobalProblemTraits;
  typedef typename GlobalProblemTraits::Scalar Scalar;

  //! Sets the maximum number of threads that the solver can use.
  /*! If \p maxThreads is zero, then it will use the current OpenMP setting.

  \warning If multi-threading is enabled without coloring,
          the result will not be deterministic, as it will depends on the
          order in which threads solve contacts.

          On the other hand, the algorithm will run much faster.
  */
  void setMaxThreads(unsigned maxThreads = 0) { m_maxThreads = maxThreads; }

  //! Sets the auto-regularization (a.k.a. proximal point) coefficient
  /*!
          The regularization works by slightly altering the local problems, so
     at each iteration we try to solve \f[ \left\{ \begin{array}{rcl}
                  y^{k+1} &=& ( M + \alpha I ) x^{k+1} - \alpha x^k + b^{k} \\
                  &s.t.& law (x^{k+1},y^{k+1})
                  \end{array}
          \right.
          \f]
          where \f$\alpha\f$ is the regularization coefficient.

          Note that as \f$ | x^{k+1} - x^{k} | \rightarrow 0 \f$ when the
     algorithm converges, we are still trying to find a solution of the same
     global problem.

          For under-determined problems, regularization might helps preventing
     \b x reaching problematically high values. Setting \f$\alpha\f$ to a too
     big value will however degrade the convergence of the global algorithm.

          \param maxRegul If greater than zero, then positive terms will be
     added to the diagonal of the local matrices so that their minimal
     eigenvalue become greater than \p maxRegul.
          */
  void setAutoRegularization(Scalar maxRegul) {
    m_autoRegularization = maxRegul;
  }

  // Debug

  //! Sets the number of iterations that should be performed between successive
  //! evaluations of the global error function
  /*! ( Those evaluations require a full matrix/vector product, and are therfore
   * quite costly ) */
  void setEvalEvery(unsigned evalEvery) { m_evalEvery = evalEvery; }
  //! Sets the minimum iteration step size under which local problems are
  //! temporarily frozen
  void setSkipTol(Scalar skipTol) { m_skipTol = skipTol; }
  //! Sets the number of iterations for temporarily freezing local problems
  void setSkipIters(unsigned skipIters) { m_skipIters = skipIters; }

  //! If true, the algorithm will reset r to zero if that would result in a
  //! lower residual ( unless otherwise specified by an explicit argument to the
  //! solve function )
  void doTryZeroAsWell(bool doTry = true) { m_tryZeroAsWell = doTry; }

  template <typename NSLaw, typename RhsT, typename ResT>
  Scalar solve(const NSLaw &law, const RhsT &b, ResT &x) const {
    return static_cast<const GaussSeidelImpl &>(*this).solve(law, b, x,
                                                             m_tryZeroAsWell);
  }

 protected:
  GaussSeidelBase()
      : m_maxThreads(0),
        m_evalEvery(25),
        m_skipIters(10),
        m_autoRegularization(0.),
        m_tryZeroAsWell(true) {
    m_tol = 1.e-6;
    m_maxIters = 250;
    m_skipTol = 1.e-6;
  }

  const BlockMatrixBase<BlockMatrixType> &explicitMatrix() const {
    // Will cause a compile error if BlockMatrixType does not derive from
    // BlockMatrixBase This is voluntary, GaussSeidel implementations does not
    // handle arbitrary expressions yet
    return m_matrix->derived();
  }
  const IterableBlockObject<BlockMatrixType> &iterableMatrix() const {
    // Will cause a compile error if BlockMatrixType does not derive from
    // IterableObjectBase This is voluntary, GaussSeidel implementations does
    // not handle arbitrary expressions yet
    return m_matrix->derived();
  }

  // Add regularization, compute scalings from diagonal block matrices
  void processLocalMatrices();

  template <typename NSLaw, typename ResT>
  Scalar evalAndKeepBest(const NSLaw &law, const ResT &x,
                         const typename GlobalProblemTraits::DynVector &y,
                         typename GlobalProblemTraits::DynVector &x_best,
                         Scalar &err_best) const;

  template <typename NSLaw, typename RhsT, typename ResT>
  bool tryZero(const NSLaw &law, const RhsT &b, ResT &x,
               typename GlobalProblemTraits::DynVector &x_best,
               Scalar &err_best) const;

  using Base::m_matrix;
  using Base::m_maxIters;
  using Base::m_scaling;
  using Base::m_tol;

  typedef typename Base::Index Index;
  typedef typename Base::BlockProblemTraits::Matrix DiagonalBlockType;
  typename ResizableSequenceContainer<DiagonalBlockType>::Type m_localMatrices;
  typename GlobalProblemTraits::DynVector m_regularization;

  //! See setMaxThreads(). Defaults to 0 .
  unsigned m_maxThreads;

  //! See setEvalEvery(). Defaults to 25
  unsigned m_evalEvery;
  //! See setSkipTol(). Defaults to 1.e-6
  Scalar m_skipTol;
  //! See setSkipIters() Defaults to 10
  unsigned m_skipIters;

  //! \sa setAutoRegularization(). Defaults to 0.
  Scalar m_autoRegularization;

  //! \sa doTryZeroAsWell(). Defaults to true.
  bool m_tryZeroAsWell;
};

}  // namespace bogus

#endif
