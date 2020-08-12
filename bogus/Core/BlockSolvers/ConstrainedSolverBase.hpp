/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_CONSTRAINED_SOLVER_BASE_HPP
#define BOGUS_CONSTRAINED_SOLVER_BASE_HPP

#include "BlockSolverBase.hpp"

namespace bogus {

template <typename Derived, typename BlockMatrixType>
class ConstrainedSolverBase : public BlockSolverBase<BlockMatrixType> {
  typedef BlockSolverBase<BlockMatrixType> Base;

 public:
  typedef typename Base::GlobalProblemTraits GlobalProblemTraits;
  typedef typename GlobalProblemTraits::Scalar Scalar;
  typedef typename BlockMatrixTraits<BlockMatrixType>::Index Index;

  typedef LocalProblemTraits<Base::BlockTraits::RowsPerBlock, Scalar>
      BlockProblemTraits;

  //! Sets whether the solver will use the infinity norm instead of the l1 one
  //! to compute the global residual from the local ones
  void useInfinityNorm(bool useInfNorm) { m_useInfinityNorm = useInfNorm; }
  bool usesInfinityNorm() const { return m_useInfinityNorm; }

  //! Eval the current global residual as a function of the local ones
  /*! \p y should be such that \p y = \ref m_matrix * \p x + rhs
          \return the current residual \c err defined as follow :
          - if \ref m_useInfinityNorm is true, then \c err \f$ := \max\limits_{1
     \leq i \leq n } law.eval(i,x_i,y_i) \f$
          - else \c err := \f$  \frac 1 {n+1} \sum\limits_{1 \leq i \leq n }
     law.eval(i,x_i,y_i) \f$

  */
  template <typename NSLaw, typename RhsT, typename ResT>
  Scalar eval(const NSLaw &law, const ResT &y, const RhsT &x) const;

  //! Projects the variable \p x on the constraints defined by \p projector
  template <typename NSLaw, typename VectorT>
  void projectOnConstraints(const NSLaw &projector, VectorT &x) const;

  //! Compute associated change of variable (see NSLaw)
  template <typename NSLaw, typename RhsT, typename ResT>
  void dualityCOV(const NSLaw &law, const RhsT &b, ResT &x) const;

  template <typename NSLaw, typename RhsT, typename ResT>
  Scalar solve(const NSLaw &law, const RhsT &b, ResT &x) const {
    return static_cast<const Derived &>(*this).solve(law, b, x);
  }

  //! Sets the system matrix and initializes internal structures
  /*! \note Derived classes should re-implement this function and call
   * updateScalings() */
  Derived &setMatrix(const BlockObjectBase<BlockMatrixType> &matrix) {
    return static_cast<Derived &>(*this).setMatrix(matrix);
  }

 protected:
  void updateScalings();

  ConstrainedSolverBase() : Base(), m_useInfinityNorm(false) {}
  using Base::m_matrix;

  typename GlobalProblemTraits::DynVector m_scaling;

  //! See useInfinityNorm(). Defaults to false.
  bool m_useInfinityNorm;
};

}  // namespace bogus

#endif
