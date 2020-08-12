/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_PROJECTED_GRADIENT_HPP
#define BOGUS_BLOCK_PROJECTED_GRADIENT_HPP

#include <vector>

#include "ConstrainedSolverBase.hpp"

namespace bogus {

//! Projected Gradient iterative solver.
template <typename BlockMatrixType>
class ProjectedGradient
    : public ConstrainedSolverBase<ProjectedGradient<BlockMatrixType>,
                                   BlockMatrixType> {
 public:
  typedef ConstrainedSolverBase<ProjectedGradient, BlockMatrixType> Base;
  typedef typename Base::Scalar Scalar;

  //! Default constructor -- you will have to call setMatrix() before using the
  //! solve() function
  ProjectedGradient() : Base() { init(); }
  //! Constructor with the system matrix
  explicit ProjectedGradient(const BlockObjectBase<BlockMatrixType> &matrix)
      : Base() {
    init();
    Base::setMatrix(matrix);
  }

  //! Finds an approximate minimum for a constrained quadratic problem
  template <typename NSLaw, typename RhsT, typename ResT>
  Scalar solve(const NSLaw &law, const RhsT &b, ResT &x) const;

  //! Finds an approximate minimum for a constrained quadratic problem
  /*!
   * Find an approximate solution to
   * \f[ \min_{x \in C}  \frac 1 2 x^T M x + x^T b \f]
   * where C is defined through the orthogonal projection operation on C
   * by NSLaw::projectOnConstraint().
   * \return an evalutation of the optimality error computed using NSLaw::eval()
   *
   * \tparam variant Which variant of the algorithm to use.
   */
  template <projected_gradient::Variant variant, typename NSLaw, typename RhsT,
            typename ResT>
  Scalar solve(const NSLaw &law, const RhsT &b, ResT &x) const;

  //! Sets the matrix M defining the quadratic objective function.
  /*! M should be symmetric, positive */
  ProjectedGradient &setMatrix(const BlockObjectBase<BlockMatrixType> &matrix) {
    m_matrix = &matrix;
    Base::updateScalings();
    return *this;
  }

  //! Sets the maximum number of line-search iterations
  void setLineSearchIterations(const unsigned lsIterations) {
    m_lsIters = lsIterations;
  }

  //! Sets the amount by which the step size will be multiplied at the beginninf
  //! of each PG iteration. Should be greater than 1
  void setLineSearchOptimisticFactor(const Scalar lsOptimisticFactor) {
    m_lsOptimisticFactor = lsOptimisticFactor;
  }

  //! Sets the amount by which the step size will be multiplied at the end of
  //! each line-search iterations. Should be in ]0,1[
  void setLineSearchPessimisticFactor(const Scalar lsPessimisticFactor) {
    m_lsPessimisticFactor = lsPessimisticFactor;
  }

  //! Sets the objective decrease coefficient for linesearchs that use an Armijo
  //! exit criterion Should be in ]0,1[
  void setLineSearchArmijoCoefficient(const Scalar lsArmijoCoefficient) {
    m_lsArmijoCoefficient = lsArmijoCoefficient;
  }

  //! Sets the variant that will be used when calling solve() without template
  //! arguments
  void setDefaultVariant(projected_gradient::Variant variant) {
    m_defaultVariant = variant;
  }

  unsigned lineSearchIterations() const { return m_lsIters; }
  Scalar lineSearchOptimisticFactor() const { return m_lsOptimisticFactor; }
  Scalar lineSearchPessimisticFactor() const { return m_lsPessimisticFactor; }
  Scalar lineSearchArmijoCoefficient() const { return m_lsArmijoCoefficient; }

 protected:
  typedef typename Base::Index Index;

  //! Sets up the default values for all parameters
  void init() {
    m_tol = 1.e-6;
    m_maxIters = 300;
    m_lsIters = 8;
    m_lsOptimisticFactor = 1.25;
    m_lsPessimisticFactor = .5;
    m_lsArmijoCoefficient = .5;
    m_defaultVariant = projected_gradient::APGD;
  }

  using Base::m_matrix;
  using Base::m_maxIters;
  using Base::m_tol;

  unsigned m_lsIters;
  Scalar m_lsOptimisticFactor;
  Scalar m_lsPessimisticFactor;
  Scalar m_lsArmijoCoefficient;

  projected_gradient::Variant m_defaultVariant;
};

}  // namespace bogus

#endif
