/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_ADMM_HPP
#define BOGUS_ADMM_HPP

#include <vector>

#include "ConstrainedSolverBase.hpp"

namespace bogus {

//! Options for ADMM solvers
namespace admm {
//! Variants of ADMM algorithm
enum Variant {
  //! Standard ADMM
  Standard,
  //! ADDM with \cite Nesterov1983 acceleration (see Goldstein 2014)
  Accelerated
};
}  // namespace admm

//! ADMM (Alternating Direction Method of Multipliers ) iterative solver.
/*!
        Minimizes J(v) with ( M v + w ) in C, C defined by the NSLaw in the
   solve() function.

        Requires ability to evaluate
   \f$ prox_{\lambda} J( x ) := \min_y J(y) + \frac 1 {2 \lambda} \Vert x - y
   \Vert^2  \f$

        \sa QuadraticProxOp
        \warning Experimental
*/
template <typename BlockMatrixType>
class ADMM
    : public ConstrainedSolverBase<ADMM<BlockMatrixType>, BlockMatrixType> {
 public:
  typedef ConstrainedSolverBase<ADMM, BlockMatrixType> Base;

  typedef typename Base::GlobalProblemTraits GlobalProblemTraits;
  typedef typename GlobalProblemTraits::Scalar Scalar;

  //! Default constructor -- you will have to call setMatrix() before using the
  //! solve() function
  ADMM() : Base() { init(); }
  //! Constructor with the system matrix
  explicit ADMM(const BlockObjectBase<BlockMatrixType> &matrix) : Base() {
    init();
    Base::setMatrix(matrix);
  }

  /*!
     Find the minimizer of J(v),  ( Mv + w in C )
     with J(v) defined through its proximal operator, \p op
     and C defined by \p law
     \param law should define projectOnConstraints()
     \param op should define a \c eval(Rhs& rhs, Res& res) function computing \p
     res as \f$ prox_{J, 1/\lambda}( rhs/\lambda ) \f$ and \c coefficient()
     returning \f$\lambda\f$

     \sa QuadraticProxOp
   */
  template <admm::Variant variant, typename NSLaw, typename ProxOp,
            typename RhsT, typename ResT>
  Scalar solve(const NSLaw &law, const ProxOp &op, const RhsT &w, ResT &v,
               ResT &r) const;

  //! Solve function using default variant
  template <typename NSLaw, typename ProxOp, typename RhsT, typename ResT>
  Scalar solve(const NSLaw &law, const ProxOp &op, const RhsT &w, ResT &v,
               ResT &r) const;

  //! Sets the problem matrix -- the one defining the constraints
  ADMM &setMatrix(const BlockObjectBase<BlockMatrixType> &matrix) {
    m_matrix = &matrix;
    Base::updateScalings();
    return *this;
  }

  //! Sets the step size for updating the dual variable (forces).
  void setStepSize(const Scalar size) { m_stepSize = size; }

  //! Sets the variant that will be used when calling solve() without template
  //! arguments
  void setDefaultVariant(admm::Variant variant) { m_defaultVariant = variant; }

  Scalar stepSize() const { return m_stepSize; }

 protected:
  typedef typename Base::Index Index;

  //! Sets up the default values for all parameters
  void init() {
    m_tol = 1.e-6;
    m_maxIters = 300;

    m_stepSize = 1.e-2;

    m_defaultVariant = admm::Accelerated;
  }

  using Base::m_matrix;
  using Base::m_maxIters;
  using Base::m_tol;

  Scalar m_stepSize;

  admm::Variant m_defaultVariant;
};

//! Evaluation of prox_{1/c} J  with J(x) = 1/2 xM'x + f'x
/*! Requires a way to solve linear systems with LHS ( M + c I)
 * \tparam ObjectType the type of the expression representing the inverse of (M
 * + c I)
 */
template <typename ObjectType>
struct QuadraticProxOp {
  typedef BlockMatrixTraits<ObjectType> BlockTraits;
  typedef typename BlockTraits::Scalar Scalar;
  typedef ProblemTraits<Scalar> GlobalProblemTraits;

  typedef BlockObjectBase<ObjectType> LinearOp;
  typedef typename GlobalProblemTraits::DynVector AffineVec;

  //! Construct the proximal evaluator of a quadratic function
  /*!
    \param lambda the inverse of the proximal function coefficient
    \param linearExpr Mutliplying by \p linearExpr should be equivalent to
    solving a linear system with LHS ( M + \p lambda * I ) \param affinePart The
    term f in the quadratic function J
  */
  //!
  QuadraticProxOp(const LinearOp &linearExpr, const Scalar lambda,
                  const AffineVec &affinePart)
      : m_coefficient(lambda),
        m_linearOp(linearExpr),
        m_affinePart(affinePart) {}

  //! Evaluates prox_{1/coeff} J ( rhs/coeff )
  /*! \warning may modify rhs */
  template <typename RhsT, typename ResT>
  void eval(RhsT &rhs, ResT &res) const;

  //! \returns coeff so that eval() computes prox_{1/coeff}
  Scalar coefficient() const { return m_coefficient; }

 private:
  const Scalar m_coefficient;
  const LinearOp &m_linearOp;
  const AffineVec &m_affinePart;
};

//! Dual AMA iterative solver (Alternating Minimization Algorithm  on dual
//! formuation of quadratic optimization problem).
/*!

        \verbatim
        Minimizes .5 r' M A^-1 M' r + r' ( w - M A^-1 f ) for r in C

        as min H(x) + G(r)  , x = M' r
        with H(x) = .5 x' A^-1 x - x' A^-1 f
                 G(r) = I_C(r) + r'w

        thanks to the identities
         prox_{l,G} (y) = prox_{l,Ic}( y - lw ) = Pi_C ( y - lw )
        and
         inf_x  H(x) - < v, x > = A v + f
        \endverbatim

        Main advantage of this solver is that it does not require inverting or
        solving linear system with the stiffness matrix A.

        \warning Experimental

*/
template <typename BlockMatrixType>
class DualAMA
    : public ConstrainedSolverBase<DualAMA<BlockMatrixType>, BlockMatrixType> {
 public:
  typedef ConstrainedSolverBase<DualAMA, BlockMatrixType> Base;

  typedef typename Base::GlobalProblemTraits GlobalProblemTraits;
  typedef typename GlobalProblemTraits::Scalar Scalar;

  //! Default constructor -- you will have to call setMatrix() before using the
  //! solve() function
  DualAMA() : Base() { init(); }
  //! Constructor with the system matrix
  explicit DualAMA(const BlockObjectBase<BlockMatrixType> &matrix) : Base() {
    init();
    Base::setMatrix(matrix);
  }

  //! Solve function using default variant
  template <typename NSLaw, typename MatrixT, typename RhsT, typename ResT>
  Scalar solve(const NSLaw &law, const BlockObjectBase<MatrixT> &A,
               const RhsT &f, const RhsT &w, ResT &v, ResT &r) const;

  //! Solve function with variant as template argument
  template <admm::Variant variant, typename NSLaw, typename MatrixT,
            typename RhsT, typename ResT>
  Scalar solve(const NSLaw &law, const BlockObjectBase<MatrixT> &A,
               const RhsT &f, const RhsT &w, ResT &v, ResT &r) const;

  //! Idem  with constraint  ( B v + k = 0 )
  template <admm::Variant variant, typename NSLaw, typename AType,
            typename BType, typename HType, typename PrecondT, typename RhsT,
            typename ORhsT, typename ResT, typename OResT>
  Scalar solveWithLinearConstraints(
      const NSLaw &law, const BlockObjectBase<AType> &A,
      const BlockObjectBase<BType> &B, const BlockObjectBase<HType> &H,
      const PrecondT &preconditioner, const RhsT &f, const ORhsT &k,
      const RhsT &w, ResT &v, OResT &p, ResT &r, Scalar stepRatio = 1) const;

  //! Sets the problem matrix -- the one defining the constraints
  DualAMA &setMatrix(const BlockObjectBase<BlockMatrixType> &matrix) {
    m_matrix = &matrix;
    Base::updateScalings();
    return *this;
  }

  //! Sets the step size for updating the dual variable (forces).
  void setFpStepSize(const Scalar size) { m_fpStepSize = size; }

  void setProjStepSize(const Scalar size) { m_projStepSize = size; }

  //! Sets the variant that will be used when calling solve() without template
  //! arguments
  void setDefaultVariant(admm::Variant variant) { m_defaultVariant = variant; }

  Scalar fpStepSize() const { return m_fpStepSize; }
  Scalar projStepSize() const { return m_projStepSize; }

  // LS

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

  unsigned lineSearchIterations() const { return m_lsIters; }
  Scalar lineSearchOptimisticFactor() const { return m_lsOptimisticFactor; }
  Scalar lineSearchPessimisticFactor() const { return m_lsPessimisticFactor; }

 protected:
  typedef typename Base::Index Index;

  //! Sets up the default values for all parameters
  void init() {
    m_tol = 1.e-6;
    m_maxIters = 300;

    m_fpStepSize = 1.e-1;
    m_projStepSize = 1.e-1;

    m_lsIters = 6;
    m_lsOptimisticFactor = 1.25;
    m_lsPessimisticFactor = .5;

    m_defaultVariant = admm::Accelerated;
  }

  using Base::m_matrix;
  using Base::m_maxIters;
  using Base::m_tol;

  Scalar m_fpStepSize;
  Scalar m_projStepSize;

  admm::Variant m_defaultVariant;

  unsigned m_lsIters;
  Scalar m_lsOptimisticFactor;
  Scalar m_lsPessimisticFactor;
};

}  // namespace bogus

#endif
