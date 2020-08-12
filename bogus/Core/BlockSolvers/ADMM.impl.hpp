/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_ADMM_IMPL_HPP
#define BOGUS_ADMM_IMPL_HPP

#include "../Block/Zero.hpp"
#include "../Utils/NumTraits.hpp"
#include "ADMM.hpp"
#include "ConstrainedSolverBase.impl.hpp"
#include "Preconditioners.impl.hpp"
#include "ProjectedGradient.impl.hpp"

namespace bogus {

//! Evaluation of prox_{1/c} (J) with J = 1/2 xM'x + f'x
/*! Requires way to solve linear system with lhs ( M + c I)
 */
template <typename ObjectType>
template <typename RhsT, typename ResT>
void QuadraticProxOp<ObjectType>::eval(RhsT &rhs, ResT &res) const {
  //  =  min J(x) + coeff/{2} || x  - rhs/coeff ||^2
  //  =  min J(x) + coeff/{2} || x ||^2 - <rhs, x>
  //  =  min 1/2 ( x'(M + coeff I)x  ) + <f - rhs, x>
  //  = (M + coeff I)^{-1} ( rhs - f )
  rhs -= m_affinePart;
  m_linearOp.template multiply<false>(rhs, res, 1, 0);
}

template <typename BlockMatrixType>
template <typename NSLaw, typename ProxOp, typename RhsT, typename ResT>
typename ADMM<BlockMatrixType>::Scalar ADMM<BlockMatrixType>::solve(
    const NSLaw &law, const ProxOp &op, const RhsT &w, ResT &v, ResT &r) const {
  switch (m_defaultVariant) {
    case admm::Standard:
      return solve<admm::Standard, NSLaw, ProxOp, RhsT, ResT>(law, op, w, v, r);
    case admm::Accelerated:
      return solve<admm::Accelerated, NSLaw, ProxOp, RhsT, ResT>(law, op, w, v,
                                                                 r);
  }

  return -1;
}

template <typename BlockMatrixType>
template <admm::Variant variant, typename NSLaw, typename ProxOp, typename RhsT,
          typename ResT>
typename ADMM<BlockMatrixType>::Scalar ADMM<BlockMatrixType>::solve(
    const NSLaw &law, const ProxOp &op, const RhsT &w, ResT &x, ResT &r) const {
  const Scalar lambda = op.coefficient();
  const Scalar gamma = stepSize();  // gamma/(lambda*lambda)
  const Scalar inv_gamma = 1. / gamma;

  Scalar res = -1;

  typename GlobalProblemTraits::DynVector ut, prox_arg(x.rows()), z;

  ut = w;
  m_matrix->template multiply<false>(x, ut, 1, 1);
  z = ut - inv_gamma * r;
  this->projectOnConstraints(law, z);

  // Nesterov acceleration
  typename GlobalProblemTraits::DynVector y, y_prev = -inv_gamma * r;
  typename GlobalProblemTraits::DynVector zacc = z, z_prev = z;
  Scalar theta_prev = 1.;  // Previous Nesterov acceleration
  Scalar res_prev = -1;

  for (unsigned adIter = 0; adIter < m_maxIters; ++adIter) {
    if (lambda < NumTraits<Scalar>::epsilon()) {
      // AMA
      m_matrix->template multiply<true>(r, prox_arg, 1, 0);
    } else {
      // ADMM
      prox_arg = x;
      ut = r +
           gamma *
               (zacc - ut);  // Re-use ut as we no-longer need its current value
      m_matrix->template multiply<true>(ut, prox_arg, 1, 1. / lambda);
    }
    op.eval(prox_arg, x);

    ut = w;
    m_matrix->template multiply<false>(x, ut, 1, 1);

    res =
        this->eval(law, r, ut) +
        (this->usesInfinityNorm() ? (ut - z).template lpNorm<Eigen::Infinity>()
                                  : (ut - z).squaredNorm() / (1 + ut.rows()));
    ;
    this->callback().trigger(adIter, res);

    if (res < this->tol()) break;

    z = ut - inv_gamma * r;
    this->projectOnConstraints(law, z);

    if (variant == admm::Accelerated && res < res_prev) {
      y = ut - z - inv_gamma * r;
      const Scalar beta = bogus::pg_impl::nesterov_inertia(theta_prev, 0.);

      r = -gamma * (y + beta * (y - y_prev));  // Over-relaxation
      y_prev = y;

      zacc = z + beta * (z - z_prev);
      z_prev = z;

    } else {
      r += gamma * (z - ut);

      zacc = z;
      z_prev = z;
      y_prev = -inv_gamma * r;
      theta_prev = 1;
    }

    //		const Scalar g = (ut - z).squaredNorm() ;
    //		std::cout << adIter << " \t gap: " << g << " \t res " << res <<
    //std::endl ;

    res_prev = res;
  }

  return res;
}

template <typename BlockMatrixType>
template <typename NSLaw, typename MatrixT, typename RhsT, typename ResT>
typename DualAMA<BlockMatrixType>::Scalar DualAMA<BlockMatrixType>::solve(
    const NSLaw &law, const BlockObjectBase<MatrixT> &A, const RhsT &f,
    const RhsT &w, ResT &v, ResT &r) const {
  switch (m_defaultVariant) {
    case admm::Accelerated:
      return solve<admm::Accelerated>(law, A, f, w, v, r);
    case admm::Standard:
      return solve<admm::Standard>(law, A, f, w, v, r);
  }

  return -1;
}

template <typename BlockMatrixType>
template <admm::Variant variant, typename NSLaw, typename MatrixT,
          typename RhsT, typename ResT>
typename DualAMA<BlockMatrixType>::Scalar DualAMA<BlockMatrixType>::solve(
    const NSLaw &law, const BlockObjectBase<MatrixT> &A, const RhsT &f,
    const RhsT &w, ResT &v, ResT &r) const {
  const typename LocalProblemTraits<0, Scalar>::Vector k;
  typename LocalProblemTraits<0, Scalar>::Vector p;
  const Zero<Scalar> zero;

  TrivialPreconditioner<BlockObjectBase<MatrixT> > precond;
  precond.setMatrix(A.derived());

  return solveWithLinearConstraints<variant>(law, A, zero, *m_matrix, precond,
                                             f, k, w, v, p, r);
}

template <typename BlockMatrixType>
template <admm::Variant variant, typename NSLaw, typename AType, typename BType,
          typename HType, typename PrecondT, typename RhsT, typename ORhsT,
          typename ResT, typename OResT>
typename DualAMA<BlockMatrixType>::Scalar
DualAMA<BlockMatrixType>::solveWithLinearConstraints(
    const NSLaw &law, const BlockObjectBase<AType> &A,
    const BlockObjectBase<BType> &B, const BlockObjectBase<HType> &H,
    const PrecondT &preconditioner, const RhsT &f, const ORhsT &k,
    const RhsT &w, ResT &v, OResT &p, ResT &r, Scalar stepRatio) const {
  typedef typename GlobalProblemTraits::DynVector DynVec;

  Scalar lambda = projStepSize();
  const Scalar gamma = fpStepSize();

  Scalar res = -1, min_res = -1;

  typename GlobalProblemTraits::DynVector r_best(r), v_best(v), p_best(p);
  typename GlobalProblemTraits::DynVector z(v.rows()), x, ut, gap,
      HrBp(v.rows()), g, g2, s(r.rows());

  const Segmenter<NSLaw::dimension, const DynVec,
                  typename BlockMatrixType::Index>
      utSegmenter(ut, m_matrix->rowOffsets());
  Segmenter<NSLaw::dimension, DynVec, typename BlockMatrixType::Index>
      sSegmenter(s, m_matrix->rowOffsets());

  // Nesterov acceleration
  DynVec y, y_prev = v;
  Scalar theta_prev = 1.;  // Previous Nesterov acceleration
  Scalar res_prev = -1;

  // Line search
  DynVec rs, ps;

  H.template multiply<true>(r, HrBp, 1, 0);
  B.template multiply<true>(p, HrBp, 1, 1);

  for (unsigned adIter = 0; adIter < m_maxIters; ++adIter) {
    x = f;
    A.template multiply<false>(v, x, 1, 1);
    gap = HrBp - x;
    preconditioner.template apply<false>(gap, z);

    ut = w;
    H.template multiply<false>(v, ut, 1, 1);

    g2 = k;
    B.template multiply<false>(v + gamma * z, g2, 1, 1);

    // Eval current reisual,  exit if small enough
    res = this->eval(law, ut, r);     // Complementarity
    res += (this->usesInfinityNorm()  // Gap
                ? z.template lpNorm<Eigen::Infinity>()
                : z.squaredNorm() / (1 + z.rows()));
    if (g2.rows() > 0)
      res += (this->usesInfinityNorm()  // Linear constraints
                  ? g2.template lpNorm<Eigen::Infinity>()
                  : g2.squaredNorm() / (1 + g2.rows()));

    this->callback().trigger(adIter, res);

    if (res < min_res || adIter == 0) {
      r_best = r;
      v_best = v;
      p_best = p;
      min_res = res;
      if (res < this->tol()) break;
    } else if (res > 1.e50)
      break;

    // Acceleration
    Scalar beta = 0;
    if (variant == admm::Accelerated && res < res_prev) {
      beta = bogus::pg_impl::nesterov_inertia(theta_prev, 0.);
    } else {
      theta_prev = 1;
    }

    this->dualityCOV(law, ut, s);
    g = ut + s;
    H.template multiply<false>(z, g, gamma, 1);

    // Line search (optional)
    if (lineSearchIterations() == 0) {
      r -= lambda * g;
      this->projectOnConstraints(law, r);

      p -= stepRatio * lambda * g2;

      H.template multiply<true>(r, HrBp, 1, 0);
      B.template multiply<true>(p, HrBp, 1, 1);

    } else {
      // TODO avoid performing extra A multiplication by keeping LS results

      const Scalar h0 = gap.squaredNorm();

      lambda *= lineSearchOptimisticFactor();
      for (unsigned lsIter = 0; lsIter < lineSearchIterations(); ++lsIter) {
        rs = r - lambda * g;
        this->projectOnConstraints(law, rs);

        ps = p - stepRatio * lambda * g2;

        H.template multiply<true>(rs, HrBp, 1, 0);
        B.template multiply<true>(ps, HrBp, 1, 1);

        gap = HrBp - x;
        preconditioner.template apply<false>(gap, z);

        y = v + gamma * z;
        y += beta * (y - y_prev);  // Over Relaxation

        gap = HrBp - f;
        A.template multiply<false>(y, gap, -1, 1);

        Scalar h = gap.squaredNorm();

        if (h < h0) {
          //					std::cout << lsIter << " [ "  << lambda<< " ] "
          //<< " \t " << h << " vs " << h0 << std::endl ;
          break;
        }

        lambda *= lineSearchPessimisticFactor();
      }

      r = rs;
      p = ps;
    }
    // (end line-search)

    gap = HrBp - x;
    preconditioner.template apply<false>(gap, z);

    y = v + gamma * z;
    v = y + beta * (y - y_prev);  // Over Relaxation

    y_prev = y;
    res_prev = res;
  }

  r = r_best;
  v = v_best;
  p = p_best;

  return min_res;
}

}  // namespace bogus

#endif
