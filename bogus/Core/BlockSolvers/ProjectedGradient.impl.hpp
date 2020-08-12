/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_PROJECTED_GRADIENT_IMPL_HPP
#define BOGUS_PROJECTED_GRADIENT_IMPL_HPP

#include "ConstrainedSolverBase.impl.hpp"
#include "ProjectedGradient.hpp"

namespace bogus {

namespace pg_impl {

template <typename Scalar, typename Mat, typename Vec1, typename Vec2>
void guess_alpha(const Mat &M, const Vec1 &x, Scalar &alpha, Vec2 &tmp) {
  // Eval optimal alpha (step size )
  tmp.setOnes();
  const Scalar nMx = (M * (x - tmp)).squaredNorm();
  if (nMx > 1.e-16)  // Don't try too big alphas
    alpha = std::min(1., (x - tmp).squaredNorm() / nMx);
}

template <typename Scalar>
Scalar nesterov_inertia(Scalar &theta, const Scalar q) {
  const Scalar theta_prev = theta;
  const Scalar theta_p2 = theta_prev * theta_prev;
  const Scalar bb = theta_p2 - q;
  const Scalar delta = std::sqrt(bb * bb + 4 * theta_p2);
  theta = .5 * (delta - bb);

  return (theta_prev - theta_p2) / (theta_p2 + theta);
}

template <typename BlockMatrixType, typename NSLaw, typename VecX,
          typename VecY, typename VecRes, typename Scalar>
bool test_residual(const ProjectedGradient<BlockMatrixType> &pg,
                   const NSLaw &law, const unsigned pgIter, const VecX &x,
                   const VecY &y, VecRes &x_best, Scalar &min_res) {
  const Scalar res = pg.eval(law, y, x);

  pg.callback().trigger(pgIter, res);
  if (0 == pgIter || res < min_res) {
    x_best = x;
    min_res = res;
  }
  return (res < pg.tol());
}

template <typename BlockMatrixType, typename NSLaw, typename MatrixT,
          typename VecX, typename VecB, typename VecBuf, typename Scalar>
void armijo_ls(const ProjectedGradient<BlockMatrixType> &pg, const NSLaw &law,
               const MatrixT &M, const VecX &x, const VecB &b,
               const VecBuf &dir, VecBuf &grad, const Scalar J, Scalar &alpha,
               VecBuf &xs, Scalar &Js, VecBuf &Mx) {
  Js = J;
  for (unsigned lsIter = 0; lsIter < pg.lineSearchIterations();
       ++lsIter, alpha *= pg.lineSearchPessimisticFactor()) {
    xs = x + alpha * dir;

    pg.projectOnConstraints(law, xs);

    Mx = M * xs;
    Js = xs.dot(.5 * Mx + b);

    const Scalar decr = grad.dot(xs - x);

    if (Js < J + pg.lineSearchArmijoCoefficient() * decr) break;
  }
}

template <projected_gradient::Variant variant>
struct PgMethod {
  // variant = Descent, APGD

  template <typename BlockMatrixType, typename NSLaw, typename MatrixT,
            typename RhsT, typename ResT>
  static typename ProjectedGradient<BlockMatrixType>::Scalar solve(
      const ProjectedGradient<BlockMatrixType> &pg, const NSLaw &law,
      const MatrixT &M, const RhsT &b, ResT &x) {
    typedef ProjectedGradient<BlockMatrixType> PgType;
    typedef typename PgType::Scalar Scalar;
    typename PgType::GlobalProblemTraits::DynVector Mx(b.rows()),
        y(b.rows()),   // = Mx +b  (gradient)
        xs(x.rows()),  // tentative new value for x
        s(b.rows()),   // Duality COV
        prev_proj, x_best;

    pg.projectOnConstraints(law, x);

    // Unconstrained objective function
    Mx = M * x;

    // Nesterov inertia
    Scalar theta_prev = 1., q = 0;

    prev_proj = x;

    Scalar alpha = 1.;
    guess_alpha(M, x, alpha, xs);

    // Best iterate storage
    Scalar min_res = -1;

    for (unsigned pgIter = 0; pgIter < pg.maxIters(); ++pgIter) {
      // y = grad J = Mx+b
      // The residual should be evaluated in prev_proj instead of x
      // However this would been one more matrix product per iterations
      y = Mx + b;
      if (test_residual(pg, law, pgIter, x, y, x_best, min_res)) break;

      pg.dualityCOV(law, y, s);
      y += s;
      const Scalar J = x.dot(.5 * Mx + b + s);

      Scalar Js = J;
      alpha *= pg.lineSearchOptimisticFactor();

      // Line-search
      for (unsigned lsIter = 0; lsIter < pg.lineSearchIterations();
           ++lsIter, alpha *= pg.lineSearchPessimisticFactor()) {
        xs = x - alpha * y;

        pg.projectOnConstraints(law, xs);

        Mx = M * xs;
        Js = xs.dot(.5 * Mx + b + s);

        const Scalar decr = y.dot(xs - x);

        if (Js < J + decr + .5 / alpha * (xs - x).squaredNorm()) break;
      }

      const Scalar decr = (xs - prev_proj).dot(y);

      if (variant == projected_gradient::APGD && decr < 0) {
        const Scalar beta = nesterov_inertia(theta_prev, q);
        x = xs + beta * (xs - prev_proj);

        Mx = M * x;

      } else {
        x = xs;

        theta_prev = 1.;  // APGD
      }

      prev_proj = xs;
    }

    x = x_best;
    return min_res;
  }
};

template <>
struct PgMethod<projected_gradient::Standard> {
  template <typename BlockMatrixType, typename NSLaw, typename MatrixT,
            typename RhsT, typename ResT>
  static typename ProjectedGradient<BlockMatrixType>::Scalar solve(
      const ProjectedGradient<BlockMatrixType> &pg, const NSLaw &law,
      const MatrixT &M, const RhsT &b, ResT &x) {
    typedef ProjectedGradient<BlockMatrixType> PgType;
    typedef typename PgType::Scalar Scalar;
    typename PgType::GlobalProblemTraits::DynVector Mx(b.rows()),
        y(b.rows()),   // = Mx +b  (gradient)
        xs(x.rows()),  // tentative new value for x
        s(b.rows()),   // Duality COV
        proj_grad(b.rows()), x_best;

    pg.projectOnConstraints(law, x);

    // Unconstrained objective function
    Mx = M * x;

    Scalar alpha = 1;

    // Best iterate storage
    Scalar min_res = -1;

    for (unsigned pgIter = 0; pgIter < pg.maxIters(); ++pgIter) {
      y = Mx + b;
      if (test_residual(pg, law, pgIter, x, y, x_best, min_res)) break;

      pg.dualityCOV(law, y, s);
      y += s;
      const Scalar J = x.dot(.5 * Mx + b + s);

      xs = x - y;
      pg.projectOnConstraints(law, xs);

      // proj_grad ~ projection of gradient on tangent cone
      proj_grad = (x - xs);
      const Scalar beta = -y.dot(proj_grad) / (proj_grad.dot(M * proj_grad));

      proj_grad *= beta;

      // Line-search
      Scalar Js;
      alpha = std::min(1., alpha * pg.lineSearchOptimisticFactor());
      armijo_ls(pg, law, M, x, b + s, proj_grad, y, J, alpha, xs, Js, Mx);

      x = xs;
    }

    x = x_best;
    return min_res;
  }
};

template <>
struct PgMethod<projected_gradient::Conjugated> {
  template <typename BlockMatrixType, typename NSLaw, typename MatrixT,
            typename RhsT, typename ResT>
  static typename ProjectedGradient<BlockMatrixType>::Scalar solve(
      const ProjectedGradient<BlockMatrixType> &pg, const NSLaw &law,
      const MatrixT &M, const RhsT &b, ResT &x) {
    typedef ProjectedGradient<BlockMatrixType> PgType;
    typedef typename PgType::Scalar Scalar;
    typename PgType::GlobalProblemTraits::DynVector Mx(b.rows()),
        y(b.rows()),   // = Mx +b  (gradient)
        xs(x.rows()),  // tentative new value for x
        s(b.rows()),   // Duality COV
        dir, prev_dir, prev_proj_grad, x_best;

    pg.projectOnConstraints(law, x);

    // Unconstrained objective function
    Mx = M * x;

    prev_proj_grad = x;

    Scalar prev_n2 = 0;
    Scalar alpha = 1.;

    // Best iterate storage
    Scalar min_res = -1;

    for (unsigned pgIter = 0; pgIter < pg.maxIters(); ++pgIter) {
      // y = grad J = Mx+b
      y = Mx + b;
      if (test_residual(pg, law, pgIter, x, y, x_best, min_res)) break;

      pg.dualityCOV(law, y, s);
      y += s;
      const Scalar J = x.dot(.5 * Mx + b + s);

      {
        // xs ~ projection of gradient on tangent cone
        xs = x - y;
        pg.projectOnConstraints(law, xs);
        xs = x - xs;

        // Conjugation
        const Scalar ng2 = xs.squaredNorm();
        if (prev_n2 == 0.) {
          dir = -y;
        } else {
          const Scalar den = prev_dir.dot(xs - prev_proj_grad);
          const Scalar beta =
              (den < 1.e-12)
                  // Polak-Ribiere
                  ? std::max(0., (ng2 - xs.dot(prev_proj_grad))) / prev_n2
                  // Hestness-Stiefel
                  : std::max(0., (ng2 - xs.dot(prev_proj_grad))) / den;

          dir = beta * prev_dir - y;

          if (dir.dot(y) >= 0.) dir = -y;
        }

        prev_proj_grad = xs;
        prev_n2 = ng2;
      }

      Scalar Js;
      alpha = alpha * pg.lineSearchOptimisticFactor();
      armijo_ls(pg, law, M, x, b + s, dir, y, J, alpha, xs, Js, Mx);

      prev_dir = (xs - x) / alpha;

      x = xs;
    }

    x = x_best;
    return min_res;
  }
};

template <>
struct PgMethod<projected_gradient::SPG> {
  template <typename BlockMatrixType, typename NSLaw, typename MatrixT,
            typename RhsT, typename ResT>
  static typename ProjectedGradient<BlockMatrixType>::Scalar solve(
      const ProjectedGradient<BlockMatrixType> &pg, const NSLaw &law,
      const MatrixT &M, const RhsT &b, ResT &x) {
    typedef ProjectedGradient<BlockMatrixType> PgType;
    typedef typename PgType::Scalar Scalar;
    typename PgType::GlobalProblemTraits::DynVector Mx(b.rows()),
        y(b.rows()),   // = Mx +b  (gradient)
        xs(x.rows()),  // tentative new value for x
        s(b.rows()),   // Duality COV
        z(x.rows()), w(x.rows()), prev_proj, x_best;

    pg.projectOnConstraints(law, x);
    prev_proj = x;

    Mx = M * x;

    Scalar theta_prev = 1, q = 0;
    Scalar a_min = 1.e-6, a_max = 1.e6;

    Scalar alpha = 1;
    guess_alpha(M, x, alpha, xs);

    // Best iterate storage
    Scalar min_res = -1;

    for (unsigned pgIter = 0; pgIter < pg.maxIters(); ++pgIter) {
      y = Mx + b;
      if (test_residual(pg, law, pgIter, x, y, x_best, min_res)) break;

      pg.dualityCOV(law, y, s);
      y += s;

      xs = x - alpha * y;
      pg.projectOnConstraints(law, xs);

      if ((xs - prev_proj).dot(y) < 0.) {
        const Scalar beta = nesterov_inertia(theta_prev, q);

        w = (xs + beta * (xs - prev_proj)) - x;
        x += w;

      } else {
        w = xs - x;
        x = xs;

        theta_prev = 1.;
      }

      z = -Mx;
      Mx = M * x;
      prev_proj = xs;

      z += Mx;
      alpha = w.dot(z) / z.squaredNorm();
      alpha = std::min(a_max, std::max(a_min, alpha));
    }

    x = x_best;
    return min_res;
  }
};

}  // namespace pg_impl

template <typename BlockMatrixType>
template <typename NSLaw, typename RhsT, typename ResT>
typename ProjectedGradient<BlockMatrixType>::Scalar
ProjectedGradient<BlockMatrixType>::solve(const NSLaw &law, const RhsT &b,
                                          ResT &x) const {
  switch (m_defaultVariant) {
    case projected_gradient::Standard:
      return solve<projected_gradient::Standard, NSLaw, RhsT, ResT>(law, b, x);
    case projected_gradient::Descent:
      return solve<projected_gradient::Descent, NSLaw, RhsT, ResT>(law, b, x);
    case projected_gradient::Conjugated:
      return solve<projected_gradient::Conjugated, NSLaw, RhsT, ResT>(law, b,
                                                                      x);
    case projected_gradient::APGD:
      return solve<projected_gradient::APGD, NSLaw, RhsT, ResT>(law, b, x);
    case projected_gradient::SPG:
      return solve<projected_gradient::SPG, NSLaw, RhsT, ResT>(law, b, x);
  }

  return -1;
}

template <typename BlockMatrixType>
template <projected_gradient::Variant variant, typename NSLaw, typename RhsT,
          typename ResT>
typename ProjectedGradient<BlockMatrixType>::Scalar
ProjectedGradient<BlockMatrixType>::solve(const NSLaw &law, const RhsT &b,
                                          ResT &x) const {
  return pg_impl::PgMethod<variant>::solve(*this, law, Base::matrix(), b, x);
}

}  // namespace bogus

#endif
