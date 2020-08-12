/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_PYRAMIDLAW_IMPL_HPP
#define BOGUS_PYRAMIDLAW_IMPL_HPP

#include "../BlockSolvers.hpp"
#include "../Utils/CppTools.hpp"
#include "../Utils/NumTraits.hpp"
#include "PyramidLaw.hpp"

namespace bogus {

template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV>
PyramidLaw<Dimension, Scalar, DeSaxceCOV>::PyramidLaw(const unsigned n,
                                                      const double *mu)
    : m_mu(mu), m_n(n) {}

template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV>
void PyramidLaw<Dimension, Scalar, DeSaxceCOV>::projectOnConstraint(
    const unsigned problemIndex, typename Traits::Vector &x) const {
  BOGUS_STATIC_ASSERT(Dimension < 4u, NOT_IMPLEMENTED);

  typedef typename LocalProblemTraits<Dimension - 1, Scalar>::Array TgComp;
  const Scalar mu = m_mu[problemIndex];

  const Scalar mun = mu * Traits::np(x);

  TgComp diff = Traits::tp(x).array().abs() - mun;
  int kMax;
  if (diff.maxCoeff(&kMax) > 0) {
    TgComp signs;
    for (unsigned k = 0; k < Dimension - 1; ++k)
      signs[k] = x[1 + k] > 0 ? 1 : -1;

    bool ok = false;

    // Check for projection on a face
    TgComp depls = signs * (mun - Traits::tp(x).array().abs());
    TgComp ss = mu * depls;

    typename Traits::Vector xproj = x;
    Traits::np(xproj) += ss[kMax] / (1 + mu * mu);
    xproj[1 + kMax] += depls[kMax] / (1 + mu * mu);

    if (mu * Traits::np(xproj) >=
        Traits::tp(xproj).template lpNorm<Eigen::Infinity>()) {
      x = xproj;
      ok = true;
    }

    // Projection on an edge
    if (!ok) {
      typename Traits::Vector n;
      Traits::np(n) = 1;
      Traits::tp(n) = mu * signs;

      x = x.dot(n) * n / n.squaredNorm();
    }

    // Inside normal cone
    if (Traits::np(x) < 0) {
      x.setZero();
    }
  }
}

template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV>
bool PyramidLaw<Dimension, Scalar, DeSaxceCOV>::solveLocal(
    const unsigned problemIndex, const typename Traits::Matrix &A,
    const typename Traits::Vector &b, typename Traits::Vector &x,
    const Scalar) const {
  BOGUS_STATIC_ASSERT(DeSaxceCOV, NOT_IMPLEMENTED);

  // Solve for
  // 0 \leq x_N \perp y_N \geq 0
  // 0 \leq \mu x_N - x^-_T \perp y^+_N \geq 0
  // 0 \leq \mu x_N - x^+_T \perp y^-_N \geq 0
  // x = x^+ - x^-

  const Scalar mu = m_mu[problemIndex];

  typename Traits::Vector xp = x.cwiseMax(0).array();
  typename Traits::Vector xm = xp - x;

  // Solve dxd problem using GS
  // TODO solve analytically or with fb, deal with not positive diag coeff
  for (unsigned it = 0; it < 2 * Dimension; ++it) {
    // Normal comp
    xm[0] = 0;
    Scalar l = A.col(0).dot(xp - xm) - A(0, 0) * xp[0] + b[0];
    xp[0] = std::max(0., -l / A(0, 0));

    Scalar mun = mu * xp[0];

    // Tangential comps
    for (unsigned k = 1; k < Dimension; ++k) {
      Scalar l = A.col(k).dot(xp - xm) - A(k, k) * (xp[k] - xm[k]) + b[k];
      xp[k] = std::min(mun, -std::min(0., l) / A(k, k));
      xm[k] = std::min(mun, std::max(0., l) / A(k, k));
    }
  }

  // Recompose x as x^+ - x^-
  x = xp - xm;

  return true;
}

template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV>
Scalar PyramidLaw<Dimension, Scalar, DeSaxceCOV>::eval(
    const unsigned problemIndex, const typename Traits::Vector &x,
    const typename Traits::Vector &ycov) const {
  typedef typename LocalProblemTraits<Dimension - 1, Scalar>::Array TgComp;
  const Scalar mu = m_mu[problemIndex];

  typename Traits::Vector y = ycov;
  if (!DeSaxceCOV) {
    // Remove controbution from change of variable
    dualityCOV(problemIndex, ycov, y);
    y = ycov - y;
  }

  // Evaluate error as summ of Fischer-Burmeister residual on each comp

  TgComp xp = Traits::tp(x).array().max(TgComp::Zero(Dimension - 1));
  TgComp yp = Traits::tp(y).array().max(TgComp::Zero(Dimension - 1));
  TgComp xm = xp - Traits::tp(x).array();
  TgComp ym = yp - Traits::tp(y).array();

  const Scalar xn = Traits::np(x);
  const Scalar yn = Traits::np(y);

  const Scalar fb = xn + yn - std::sqrt(xn * xn + yn * yn);
  Scalar err = fb * fb;

  xm = mu * Traits::np(x) - xm;
  xp = mu * Traits::np(x) - xp;

  TgComp fbm = xm + yp - (xm * xm + yp * yp).sqrt();
  TgComp fbp = xp + ym - (xp * xp + ym * ym).sqrt();

  err += fbm.matrix().squaredNorm() + fbp.matrix().squaredNorm();
  return err;
}

}  // namespace bogus

#endif
