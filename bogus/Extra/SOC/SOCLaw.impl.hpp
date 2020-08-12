/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * So-bogus is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * So-bogus is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with So-bogus.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BOGUS_SOCLAW_IMPL_HPP
#define BOGUS_SOCLAW_IMPL_HPP

#include "../../Core/Utils/NumTraits.hpp"
#include "FischerBurmeister.hpp"
#include "LocalSOCSolver.hpp"
#include "LocalSOCSolver.impl.hpp"
#include "SOCLaw.hpp"

namespace bogus {

template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV,
          local_soc_solver::Strategy Strat>
SOCLaw<Dimension, Scalar, DeSaxceCOV, Strat>::SOCLaw(const unsigned n,
                                                     const double *mu)
    : m_mu(mu),
      m_n(n),
      m_localTol(std::pow(NumTraits<Scalar>::epsilon(), .75)) {}

template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV,
          local_soc_solver::Strategy Strat>
bool SOCLaw<Dimension, Scalar, DeSaxceCOV, Strat>::solveLocal(
    const unsigned problemIndex, const typename Traits::Matrix &A,
    const typename Traits::Vector &b, typename Traits::Vector &xm,
    const Scalar scaling) const {
  typedef LocalSOCSolver<Traits::dimension, typename Traits::Scalar, DeSaxceCOV,
                         Strat>
      LocalSolver;
  return m_localTol >
         LocalSolver::solve(A, b, xm, m_mu[problemIndex], m_localTol, scaling);
}

template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV,
          local_soc_solver::Strategy Strat>
void SOCLaw<Dimension, Scalar, DeSaxceCOV, Strat>::projectOnConstraint(
    const unsigned problemIndex, typename Traits::Vector &x) const {
  const Scalar nxt = Traits::tp(x).norm();
  const Scalar mu = m_mu[problemIndex];
  const Scalar xn = Traits::np(x);

  if (mu < 0) {
    x.setZero();
    return;
  }

  // x inside cone
  if (nxt <= mu * xn) return;

  // x inside dual cone
  if (mu * nxt <= -xn) {
    x.setZero();
    return;
  }

  const Scalar den2 = (1 + mu * mu);
  const Scalar proj = (xn + mu * nxt) / den2;

  Traits::np(x) = proj;
  Traits::tp(x) *= proj * mu / nxt;
}

}  // namespace bogus

#endif
