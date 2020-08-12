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

#ifndef BOGUS_HBLAW_IMPL_HPP
#define BOGUS_HBLAW_IMPL_HPP

#include "../../Core/Utils/NumTraits.hpp"
#include "FischerBurmeister.hpp"
#include "HBLaw.hpp"
#include "LocalHBSolver.hpp"
#include "LocalHBSolver.impl.hpp"

namespace bogus {

template <DenseIndexType Dimension, typename Scalar,
          local_soc_solver::Strategy Strat>
HBLaw<Dimension, Scalar, Strat>::HBLaw(const unsigned n, const double *yields,
                                       const double *etas, const double *powers)
    : m_yields(yields),
      m_etas(etas),
      m_powers(powers),
      m_n(n),
      m_localTol(std::pow(NumTraits<Scalar>::epsilon(), .75)) {}

template <DenseIndexType Dimension, typename Scalar,
          local_soc_solver::Strategy Strat>
bool HBLaw<Dimension, Scalar, Strat>::solveLocal(
    const unsigned problemIndex, const typename Traits::Matrix &A,
    const typename Traits::Vector &b, typename Traits::Vector &xm,
    const Scalar scaling) const {
  typedef LocalHBSolver<Traits::dimension, typename Traits::Scalar, Strat>
      LocalSolver;
  return m_localTol > LocalSolver::solve(A, b, xm, m_yields[problemIndex],
                                         m_etas[problemIndex],
                                         m_powers[problemIndex], m_localTol,
                                         scaling);
}

template <DenseIndexType Dimension, typename Scalar,
          local_soc_solver::Strategy Strat>
void HBLaw<Dimension, Scalar, Strat>::projectOnConstraint(
    const unsigned problemIndex, typename Traits::Vector &x) const {
  double nx = Traits::tp(x).norm();
  if (nx < 1e-20) {
    // Degenerated case, ignore
    return;
  }

  typedef HerschelBulkleyFischerBurmeister<Dimension, Scalar> FBFunc;
  const Scalar yield_func =
      FBFunc::yieldFunction(m_yields[problemIndex], m_etas[problemIndex],
                            m_powers[problemIndex], x, Traits::Vector::Zero());

  if (yield_func >= 0.) return;

  double yield_criteria = yield_func + nx;
  Traits::tp(x) *= yield_criteria / nx;
}
}  // namespace bogus

#endif
