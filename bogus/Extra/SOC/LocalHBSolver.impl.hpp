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

#ifndef BOGUS_LOCAL_HB_SOLVER_IMPL_HPP
#define BOGUS_LOCAL_HB_SOLVER_IMPL_HPP

#include "../../Core/Utils/LinearSolverBase.hpp"
#include "../../Core/Utils/NonSmoothNewton.impl.hpp"
#include "../../Core/Utils/Polynomial.impl.hpp"
#include "FischerBurmeister.impl.hpp"
#include "LocalHBSolver.hpp"

namespace bogus {

template <DenseIndexType Dimension, typename Scalar,
          local_soc_solver::Strategy Strat>
Scalar LocalHBSolver<Dimension, Scalar, Strat>::solve(
    const typename Traits::Matrix &A, const typename Traits::Vector &b,
    typename Traits::Vector &x, const Scalar yield, const Scalar eta,
    const Scalar n, const Scalar tol, const Scalar scaling) {
  typedef HerschelBulkleyFischerBurmeister<Dimension, Scalar> FBFunc;
  typedef typename Traits::LUType LUType;

  FBFunc fb(yield, eta, n, A, b, scaling);
  NonSmoothNewton<FBFunc> nsNewton(fb, tol);

  if (Strat == local_soc_solver::PureNewton) {
    return nsNewton.solve(x);
  }

  double res = 0.;

  if (Strat == local_soc_solver::Hybrid) {
    res = nsNewton.solve(x);
    if (res < tol) return res;
  }

  Vector x0 = x;
  LUType(A).solve(-b, x);

  double nx = Traits::tp(x).norm();
  if (nx < 1e-20) {
    // Degenerated case, ignore
    return 0.;
  }

  double yield_func = fb.yieldFunction(x);
  if (yield_func >= 0.) {
    // Sticking case
    return 0.;
  }

  // Sliding case
  double yield_criteria = yield_func + nx;
  Traits::tp(x) *= yield_criteria / nx;

  // Refinement of final solution
  if (Strat == local_soc_solver::RevHybrid) {
    res = nsNewton.solve(x);
  } else if (Strat == local_soc_solver::Hybrid) {
    const double refinedRes = nsNewton.solve(x);
    if (refinedRes <= res) return refinedRes;

    // This can happen if the solver returned a very bad value, like an
    // unreastically huge alpha
    x = x0;
  }

  return res;
}

}  // namespace bogus

#endif
