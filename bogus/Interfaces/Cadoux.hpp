/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2015 Gilles Daviet <gdaviet@gmail.com>
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

#ifndef BOGUS_CADOUX_HPP
#define BOGUS_CADOUX_HPP

#include "../Core/BlockSolvers/ConstrainedSolverBase.impl.hpp"
#include "../Core/Utils/CppTools.hpp"
#include "../Extra/SecondOrder.impl.hpp"

namespace bogus {

//! Solves a Coulomb friction problem using the Cadoux algorithm ( with
//! fixed-point iteration )
/*!
          See \cite ACLM11
          \param W  The Delassus operator
          \param b  The affine term of the linear system u = Wr + b
          \param mu The vector of friction coefficients
          \param minimizer The minimizer to use for each inner SOCQP problem
          \param r  Both the initial guess and the result
          \param cadouxIterations Number of fixed-point iterations
          \param callback NULL, or a pointer to a user-defined function that
   takes ( unsigned iteration, double residual ) as arguments \param tolTighten
   How much should the tolerance of the inner solver be tightened \returns the
   error as returned by the minimizer eval() function
          */
template <unsigned Dimension, typename WType, typename Method, typename MatrixT>
static typename WType::Scalar solveCadoux(
    const WType& W, const typename WType::Scalar* b,
    const typename WType::Scalar* mu,
    ConstrainedSolverBase<Method, MatrixT>& minimizer,
    typename WType::Scalar* r, const unsigned cadouxIterations,
    const Signal<unsigned, typename WType::Scalar>* callback =
        BOGUS_NULL_PTR(const void),
    const typename WType::Scalar tolTighten = 1.e-1) {
  // We might experience slow convergence is inner solve not precise enough

  typedef typename WType::Scalar Scalar;
  const std::ptrdiff_t n = W.rowsOfBlocks();

  SOCLaw<Dimension, Scalar, true> coulombLaw(n, mu);
  SOCLaw<Dimension, Scalar, false> socLaw(n, mu);

  Eigen::Map<Eigen::VectorXd> r_map(r, W.rows());
  Eigen::Map<const Eigen::VectorXd> b_map(b, W.rows());

  Eigen::VectorXd s(W.rows());

  Scalar res = -1;
  const Scalar tol = minimizer.tol();

  // Evaluate initial error
  s = W * r_map + b_map;
  res = minimizer.eval(coulombLaw, s, r_map);
  if (callback) callback->trigger(0, res);

  for (unsigned cdxIter = 0; cdxIter < cadouxIterations; ++cdxIter) {
    minimizer.setTol(tolTighten * std::max(tol, std::min(res, 1.)));

    minimizer.dualityCOV(coulombLaw, s, s);

    s += b_map;
    minimizer.solve(socLaw, s, r_map);

    // Evaluate current error
    s = W * r_map + b_map;
    res = minimizer.eval(coulombLaw, s, r_map);

    if (callback) callback->trigger(cdxIter + 1, res);
    if (res < tol) break;
  }

  minimizer.setTol(tol);

  return res;
}

//! Same as solveCadoux, with r = W*u +b
/*! \warning requires mu > 0 */
template <unsigned Dimension, typename WType, typename Method, typename MatrixT>
static double solveCadouxVel(const WType& W, const typename WType::Scalar* b,
                             const typename WType::Scalar* mu,
                             ConstrainedSolverBase<Method, MatrixT>& minimizer,
                             typename WType::Scalar* u,
                             const unsigned cadouxIterations,
                             const Signal<unsigned, typename WType::Scalar>*
                                 callback = BOGUS_NULL_PTR(const void),
                             const typename WType::Scalar tolTighten = 1.e-1) {
  // Wu + b = r
  // u* = u + s n
  // Wu* + b - W(s n) = r

  typedef typename WType::Scalar Scalar;
  const std::ptrdiff_t n = W.rowsOfBlocks();

  SOCLaw<Dimension, Scalar, false> socLaw(n, mu);

  Eigen::Map<Eigen::VectorXd> u_map(u, W.rows());
  Eigen::Map<const Eigen::VectorXd> b_map(b, W.rows());

  Eigen::VectorXd s(W.rows()), Wsb(W.rows()), ustar(u_map), r(W.rows());

  double res = -1;
  const double tol = minimizer.tol();

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
  for (std::ptrdiff_t i = 0; i < n; ++i) {
    s[Dimension * i] = u_map.segment<Dimension - 1>(Dimension * i + 1).norm() *
                       std::max(0., 1. / mu[i]);
    s.segment<Dimension - 1>(Dimension * i + 1).setZero();
  }
  ustar = u_map + s;

  // Evaluate intial error
  r = W * u_map + b_map;
  res = minimizer.eval(socLaw, r, ustar);
  if (callback) callback->trigger(0, res);

  for (unsigned cdxIter = 0; cdxIter < cadouxIterations; ++cdxIter) {
    minimizer.setTol(tolTighten * std::max(tol, std::min(res, 1.)));

    Wsb = b_map - W * s;

    minimizer.solve(socLaw, Wsb, ustar);

    u_map = ustar - s;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
    for (std::ptrdiff_t i = 0; i < n; ++i) {
      s[Dimension * i] =
          ustar.segment<Dimension - 1>(Dimension * i + 1).norm() *
          std::max(0., 1. / mu[i]);
      s.segment<Dimension - 1>(Dimension * i + 1).setZero();
    }

    ustar = u_map + s;

    // Evaluate current error
    r = W * u_map + b_map;
    res = minimizer.eval(socLaw, r, ustar);

    if (callback) callback->trigger(cdxIter + 1, res);
    if (cdxIter > 0 && res < tol) break;
  }

  minimizer.setTol(tol);

  return res;
}

}  // namespace bogus

#endif
