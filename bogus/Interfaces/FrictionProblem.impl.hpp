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

#include "../Core/Block.impl.hpp"
#include "../Extra/SecondOrder.impl.hpp"
#include "Cadoux.hpp"
#include "FrictionProblem.hpp"

namespace bogus {

namespace friction_problem {

template <unsigned Dimension, typename EigenDerived, typename Index>
void applyPermutation(const std::vector<std::size_t> &permutation,
                      Eigen::MatrixBase<EigenDerived> &vec,
                      const Index *offsets) {
  Segmenter<Dimension, EigenDerived, Index> segmenter(vec.derived(), offsets);
  bogus::applyPermutation(permutation.size(), data_pointer(permutation),
                          segmenter);
}

template <unsigned Dimension, typename Method, typename MatrixT>
static double solve(const DualFrictionProblem<Dimension> &dual,
                    const ConstrainedSolverBase<Method, MatrixT> &gs, double *r,
                    const bool staticProblem,
                    const bool herschelBulkleyProblem) {
  typename Eigen::VectorXd::MapType r_map(r, dual.W.rows());

  if (dual.permuted())
    applyPermutation<Dimension>(dual.permutation(), r_map,
                                dual.W.majorIndex().innerOffsetsData());

  double res;

  if (herschelBulkleyProblem) {
    res = gs.solve(typename DualFrictionProblem<Dimension>::HBLawType(
                       dual.W.rowsOfBlocks(), dual.yields.data(),
                       dual.etas.data(), dual.powers.data()),
                   dual.b, r_map);
  } else {
    res =
        staticProblem
            ? gs.solve(typename DualFrictionProblem<Dimension>::SOCLawType(
                           dual.W.rowsOfBlocks(), dual.mu.data()),
                       dual.b, r_map)
            : gs.solve(typename DualFrictionProblem<Dimension>::CoulombLawType(
                           dual.W.rowsOfBlocks(), dual.mu.data()),
                       dual.b, r_map);
  }

  if (dual.permuted())
    applyPermutation<Dimension>(dual.invPermutation(), r_map,
                                dual.W.majorIndex().innerOffsetsData());

  return res;
}

template <unsigned Dimension, typename Method, typename MatrixT>
static double eval(const DualFrictionProblem<Dimension> &dual,
                   const ConstrainedSolverBase<Method, MatrixT> &gs,
                   const double *r_data, const bool staticProblem,
                   const bool herschelBulkleyProblem) {
  Eigen::VectorXd r = Eigen::VectorXd::Map(r_data, dual.W.rows());

  if (dual.permuted())
    applyPermutation<Dimension>(dual.permutation(), r,
                                dual.W.majorIndex().innerOffsetsData());

  const Eigen::VectorXd u = dual.W * r + dual.b;

  double res;

  if (herschelBulkleyProblem) {
    res = gs.eval(typename DualFrictionProblem<Dimension>::HBLawType(
                      dual.W.rowsOfBlocks(), dual.yields.data(),
                      dual.etas.data(), dual.powers.data()),
                  u, r);
  } else {
    res = staticProblem
              ? gs.eval(typename DualFrictionProblem<Dimension>::SOCLawType(
                            dual.W.rowsOfBlocks(), dual.mu.data()),
                        u, r)
              : gs.eval(typename DualFrictionProblem<Dimension>::CoulombLawType(
                            dual.W.rowsOfBlocks(), dual.mu.data()),
                        u, r);
  }

  return res;
}

template <unsigned Dimension, typename Method, typename MatrixT>
static double solveCadoux(const DualFrictionProblem<Dimension> &problem,
                          ConstrainedSolverBase<Method, MatrixT> &minimizer,
                          double *r, const unsigned cadouxIterations,
                          const Signal<unsigned, double> *callback) {
  Eigen::Map<Eigen::VectorXd> r_map(r, problem.W.rows());

  if (problem.permuted())
    applyPermutation<Dimension>(problem.permutation(), r_map,
                                problem.W.majorIndex().innerOffsetsData());

  const double res =
      solveCadoux<Dimension>(problem.W, problem.b.data(), problem.mu.data(),
                             minimizer, r, cadouxIterations, callback);

  if (problem.permuted())
    applyPermutation<Dimension>(problem.invPermutation(), r_map,
                                problem.W.majorIndex().innerOffsetsData());

  return res;
}

template <unsigned Dimension, typename Method, typename MatrixT>
static double solveCadouxVel(const DualFrictionProblem<Dimension> &problem,
                             ConstrainedSolverBase<Method, MatrixT> &minimizer,
                             double *u, const unsigned cadouxIterations,
                             const Signal<unsigned, double> *callback) {
  Eigen::Map<Eigen::VectorXd> u_map(u, problem.W.rows());

  if (problem.permuted())
    applyPermutation<Dimension>(problem.permutation(), u_map,
                                problem.W.majorIndex().innerOffsetsData());

  const double res =
      solveCadouxVel<Dimension>(problem.W, problem.b.data(), problem.mu.data(),
                                minimizer, u, cadouxIterations, callback);

  if (problem.permuted())
    applyPermutation<Dimension>(problem.invPermutation(), u_map,
                                problem.W.majorIndex().innerOffsetsData());

  return res;
}

}  // namespace friction_problem

}  // namespace bogus
