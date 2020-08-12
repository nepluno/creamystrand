/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "../Core/BlockSolvers/GaussSeidel.impl.hpp"
#include "../Core/BlockSolvers/ProjectedGradient.impl.hpp"
#include "FrictionProblem.impl.hpp"

namespace bogus {

// DualFrictionProblem

template <unsigned Dimension>
void DualFrictionProblem<Dimension>::computeFrom(
    const PrimalFrictionProblem<Dimension> &primal,
    const bool diagonalProblem) {
  if (diagonalProblem) {
    // W
    typename PrimalFrictionProblem<Dimension>::HType SH =
        primal.H * primal.DiagMInv;

    W = primal.H * SH.transpose();

    // M^-1 f, b
    b = primal.E.transpose() * primal.w -
        primal.H * (primal.DiagMInv * primal.f) - W * primal.rc;
  } else {
    // W
    W = primal.H * (primal.MInv * primal.H.transpose());

    // M^-1 f, b
    b = primal.E.transpose() * primal.w - primal.H * (primal.MInv * primal.f) -
        W * primal.rc;
  }

  mu = primal.mu;

  yields = primal.yields;

  etas = primal.etas;

  powers = primal.powers;
}

template <unsigned Dimension>
double DualFrictionProblem<Dimension>::solveWith(
    GaussSeidelType &gs, double *r, const bool staticProblem,
    const bool herschelBulkleyProblem) const {
  gs.setMatrix(W);

  return friction_problem::solve(*this, gs, r, staticProblem,
                                 herschelBulkleyProblem);
}

template <unsigned Dimension>
double DualFrictionProblem<Dimension>::solveWith(
    ProjectedGradientType &pg, double *r, const bool staticProblem,
    const bool herschelBulkleyProblem) const {
  pg.setMatrix(W);

  return friction_problem::solve(*this, pg, r, staticProblem,
                                 herschelBulkleyProblem);
}

template <unsigned Dimension>
double DualFrictionProblem<Dimension>::evalWith(
    const GaussSeidelType &gs, const double *r, const bool staticProblem,
    const bool herschelBulkleyProblem) const {
  return friction_problem::eval(*this, gs, r, staticProblem,
                                herschelBulkleyProblem);
}

template <unsigned Dimension>
double DualFrictionProblem<Dimension>::evalWith(
    const ProjectedGradientType &gs, const double *r, const bool staticProblem,
    const bool herschelBulkleyProblem) const {
  return friction_problem::eval(*this, gs, r, staticProblem,
                                herschelBulkleyProblem);
}

template <unsigned Dimension>
double DualFrictionProblem<Dimension>::solveCadoux(
    GaussSeidelType &gs, double *r, const unsigned cadouxIterations,
    const Signal<unsigned, double> *callback) const {
  gs.setMatrix(W);

  return friction_problem::solveCadoux(*this, gs, r, cadouxIterations,
                                       callback);
}

template <unsigned Dimension>
double DualFrictionProblem<Dimension>::solveCadoux(
    ProjectedGradientType &pg, double *r, const unsigned cadouxIterations,
    const Signal<unsigned, double> *callback) const {
  pg.setMatrix(W);

  return friction_problem::solveCadoux(*this, pg, r, cadouxIterations,
                                       callback);
}

template <unsigned Dimension>
double DualFrictionProblem<Dimension>::solveCadouxVel(
    GaussSeidelType &gs, double *u, const unsigned cadouxIterations,
    const Signal<unsigned, double> *callback) const {
  gs.setMatrix(W);

  return friction_problem::solveCadouxVel(*this, gs, u, cadouxIterations,
                                          callback);
}

template <unsigned Dimension>
double DualFrictionProblem<Dimension>::solveCadouxVel(
    ProjectedGradientType &pg, double *u, const unsigned cadouxIterations,
    const Signal<unsigned, double> *callback) const {
  pg.setMatrix(W);

  return friction_problem::solveCadouxVel(*this, pg, u, cadouxIterations,
                                          callback);
}

template <unsigned Dimension>
void DualFrictionProblem<Dimension>::applyPermutation(
    const std::vector<std::size_t> &permutation) {
  assert(!permuted());

  m_permutation = permutation;

  m_invPermutation.resize(m_permutation.size());
  for (std::size_t i = 0; i < m_permutation.size(); ++i)
    m_invPermutation[m_permutation[i]] = i;

  W.applyPermutation(data_pointer(m_permutation));
  friction_problem::applyPermutation<Dimension>(m_permutation, b,
                                                W.colOffsets());
  bogus::applyPermutation(m_permutation.size(), data_pointer(m_permutation),
                          mu);
}

template <unsigned Dimension>
void DualFrictionProblem<Dimension>::undoPermutation() {
  if (!permuted()) return;

  W.applyPermutation(data_pointer(m_invPermutation));
  friction_problem::applyPermutation<Dimension>(m_invPermutation, b,
                                                W.colOffsets());
  bogus::applyPermutation(m_invPermutation.size(),
                          data_pointer(m_invPermutation), mu);

  m_permutation.clear();
}

#ifdef BOGUS_INSTANTIATE_2D_SOC
template struct DualFrictionProblem<2u>;
#endif

#ifdef BOGUS_INSTANTIATE_3D_SOC
template struct DualFrictionProblem<3u>;
#endif

#ifdef BOGUS_INSTANTIATE_DYNAMIC_SOC
template struct DualFrictionProblem<Eigen::Dynamic>;
#endif

}  // namespace bogus
