/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "../Core/BlockSolvers/ADMM.impl.hpp"
#include "../Core/BlockSolvers/ProductGaussSeidel.impl.hpp"
#include "FrictionProblem.impl.hpp"

namespace bogus {

template <unsigned Dimension>
void PrimalFrictionProblem<Dimension>::computeMInv() {
  // M^-1
  MInv.cloneStructure(M);
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
  for (std::ptrdiff_t i = 0; i < (std::ptrdiff_t)M.nBlocks(); ++i) {
    MInv.block(i).compute(M.block(i));
  }
}

template <unsigned Dimension>
void PrimalFrictionProblem<Dimension>::computeDiagMInv() {
  // M^-1
  DiagMInv.cloneStructure(DiagM);
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
  for (std::ptrdiff_t i = 0; i < (std::ptrdiff_t)DiagM.nBlocks(); ++i) {
    DiagMInv.block(i) = DiagM.block(i).inverse();
  }
}

template <unsigned Dimension>
double PrimalFrictionProblem<Dimension>::solveWith(
    DiagProductGaussSeidelType& pgs, double* r,
    const bool staticProblem) const {
  // b = E' w - H M^-1 f
  Eigen::VectorXd b;

  b = E.transpose() * w - H * (DiagMInv * (f + H.transpose() * rc));

  Eigen::VectorXd::MapType r_map(r, H.rows());

  pgs.setMatrix(H);

  pgs.setDiagonal(DiagMInv);

  if (staticProblem) {
    bogus::SOCLaw<Dimension, double, false> law(H.rowsOfBlocks(), mu.data());
    return pgs.solve(law, b, r_map);
  } else {
    bogus::SOCLaw<Dimension, double, true> law(H.rowsOfBlocks(), mu.data());
    return pgs.solve(law, b, r_map);
  }
}

template <unsigned Dimension>
double PrimalFrictionProblem<Dimension>::solveWith(
    ProductGaussSeidelType& pgs, double* r, const bool staticProblem) const {
  // b = E' w - H M^-1 f
  Eigen::VectorXd b;

  b = E.transpose() * w - H * (MInv * (f + H.transpose() * rc));

  Eigen::VectorXd::MapType r_map(r, H.rows());

  pgs.setMatrix(H);

  pgs.setDiagonal(MInv);

  if (staticProblem) {
    bogus::SOCLaw<Dimension, double, false> law(H.rowsOfBlocks(), mu.data());
    return pgs.solve(law, b, r_map);
  } else {
    bogus::SOCLaw<Dimension, double, true> law(H.rowsOfBlocks(), mu.data());
    return pgs.solve(law, b, r_map);
  }
}

template <unsigned Dimension>
double PrimalFrictionProblem<Dimension>::solveWith(ADMMType& admm,
                                                   double lambda, double* v,
                                                   double* r) const {
  const Eigen::VectorXd fc = f + H.transpose() * rc;
  const Eigen::VectorXd Ew = E.transpose() * w;

  Eigen::VectorXd::MapType r_map(r, H.rows());
  Eigen::VectorXd::MapType v_map(v, H.cols());

  bogus::QuadraticProxOp<MInvType> prox(MInv, lambda, fc);

  Eigen::ArrayXd inv_mu = 1. / mu.array();
  bogus::SOCLaw<Dimension, double, false> law(inv_mu.rows(), inv_mu.data());

  admm.setMatrix(H);
  return admm.solve(law, prox, Ew, v_map, r_map);
}

template <unsigned Dimension>
double PrimalFrictionProblem<Dimension>::solveWith(
    DualAMAType& dama, double* v, double* r, const bool staticProblem) const {
  const Eigen::VectorXd fc = f + H.transpose() * rc;
  const Eigen::VectorXd Ew = E.transpose() * w;

  Eigen::VectorXd::MapType r_map(r, H.rows());
  Eigen::VectorXd::MapType v_map(v, H.cols());

  Eigen::ArrayXd inv_mu = 1. / mu.array();
  bogus::SOCLaw<Dimension, double, false> law(inv_mu.rows(), inv_mu.data());

  dama.setMatrix(H);

  if (staticProblem) {
    bogus::SOCLaw<Dimension, double, false> law(H.rowsOfBlocks(), mu.data());
    return dama.solve(law, M, fc, Ew, v_map, r_map);
  } else {
    bogus::SOCLaw<Dimension, double, true> law(H.rowsOfBlocks(), mu.data());
    return dama.solve(law, M, fc, Ew, v_map, r_map);
  }
}

#ifdef BOGUS_INSTANTIATE_2D_SOC
template struct PrimalFrictionProblem<2u>;
#endif

#ifdef BOGUS_INSTANTIATE_3D_SOC
template struct PrimalFrictionProblem<3u>;
#endif

#ifdef BOGUS_INSTANTIATE_DYNAMIC_SOC
template struct PrimalFrictionProblem<Eigen::Dynamic>;
#endif

}  // namespace bogus
