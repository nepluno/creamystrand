/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_GAUSS_SEIDEL_BASE_IMPL_HPP
#define BOGUS_BLOCK_GAUSS_SEIDEL_BASE_IMPL_HPP

#include "../Utils/Threads.hpp"
#include "ConstrainedSolverBase.impl.hpp"
#include "GaussSeidelBase.hpp"

namespace bogus {

template <typename GaussSeidelImpl, typename BlockMatrixType>
void GaussSeidelBase<GaussSeidelImpl, BlockMatrixType>::processLocalMatrices() {
  Base::updateScalings();

  if (!m_matrix) return;

  const Index n = m_matrix->rowsOfBlocks();
  m_regularization.resize(n);

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
  for (Index i = 0; i < n; ++i) {
    m_scaling[i] =
        block_solvers_impl::estimate_block_scaling(m_localMatrices[i]);

    if (m_autoRegularization > 0.) {
      m_regularization(i) =
          std::max(0., m_autoRegularization -
                           m_localMatrices[i].eigenvalues().real().minCoeff());
      m_localMatrices[i].diagonal().array() += m_regularization(i);
    } else
      m_regularization(i) = 0.;
  }
}

template <typename GaussSeidelImpl, typename BlockMatrixType>
template <typename NSLaw, typename ResT>
typename GaussSeidelBase<GaussSeidelImpl, BlockMatrixType>::Scalar
GaussSeidelBase<GaussSeidelImpl, BlockMatrixType>::evalAndKeepBest(
    const NSLaw &law, const ResT &x,
    const typename GlobalProblemTraits::DynVector &y,
    typename GlobalProblemTraits::DynVector &x_best, Scalar &err_best) const {
  const Scalar err = Base::eval(law, y, x);

  if (err < err_best) {
    x_best = x;
    err_best = err;
  }

  return err;
}

template <typename GaussSeidelImpl, typename BlockMatrixType>
template <typename NSLaw, typename RhsT, typename ResT>
bool GaussSeidelBase<GaussSeidelImpl, BlockMatrixType>::tryZero(
    const NSLaw &law, const RhsT &b, ResT &x,
    typename GlobalProblemTraits::DynVector &x_best, Scalar &err_best) const {
  x.setZero();
  const Scalar err_zero = Base::eval(law, b, x);

  if (err_zero < err_best) {
    err_best = err_zero;
    x_best.setZero();
    return true;
  }

  x = x_best;
  return false;
}

}  // namespace bogus

#endif
