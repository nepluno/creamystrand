/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_PRODUCT_GAUSS_SEIDEL_IMPL_HPP
#define BOGUS_PRODUCT_GAUSS_SEIDEL_IMPL_HPP

#include "GaussSeidelBase.impl.hpp"
#include "ProductGaussSeidel.hpp"

#ifndef BOGUS_DONT_PARALLELIZE
#include <omp.h>
#endif

namespace bogus {

namespace block_solvers_impl {

// Utils

template <typename Type>
void DiagonalMatrixWrapper<Type, false>::computeBlockIndices() {
  // Stores the index of each diagonal block of (*matrixPtr)
  const Index rows = m_matrixPtr->rowsOfBlocks();
  m_blockIndices.resize(rows);

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
  for (Index i = 0; i < rows; ++i) {
    m_blockIndices[i] = m_matrixPtr->diagonalBlockPtr(i);
  }
}

// DMtStorage w/o precomputing

template <typename MType, typename DType, bool Precompute>
template <typename Rhs, typename Intermediate, typename Res>
void DMtStorage<MType, DType, Precompute>::multiply(const Rhs& rhs,
                                                    Intermediate& itm,
                                                    Res& res) const {
  m_M->template multiply<true>(rhs, itm, 1, 0);
  res += (*m_M) * m_D->get() * itm;
}

template <typename MType, typename DType, bool Precompute>
template <typename Rhs, typename Res>
void DMtStorage<MType, DType, Precompute>::colMultiply(
    typename MType::Index col, const Rhs& rhs, Res& itm) const {
  m_M->template colMultiply<true>(col, rhs, itm);
}

template <typename MType, typename DType, bool Precompute>
template <typename Rhs, typename Res>
void DMtStorage<MType, DType, Precompute>::rowMultiply(
    typename MType::Index row, const Rhs& itm, Res& res) const {
  m_M->template rowMultiplyPrecompose<false>(row, itm, res, m_D->asArray());
}

// DMtStorgae w/ precomputing

template <typename MType, typename DType>
void DMtStorage<MType, DType, true>::compute(const MType& M, const DType& D) {
  m_DMt = D.get() * M.transpose();
  m_M = &M;
}

template <typename MType, typename DType>
template <typename Rhs, typename Intermediate, typename Res>
void DMtStorage<MType, DType, true>::multiply(const Rhs& rhs, Intermediate& itm,
                                              Res& res) const {
  m_DMt.template multiply<false>(rhs, itm, 1, 0);
  res += (*m_M) * itm;
}

template <typename MType, typename DType>
template <typename Rhs, typename Res>
void DMtStorage<MType, DType, true>::colMultiply(typename MType::Index col,
                                                 const Rhs& rhs,
                                                 Res& itm) const {
  m_DMt.template colMultiply<false>(col, rhs, itm);
}

template <typename MType, typename DType>
template <typename Rhs, typename Res>
void DMtStorage<MType, DType, true>::rowMultiply(typename MType::Index row,
                                                 const Rhs& itm,
                                                 Res& res) const {
  m_M->template rowMultiply<false>(row, itm, res);
}

template <typename D, typename T>
struct SelfProductAccumulator {
  const D& diag;
  T& acc;
  SelfProductAccumulator(const D& diag_, T& mat) : diag(diag_), acc(mat) {}

  template <typename Index, typename Matrix>
  void operator()(const Index i, const Matrix& block) {
    acc += block * (diag[i] * block.transpose());
  }
};

template <typename D, typename T>
SelfProductAccumulator<D, T> accumulate(const D& diag, T& res) {
  return SelfProductAccumulator<D, T>(diag, res);
}

}  // namespace block_solvers_impl

template <typename BlockMatrixType, typename DiagonalType, bool PrecomputeDMt>
ProductGaussSeidel<BlockMatrixType, DiagonalType, PrecomputeDMt>&
ProductGaussSeidel<BlockMatrixType, DiagonalType, PrecomputeDMt>::setMatrix(
    const BlockObjectBase<BlockMatrixType>& M) {
  m_matrix = &M;

  updateLocalMatrices();

  return *this;
}

template <typename BlockMatrixType, typename DiagonalType, bool PrecomputeDMt>
ProductGaussSeidel<BlockMatrixType, DiagonalType, PrecomputeDMt>&
ProductGaussSeidel<BlockMatrixType, DiagonalType, PrecomputeDMt>::setDiagonal(
    const DiagonalType& diagonal) {
  m_diagonal = DiagWrapper(diagonal);

  updateLocalMatrices();

  return *this;
}

template <typename BlockMatrixType, typename DiagonalType, bool PrecomputeDMt>
void ProductGaussSeidel<BlockMatrixType, DiagonalType,
                        PrecomputeDMt>::updateLocalMatrices() {
  if (!(m_matrix && m_diagonal.valid())) return;

  m_DMt.compute(m_matrix->derived(), m_diagonal);

  const Index n = m_matrix->rowsOfBlocks();
  m_localMatrices.resize(n);

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
  for (Index i = 0; i < n; ++i) {
    m_localMatrices[i].setZero();

    // TODO: use precomputed DMt if available
    Base::iterableMatrix().eachBlockOfRow(
        i, block_solvers_impl::accumulate(m_diagonal.asArray(),
                                          m_localMatrices[i]));
  }

  Base::processLocalMatrices();
}

template <typename BlockMatrixType, typename DiagonalType, bool PrecomputeDMt>
template <typename NSLaw, typename VecT, typename ResT>
void ProductGaussSeidel<BlockMatrixType, DiagonalType, PrecomputeDMt>::
    innerLoop(bool parallelize, const NSLaw& law, const VecT& b,
              std::vector<unsigned char>& skip, Scalar& ndxRef, VecT& Mx,
              ResT& x) const {
  typedef typename NSLaw::Traits LocalProblemTraits;

  Segmenter<NSLaw::dimension, ResT, typename BlockMatrixType::Index> xSegmenter(
      x, m_matrix->rowOffsets());
  const Segmenter<NSLaw::dimension, const VecT, typename BlockMatrixType::Index>
      bSegmenter(b, m_matrix->rowOffsets());

  const Scalar absSkipTol = std::min(m_skipTol, m_tol);
  const Scalar absSkipIters =
      std::min(m_skipIters, (unsigned)std::sqrt((Scalar)skip.size()));

  const std::ptrdiff_t n = m_matrix->rowsOfBlocks();

#ifdef BOGUS_DONT_PARALLELIZE
  (void)parallelize;
#else
#pragma omp parallel if (parallelize)
  {
#endif
  typename LocalProblemTraits::Vector lb, lx, ldx;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp for
#endif
  for (std::ptrdiff_t i = 0; i < n; ++i) {
    if (skip[i]) {
      --skip[i];
      continue;
    }

    lx = xSegmenter[i];
    lb = bSegmenter[i] - m_localMatrices[i] * lx;

    m_DMt.rowMultiply(i, Mx, lb);

    ldx = -lx;

    const bool ok = law.solveLocal(i, m_localMatrices[i], lb, lx, m_scaling[i]);
    ldx += lx;

    if (!ok) {
      ldx *= .5;
    }
    xSegmenter[i] += ldx;

    m_DMt.colMultiply(i, ldx, Mx);

    const Scalar nx2 = m_scaling[i] * m_scaling[i] * lx.squaredNorm();
    const Scalar ndx2 = m_scaling[i] * m_scaling[i] * ldx.squaredNorm();
    // Thread-safety:
    // we do not care if we lose an update of ndxRef to a data race,
    // but we need its bytes to be consistent.
    // Ok on x86(_64) if ndxRef is aligned as assignment is atomic
    if (ndx2 > ndxRef) ndxRef = ndx2;

    if (std::min(nx2, ndx2) < absSkipTol ||
        ndx2 < m_skipTol * std::min(nx2, ndxRef)) {
      skip[i] = absSkipIters;
    }
  }

#ifndef BOGUS_DONT_PARALLELIZE
}
#endif

}  // namespace bogus

template <typename BlockMatrixType, typename DiagonalType, bool PrecomputeDMt>
template <typename NSLaw, typename RhsT, typename ResT>
typename ProductGaussSeidel<BlockMatrixType, DiagonalType,
                            PrecomputeDMt>::Scalar
ProductGaussSeidel<BlockMatrixType, DiagonalType, PrecomputeDMt>::solve(
    const NSLaw& law, const RhsT& b, ResT& x, bool tryZeroAsWell) const {
  const bogus::Zero<Scalar> zero;
  return solveWithLinearConstraints(law, zero, b, x, tryZeroAsWell, 0);
}

template <typename BlockMatrixType, typename DiagonalType, bool PrecomputeDMt>
template <typename NSLaw, typename RhsT, typename ResT, typename LSDerived,
          typename HDerived>
typename ProductGaussSeidel<BlockMatrixType, DiagonalType,
                            PrecomputeDMt>::Scalar
ProductGaussSeidel<BlockMatrixType, DiagonalType, PrecomputeDMt>::
    solveWithLinearConstraints(const NSLaw& law,
                               const BlockObjectBase<LSDerived>& Cinv,
                               const BlockObjectBase<HDerived>& H,
                               const Scalar alpha, const RhsT& b, const RhsT& c,
                               ResT& x, bool tryZeroAsWell,
                               unsigned solveEvery) const {
  const typename GlobalProblemTraits::DynVector k = b + H * Cinv * c;

  return solveWithLinearConstraints(law, -alpha * H * Cinv * H.transpose(), k,
                                    x, tryZeroAsWell, solveEvery);
}

template <typename BlockMatrixType, typename DiagonalType, bool PrecomputeDMt>
template <typename NSLaw, typename RhsT, typename ResT, typename WDerived>
typename ProductGaussSeidel<BlockMatrixType, DiagonalType,
                            PrecomputeDMt>::Scalar
ProductGaussSeidel<BlockMatrixType, DiagonalType, PrecomputeDMt>::
    solveWithLinearConstraints(const NSLaw& law,
                               const BlockObjectBase<WDerived>& W,
                               const RhsT& b, ResT& x, bool tryZeroAsWell,
                               unsigned solveEvery) const {
  assert(m_matrix);
  assert(m_diagonal.valid());
  assert(solveEvery == 0 || 0 == m_evalEvery % solveEvery);

  typename GlobalProblemTraits::DynVector Mx(m_matrix->cols()), y, x_best;

  typename GlobalProblemTraits::DynVector w = b;
  W.template multiply<false>(x, w, 1, 1);

  Scalar err_best = std::numeric_limits<Scalar>::max();

  y = w;
  m_DMt.multiply(x, Mx, y);
  Base::evalAndKeepBest(law, x, y, x_best, err_best);

  if (tryZeroAsWell && Base::tryZero(law, b, x, x_best, err_best)) {
    w = b;
    Mx.setZero();
  }

  this->m_callback.trigger(0, err_best);

  const Index n = m_matrix->rowsOfBlocks();

  WithMaxThreads wmt(m_maxThreads);
  const int newMaxThreads = wmt.nThreads();
  const bool parallelize =
      (newMaxThreads != 1 && n > newMaxThreads * newMaxThreads);

  std::vector<unsigned char> skip(n, 0);
  Scalar ndxRef = 0;  // Reference step size

  unsigned GSIter;
  for (GSIter = 1; GSIter <= m_maxIters; ++GSIter) {
    innerLoop(parallelize, law, w, skip, ndxRef, Mx, x);

    if (solveEvery > 0 && 0 == (GSIter % solveEvery)) {
      w = b;
      W.template multiply<false>(x, w, 1, 1);
    }

    if (0 == (GSIter % m_evalEvery)) {
      y = w;
      m_DMt.multiply(x, Mx, y);
      const Scalar err = Base::evalAndKeepBest(law, x, y, x_best, err_best);

      this->m_callback.trigger(GSIter, err);

      if (err < m_tol) {
        break;
      }

      ndxRef /= m_evalEvery;
    }
  }

  if (GSIter > m_maxIters) x = x_best;

  return err_best;
}

}  // namespace bogus

#endif
