/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_CONSTRAINED_SOLVER_BASE_IMPL_HPP
#define BOGUS_CONSTRAINED_SOLVER_BASE_IMPL_HPP

#include "../Block/Access.hpp"
#include "../Block/BlockMatrixBase.hpp"
#include "ConstrainedSolverBase.hpp"

namespace bogus {

template <typename Derived, typename BlockMatrixType>
template <typename NSLaw, typename RhsT, typename ResT>
typename ConstrainedSolverBase<Derived, BlockMatrixType>::Scalar
ConstrainedSolverBase<Derived, BlockMatrixType>::eval(const NSLaw &law,
                                                      const ResT &y,
                                                      const RhsT &x) const {
  const Index dimension = BlockProblemTraits::dimension;

  const Segmenter<dimension, const RhsT, typename BlockMatrixType::Index>
      xSegmenter(x, m_matrix->rowOffsets());
  const Segmenter<dimension, const ResT, typename BlockMatrixType::Index>
      ySegmenter(y, m_matrix->rowOffsets());

  typedef typename BlockMatrixTraits<BlockMatrixType>::Index Index;

  const Index n = m_matrix->rowsOfBlocks();

  Scalar err = 0., lres;
  typename NSLaw::Traits::Vector lx, ly;

  if (m_useInfinityNorm) {
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel private(lx, ly, lres)
#endif
    {
      lres = 0.;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp for
#endif
      for (Index i = 0; i < n; ++i) {
        lx = xSegmenter[i];
        ly = ySegmenter[i];
        lres = std::max(law.eval(i, lx, ly, m_scaling[i]), lres);
      }

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp critical
#endif
      err = std::max(err, lres);
    }

    return err;

  } else {
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private(lx, ly, lres) reduction(+ : err)
#endif
    for (Index i = 0; i < n; ++i) {
      lx = xSegmenter[i];
      ly = ySegmenter[i];
      lres = law.eval(i, lx, ly, m_scaling[i]);
      err += lres;
    }

    return err / (1 + n);
  }
}

namespace block_solvers_impl {

template <typename MatrixT>
typename MatrixT::Scalar estimate_block_scaling(const MatrixT &block) {
  return std::max((typename MatrixT::Scalar)1, block.trace() / block.rows());
}

template <typename Derived>
void estimate_row_scaling(const BlockObjectBase<Derived> &,
                          typename BlockObjectBase<Derived>::Scalar *) {}

template <typename Derived>
void estimate_row_scaling(const BlockMatrixBase<Derived> &mat,
                          typename Derived::Scalar *scalings) {
  typedef BlockMatrixTraits<Derived> BlockTraits;
  typedef typename BlockTraits::BlockType LocalMatrixType;
  typedef MatrixTraits<LocalMatrixType> Traits;

  // For square matrices, estimate diag block
  if (mat.rows() == mat.cols() && mat.rowsOfBlocks() == mat.colsOfBlocks()) {
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
    for (typename BlockTraits::Index row = 0; row < mat.rowsOfBlocks(); ++row) {
      const typename BlockTraits::BlockPtr ptr = mat.diagonalBlockPtr(row);
      if (ptr != Derived::InvalidBlockPtr) {
        scalings[row] =
            estimate_block_scaling(Traits::asConstMatrix(mat.block(ptr)));
      }
    }
  }
}

}  // namespace block_solvers_impl

template <typename Derived, typename BlockMatrixType>
void ConstrainedSolverBase<Derived, BlockMatrixType>::updateScalings() {
  if (!m_matrix) {
    return;
  }

  const Index n = m_matrix->rowsOfBlocks();
  m_scaling.setOnes(n);

  block_solvers_impl::estimate_row_scaling(m_matrix->derived(),
                                           m_scaling.data());
}

template <typename Derived, typename BlockMatrixType>
template <typename NSLaw, typename VectorT>
void ConstrainedSolverBase<Derived, BlockMatrixType>::projectOnConstraints(
    const NSLaw &law, VectorT &x) const {
  Segmenter<NSLaw::dimension, VectorT, typename BlockMatrixType::Index>
      xSegmenter(x, m_matrix->rowOffsets());

  const Index n = m_matrix->rowsOfBlocks();
  typename NSLaw::Traits::Vector lx;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private(lx)
#endif
  for (Index i = 0; i < n; ++i) {
    lx = xSegmenter[i];
    law.projectOnConstraint(i, lx);
    xSegmenter[i] = lx;
  }
}

template <typename Derived, typename BlockMatrixType>
template <typename NSLaw, typename RhsT, typename ResT>
void ConstrainedSolverBase<Derived, BlockMatrixType>::dualityCOV(
    const NSLaw &law, const RhsT &u, ResT &s) const {
  const Segmenter<NSLaw::dimension, const RhsT, typename BlockMatrixType::Index>
      uSegmenter(u, m_matrix->rowOffsets());
  Segmenter<NSLaw::dimension, ResT, typename BlockMatrixType::Index> sSegmenter(
      s, m_matrix->rowOffsets());

  const Index n = m_matrix->rowsOfBlocks();
  typename NSLaw::Traits::Vector ls;

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for private(ls)
#endif
  for (Index i = 0; i < n; ++i) {
    law.dualityCOV(i, uSegmenter[i], ls);
    sSegmenter[i] = ls;
  }
}

}  // namespace bogus

#endif
