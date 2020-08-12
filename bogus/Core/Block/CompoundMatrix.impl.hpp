/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_COMPOUND_BLOCK_MATRIX_IMPL_HPP
#define BOGUS_COMPOUND_BLOCK_MATRIX_IMPL_HPP

#include "CompoundMatrix.hpp"

namespace bogus {

namespace mv_impl {
template <bool ColWiseMult, bool DoTranspose>
struct CompoundMultiplier {  // colwise compound not transposed or rowwise
                             // transposed

  template <typename MatrixT1, typename MatrixT2, typename RhsT, typename ResT,
            typename Scalar>
  static void multiply(
      const CompoundBlockMatrix<!DoTranspose, MatrixT1, MatrixT2>& mat,
      const RhsT& rhs, ResT& res, Scalar alpha = 1, Scalar beta = 0) {
    const Segmenter<
        internal::DYNAMIC, const RhsT,
        typename CompoundBlockMatrix<!DoTranspose, MatrixT1, MatrixT2>::Index>
        rhsSeg(rhs, mat.compoundOffsets());
    mat.first().template multiply<DoTranspose>(rhsSeg[0], res, alpha, beta);
    mat.second().template multiply<DoTranspose>(rhsSeg[1], res, alpha, 1);
  }

  template <typename MatrixT1, typename MatrixT2, typename RhsT, typename ResT,
            typename Op>
  static void rowMultiplyPrecompose(
      const CompoundBlockMatrix<!DoTranspose, MatrixT1, MatrixT2>& mat,
      typename MatrixT1::Index row, const RhsT& rhs, ResT& res, const Op& op) {
    const Segmenter<
        internal::DYNAMIC, const RhsT,
        typename CompoundBlockMatrix<!DoTranspose, MatrixT1, MatrixT2>::Index>
        rhsSeg(rhs, mat.compoundOffsets());
    mat.first().template rowMultiplyPrecompose<DoTranspose>(row, rhsSeg[0], res,
                                                            op);
    mat.second().template rowMultiplyPrecompose<DoTranspose>(row, rhsSeg[1],
                                                             res, op);
  }

  template <typename MatrixT1, typename MatrixT2, typename RhsT, typename ResT,
            typename Op>
  static void colMultiplyPostcompose(
      const CompoundBlockMatrix<!DoTranspose, MatrixT1, MatrixT2>& mat,
      typename MatrixT1::Index col, const RhsT& rhs, ResT& res, const Op& op) {
    typename MatrixT1::Index off = col - mat.secondBegin();
    if (off < 0)
      mat.first().template colMultiplyPostcompose<DoTranspose>(col, rhs, res,
                                                               op);
    else
      mat.second().template colMultiplyPostcompose<DoTranspose>(off, rhs, res,
                                                                op);
  }
};
template <bool DoTranspose>
struct CompoundMultiplier<false, DoTranspose> {  // colwise compound transposed
                                                 // or rowwise not transposed

  template <typename MatrixT1, typename MatrixT2, typename RhsT, typename ResT,
            typename Scalar>
  static void multiply(
      const CompoundBlockMatrix<DoTranspose, MatrixT1, MatrixT2>& mat,
      const RhsT& rhs, ResT& res, Scalar alpha = 1, Scalar beta = 0) {
    typedef Segmenter<
        internal::DYNAMIC, ResT,
        typename CompoundBlockMatrix<DoTranspose, MatrixT1, MatrixT2>::Index>
        SegmenterType;
    SegmenterType resSeg(res, mat.compoundOffsets());

    typename SegmenterType::ReturnType first(resSeg[0]);
    mat.first().template multiply<DoTranspose>(rhs, first, alpha, beta);
    typename SegmenterType::ReturnType second(resSeg[1]);
    mat.second().template multiply<DoTranspose>(rhs, second, alpha, beta);
  }

  template <typename MatrixT1, typename MatrixT2, typename RhsT, typename ResT,
            typename Op>
  static void rowMultiplyPrecompose(
      const CompoundBlockMatrix<DoTranspose, MatrixT1, MatrixT2>& mat,
      typename MatrixT1::Index row, const RhsT& rhs, ResT& res, const Op& op) {
    typename MatrixT1::Index off = row - mat.secondBegin();
    if (off < 0)
      mat.first().template rowMultiplyPrecompose<DoTranspose>(row, rhs, res,
                                                              op);
    else
      mat.second().template rowMultiplyPrecompose<DoTranspose>(off, rhs, res,
                                                               op);
  }

  template <typename MatrixT1, typename MatrixT2, typename RhsT, typename ResT,
            typename Op>
  static void colMultiplyPostcompose(
      const CompoundBlockMatrix<DoTranspose, MatrixT1, MatrixT2>& mat,
      typename MatrixT1::Index col, const RhsT& rhs, ResT& res, const Op& op) {
    typedef Segmenter<
        internal::DYNAMIC, ResT,
        typename CompoundBlockMatrix<DoTranspose, MatrixT1, MatrixT2>::Index>
        SegmenterType;
    SegmenterType resSeg(res, mat.compoundOffsets());

    typename SegmenterType::ReturnType first(resSeg[0]);
    mat.first().template colMultiplyPostcompose<DoTranspose>(col, rhs, first,
                                                             op);
    typename SegmenterType::ReturnType second(resSeg[1]);
    mat.second().template colMultiplyPostcompose<DoTranspose>(col, rhs, second,
                                                              op);
  }
};
}  // namespace mv_impl

// Column-wise compound matrix

template <bool ColWise, typename MatrixT1, typename MatrixT2>
CompoundBlockMatrix<ColWise, MatrixT1, MatrixT2>::CompoundBlockMatrix(
    const IterableBlockObject<MatrixT1>& first,
    const IterableBlockObject<MatrixT2>& second)
    : m_first(first), m_second(second) {
  // Asserts
  // TODO check compatibility of block sizes

  // Macro offsets
  m_compoundOffsets[0] = 0;
  m_compoundOffsets[1] = m_first.cols();
  m_compoundOffsets[2] = cols();

  // Underlying Offsets
  m_offsets.resize(colsOfBlocks() + 1);
  std::copy(m_first.colOffsets(),
            m_first.colOffsets() + m_first.colsOfBlocks() + 1,
            m_offsets.begin());

  const Index secStart = m_first.colsOfBlocks() + 1;
  for (Index col = 0; col < m_second.colsOfBlocks(); ++col) {
    m_offsets[secStart + col] =
        m_second.colOffsets()[col + 1] + m_compoundOffsets[1];
  }
}

template <bool ColWise, typename MatrixT1, typename MatrixT2>
template <bool DoTranspose, typename RhsT, typename ResT>
void CompoundBlockMatrix<ColWise, MatrixT1, MatrixT2>::multiply(
    const RhsT& rhs, ResT& res, Scalar alpha, Scalar beta) const {
  mv_impl::CompoundMultiplier<!DoTranspose, DoTranspose>::multiply(
      *this, rhs, res, alpha, beta);
}

template <bool ColWise, typename MatrixT1, typename MatrixT2>
template <bool DoTranspose, typename RhsT, typename ResT, typename PreOp>
void CompoundBlockMatrix<ColWise, MatrixT1, MatrixT2>::rowMultiplyPrecompose(
    const Index row, const RhsT& rhs, ResT& res, const PreOp& op) const {
  mv_impl::CompoundMultiplier<!DoTranspose, DoTranspose>::rowMultiplyPrecompose(
      *this, row, rhs, res, op);
}

template <bool ColWise, typename MatrixT1, typename MatrixT2>
template <bool DoTranspose, typename RhsT, typename ResT, typename PostOp>
void CompoundBlockMatrix<ColWise, MatrixT1, MatrixT2>::colMultiplyPostcompose(
    const Index col, const RhsT& rhs, ResT& res, const PostOp& op) const {
  mv_impl::CompoundMultiplier<!DoTranspose,
                              DoTranspose>::colMultiplyPostcompose(*this, col,
                                                                   rhs, res,
                                                                   op);
}

// Colwise = false
// Row-wise compound matrix

template <typename MatrixT1, typename MatrixT2>
CompoundBlockMatrix<false, MatrixT1, MatrixT2>::CompoundBlockMatrix(
    const IterableBlockObject<MatrixT1>& first,
    const IterableBlockObject<MatrixT2>& second)
    : m_first(first), m_second(second) {
  // Macro offsets
  m_compoundOffsets[0] = 0;
  m_compoundOffsets[1] = m_first.rows();
  m_compoundOffsets[2] = rows();

  // Underlying Offsets
  m_offsets.resize(rowsOfBlocks() + 1);
  std::copy(m_first.rowOffsets(),
            m_first.rowOffsets() + m_first.rowsOfBlocks() + 1,
            m_offsets.begin());

  const Index secStart = m_first.rowsOfBlocks() + 1;
  for (Index row = 0; row < m_second.rowsOfBlocks(); ++row) {
    m_offsets[secStart + row] =
        m_second.rowOffsets()[row + 1] + m_compoundOffsets[1];
  }
}

template <typename MatrixT1, typename MatrixT2>
template <bool DoTranspose, typename RhsT, typename ResT>
void CompoundBlockMatrix<false, MatrixT1, MatrixT2>::multiply(
    const RhsT& rhs, ResT& res, Scalar alpha, Scalar beta) const {
  mv_impl::CompoundMultiplier<DoTranspose, DoTranspose>::multiply(
      *this, rhs, res, alpha, beta);
}

template <typename MatrixT1, typename MatrixT2>
template <bool DoTranspose, typename RhsT, typename ResT, typename PreOp>
void CompoundBlockMatrix<false, MatrixT1, MatrixT2>::rowMultiplyPrecompose(
    const Index row, const RhsT& rhs, ResT& res, const PreOp& op) const {
  mv_impl::CompoundMultiplier<DoTranspose, DoTranspose>::rowMultiplyPrecompose(
      *this, row, rhs, res, op);
}

template <typename MatrixT1, typename MatrixT2>
template <bool DoTranspose, typename RhsT, typename ResT, typename PostOp>
void CompoundBlockMatrix<false, MatrixT1, MatrixT2>::colMultiplyPostcompose(
    const Index col, const RhsT& rhs, ResT& res, const PostOp& op) const {
  mv_impl::CompoundMultiplier<DoTranspose, DoTranspose>::colMultiplyPostcompose(
      *this, col, rhs, res, op);
}

}  // namespace bogus

#endif
