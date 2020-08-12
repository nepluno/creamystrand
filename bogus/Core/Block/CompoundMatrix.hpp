/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_COMPOUND_BLOCK_MATRIX_HPP
#define BOGUS_COMPOUND_BLOCK_MATRIX_HPP

#include "Access.hpp"
#include "IterableBlockObject.hpp"

namespace bogus {

//! A matrix made by concatenating two other matrices of possibly different
//! types
/*! \tparam ColWise If true, the two matrices are side-by-side (column
 * concatenation), otherwise they are on top of another
 *
 *  If ColWise, the two matrices must have the same number and size of rows of
 * blocks. Otherwise, they must have the same number and size of columns of
 * blocks.
 *
 *  Implements IterableBlockObject
 */
template <bool ColWise, typename MatrixT1, typename MatrixT2>
class CompoundBlockMatrix;

namespace internal {
// Inject a few traits for Compounds with similar blocks
template <bool ColWise, typename MatrixT1, typename MatrixT2>
struct CompoundBlockMatrixBase
    : public IterableBlockObject<
          CompoundBlockMatrix<ColWise, MatrixT1, MatrixT2> > {};

template <bool ColWise, typename MatrixT1>
struct CompoundBlockMatrixBase<ColWise, MatrixT1, MatrixT1>
    : public IterableBlockObject<
          CompoundBlockMatrix<ColWise, MatrixT1, MatrixT1> > {
  template <typename OtherBlockType, bool PreserveSymmetry = true,
            bool SwitchDirection = false>
  struct MutableImpl {
    typedef typename MatrixT1::template MutableImpl<
        OtherBlockType, PreserveSymmetry, SwitchDirection>::Type Type;
  };
};

template <typename MatrixT1, typename MatrixT2>
struct SameBlockMatrixTraits {};
template <typename MatrixT1>
struct SameBlockMatrixTraits<MatrixT1, MatrixT1> {
  typedef BlockMatrixTraits<MatrixT1> CommonTraits;
  typedef typename CommonTraits::BlockType BlockType;
};
}  // namespace internal

template <bool ColWise, typename MatrixT1, typename MatrixT2>
class CompoundBlockMatrix
    : public internal::CompoundBlockMatrixBase<ColWise, MatrixT1, MatrixT2> {
 public:
  // side-by-side (ColWise = true)
  typedef IterableBlockObject<CompoundBlockMatrix<ColWise, MatrixT1, MatrixT2> >
      Base;

  typedef typename Base::Index Index;
  typedef typename Base::Scalar Scalar;
  typedef typename Base::ConstTransposeReturnType ConstTransposeReturnType;

  CompoundBlockMatrix(const IterableBlockObject<MatrixT1>& first,
                      const IterableBlockObject<MatrixT2>& second);

  // BlockObjectBase

  Index rows() const { return m_first.rows(); }
  Index cols() const { return m_first.cols() + m_second.cols(); }

  Index blockRows(Index row) const { return m_first.blockRows(row); }
  Index blockCols(Index col) const {
    const Index off = col - secondBegin();
    return off < 0 ? m_first.blockCols(col) : m_second.blockCols(off);
  }

  Index rowsOfBlocks() const { return m_first.rowsOfBlocks(); }
  Index colsOfBlocks() const {
    return m_first.colsOfBlocks() + m_second.colsOfBlocks();
  }

  const Index* rowOffsets() const { return m_first.rowOffsets(); }
  const Index* colOffsets() const { return data_pointer(m_offsets); }

  ConstTransposeReturnType transpose() const {
    return Transpose<CompoundBlockMatrix>(*this);
  }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT& rhs, ResT& res, Scalar alpha = 1,
                Scalar beta = 0) const;

  // Iterable Block Object

  Index size() const { return m_first.size() + m_second.size(); }

  template <typename Func>
  void eachBlockOfRow(const Index row, Func func) const {
    m_first.derived().eachBlockOfRow(row, func);
    m_second.derived().eachBlockOfRow(row, func);
  }
  template <typename Func>
  void eachBlockOfCol(const Index col, Func func) const {
    const Index off = col - secondBegin();
    if (off < 0)
      m_first.derived().eachBlockOfCol(col, func);
    else
      m_second.derived().eachBlockOfCol(off, func);
  }

  template <bool DoTranspose, typename RhsT, typename ResT, typename PreOp>
  void rowMultiplyPrecompose(const Index row, const RhsT& rhs, ResT& res,
                             const PreOp& op) const;

  template <bool DoTranspose, typename RhsT, typename ResT, typename PostOp>
  void colMultiplyPostcompose(const Index col, const RhsT& rhs, ResT& res,
                              const PostOp& op) const;

  // Accessors

  const IterableBlockObject<MatrixT1>& first() const { return m_first; }
  const IterableBlockObject<MatrixT2>& second() const { return m_second; }
  const Index* compoundOffsets() const { return m_compoundOffsets; }
  Index secondBegin() const { return m_first.colsOfBlocks(); }

 private:
  const IterableBlockObject<MatrixT1>& m_first;
  const IterableBlockObject<MatrixT2>& m_second;
  Index m_compoundOffsets[3];
  std::vector<Index> m_offsets;
};

template <bool ColWise, typename MatrixT1, typename MatrixT2>
struct BlockMatrixTraits<CompoundBlockMatrix<ColWise, MatrixT1, MatrixT2> >
    : public BlockMatrixTraits<
          BlockObjectBase<CompoundBlockMatrix<ColWise, MatrixT1, MatrixT2> > >,
      internal::SameBlockMatrixTraits<MatrixT1, MatrixT2> {
  typedef BlockMatrixTraits<MatrixT1> OrigTraits;
  typedef BlockMatrixTraits<MatrixT2> OtherTraits;

  typedef typename OrigTraits::Scalar Scalar;
  typedef typename OrigTraits::Index Index;

  enum {
    RowsPerBlock = SwapIf < ColWise || ((int)OrigTraits::RowsPerBlock) ==
                                           (int)OtherTraits::RowsPerBlock,
    internal::DYNAMIC,
    OrigTraits::RowsPerBlock > ::First,
    ColsPerBlock = SwapIf < !ColWise || ((int)OrigTraits::ColsPerBlock) ==
                                            (int)OtherTraits::ColsPerBlock,
    internal::DYNAMIC,
    OrigTraits::ColsPerBlock > ::First
  };
};

template <typename MatrixT1, typename MatrixT2>
class CompoundBlockMatrix<false, MatrixT1, MatrixT2>
    : public internal::CompoundBlockMatrixBase<false, MatrixT1, MatrixT2> {
 public:
  // one-atop-the-other (ColWise = false)
  typedef IterableBlockObject<CompoundBlockMatrix<false, MatrixT1, MatrixT2> >
      Base;

  typedef typename Base::Index Index;
  typedef typename Base::Scalar Scalar;
  typedef typename Base::ConstTransposeReturnType ConstTransposeReturnType;

  CompoundBlockMatrix(const IterableBlockObject<MatrixT1>& first,
                      const IterableBlockObject<MatrixT2>& second);

  // BlockObjectBase

  Index cols() const { return m_first.cols(); }
  Index rows() const { return m_first.rows() + m_second.rows(); }

  Index blockCols(Index col) const { return m_first.blockCols(col); }
  Index blockRows(Index row) const {
    const Index off = row - secondBegin();
    return off < 0 ? m_first.blockRows(row) : m_second.blockRows(off);
  }

  Index colsOfBlocks() const { return m_first.colsOfBlocks(); }
  Index rowsOfBlocks() const {
    return m_first.rowsOfBlocks() + m_second.rowsOfBlocks();
  }

  const Index* colOffsets() const { return m_first.colOffsets(); }
  const Index* rowOffsets() const { return data_pointer(m_offsets); }

  ConstTransposeReturnType transpose() const {
    return Transpose<CompoundBlockMatrix>(*this);
  }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT& rhs, ResT& res, Scalar alpha = 1,
                Scalar beta = 0) const;

  // Iterable Block Object

  Index size() const { return m_first.size() + m_second.size(); }

  template <typename Func>
  void eachBlockOfCol(const Index col, Func func) const {
    m_first.derived().eachBlockOfRow(col, func);
    m_second.derived().eachBlockOfRow(col, func);
  }
  template <typename Func>
  void eachBlockOfRow(const Index row, Func func) const {
    const Index off = row - secondBegin();
    if (off < 0)
      m_first.derived().eachBlockOfCol(row, func);
    else
      m_second.derived().eachBlockOfCol(off, func);
  }

  template <bool DoTranspose, typename RhsT, typename ResT, typename PreOp>
  void rowMultiplyPrecompose(const Index row, const RhsT& rhs, ResT& res,
                             const PreOp& op) const;

  template <bool DoTranspose, typename RhsT, typename ResT, typename PostOp>
  void colMultiplyPostcompose(const Index col, const RhsT& rhs, ResT& res,
                              const PostOp& op) const;

  // Accessors

  const IterableBlockObject<MatrixT1>& first() const { return m_first; }
  const IterableBlockObject<MatrixT2>& second() const { return m_second; }
  const Index* compoundOffsets() const { return m_compoundOffsets; }
  Index secondBegin() const { return m_first.rowsOfBlocks(); }

 private:
  const IterableBlockObject<MatrixT1>& m_first;
  const IterableBlockObject<MatrixT2>& m_second;
  Index m_compoundOffsets[3];
  std::vector<Index> m_offsets;
};

}  // namespace bogus

#endif
