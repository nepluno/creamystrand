/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_ZERO_HPP
#define BOGUS_BLOCK_ZERO_HPP

#include "IterableBlockObject.hpp"

namespace bogus {

//! Representation of the null matrix.
/*! Useful as argument to functions expecting a BlockObjectBase */
template <typename Scalar>
class Zero : public IterableBlockObject<Zero<Scalar> > {
 public:
  typedef IterableBlockObject<Zero<Scalar> > Base;

  typedef typename Base::Index Index;
  typedef typename Base::ConstTransposeReturnType ConstTransposeReturnType;

  explicit Zero(Index rows = 0, Index cols = 0) : m_rows(rows), m_cols(cols) {
    m_rowOffsets[0] = 0;
    m_rowOffsets[1] = rows;
    m_colOffsets[0] = 0;
    m_colOffsets[1] = cols;
  }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  Index blockRows(Index) const { return rows(); }
  Index blockCols(Index) const { return cols(); }

  Index rowsOfBlocks() const { return 1; }
  Index colsOfBlocks() const { return 1; }

  const Index *rowOffsets() const { return &m_rowOffsets; }
  const Index *colOffsets() const { return &m_colOffsets; }

  ConstTransposeReturnType transpose() const { *this; }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT &, ResT &res, Scalar = 1, Scalar beta = 0) const;

  // IterableBlockObject

  Index size() const { return 0; }

  template <typename Func>
  void eachBlockOfRow(const Index, Func) {}
  template <typename Func>
  void eachBlockOfCol(const Index, Func) {}

  template <bool DoTranspose, typename RhsT, typename ResT, typename PreOp>
  void rowMultiplyPrecompose(const Index, const RhsT &, ResT &,
                             const PreOp &) const {}
  template <bool DoTranspose, typename RhsT, typename ResT, typename PostOp>
  void colMultiplyPostcompose(const Index, const RhsT &, ResT &,
                              const PostOp &) const {}

 private:
  const Index m_rows;
  const Index m_cols;
  Index m_rowOffsets[2];
  Index m_colOffsets[2];
};

template <typename Scalar_>
struct BlockMatrixTraits<Zero<Scalar_> >
    : public BlockMatrixTraits<BlockObjectBase<Zero<Scalar_> > > {
  typedef Scalar_ Scalar;

  enum {
    is_symmetric = 1,
  };

  typedef Zero<Scalar> PlainObjectType;
  typedef const PlainObjectType &ConstTransposeReturnType;
  typedef PlainObjectType TransposeObjectType;
};

}  // namespace bogus

#endif
