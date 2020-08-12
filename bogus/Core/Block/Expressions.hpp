/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_EXPRESSIONS_HPP
#define BOGUS_BLOCK_EXPRESSIONS_HPP

#include "BlockObjectBase.hpp"

namespace bogus {

//! Base class for Transpose views of a BlockObjectBase
template <typename MatrixT>
struct Transpose : public BlockObjectBase<Transpose<MatrixT> > {
  typedef BlockObjectBase<Transpose<MatrixT> > Base;
  typedef BlockMatrixTraits<Transpose<MatrixT> > Traits;
  typedef typename Traits::PlainObjectType PlainObjectType;
  typedef typename Traits::Index Index;
  typedef typename Base::Scalar Scalar;

  const PlainObjectType& matrix;

  Transpose(const PlainObjectType& m) : matrix(m.derived()) {}

  typename Base::ConstTransposeReturnType transpose() const { return matrix; }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT& rhs, ResT& res, Scalar alpha = 1,
                Scalar beta = 0) const {
    matrix.template multiply<!DoTranspose>(rhs, res, alpha, beta);
  }

  Index rows() const { return matrix.cols(); }
  Index cols() const { return matrix.rows(); }

  Index rowsOfBlocks() const { return matrix.colsOfBlocks(); }
  Index colsOfBlocks() const { return matrix.rowsOfBlocks(); }
  Index blockRows(Index row) const { return matrix.blockCols(row); }
  Index blockCols(Index col) const { return matrix.blockRows(col); }
  const Index* rowOffsets() const { return matrix.colOffsets(); }
  const Index* colOffsets() const { return matrix.rowOffsets(); }
};

template <typename MatrixT>
struct BlockMatrixTraits<Transpose<MatrixT> > {
  typedef BlockMatrixTraits<MatrixT> OrigTraits;
  typedef typename OrigTraits::Index Index;
  typedef typename OrigTraits::Scalar Scalar;

  enum {
    is_transposed = 1,
    is_temporary = 1,
    is_symmetric = OrigTraits::is_symmetric
  };
  enum {
    RowsPerBlock = OrigTraits::ColsPerBlock,
    ColsPerBlock = OrigTraits::RowsPerBlock
  };

  typedef typename OrigTraits::PlainObjectType PlainObjectType;
  typedef const PlainObjectType& ConstTransposeReturnType;
  typedef PlainObjectType TransposeObjectType;
};

template <typename ObjectT, bool IsTemporary>
struct BlockStorage {
  typedef const ObjectT& ConstValue;
};
template <typename ObjectT>
struct BlockStorage<ObjectT, true> {
  typedef const ObjectT ConstValue;
};

template <typename ObjectT>
struct BlockOperand {
  typedef ObjectT ObjectType;
  typedef typename ObjectT::PlainObjectType PlainObjectType;

  typedef BlockMatrixTraits<ObjectT> Traits;
  enum { do_transpose = Traits::is_transposed };
  typedef typename Traits::Scalar Scalar;

  typename BlockStorage<ObjectT, Traits::is_temporary>::ConstValue object;
  Scalar scaling;

  BlockOperand(const ObjectT& o, Scalar s = 1) : object(o), scaling(s) {}
};

template <template <typename, typename> class BlockOp, typename LhsMatrixT,
          typename RhsMatrixT>
struct BinaryBlockOp
    : public BlockObjectBase<BlockOp<LhsMatrixT, RhsMatrixT> > {
  typedef BlockObjectBase<BlockOp<LhsMatrixT, RhsMatrixT> > Base;
  typedef typename Base::PlainObjectType PlainObjectType;

  typedef BlockOperand<LhsMatrixT> Lhs;
  typedef BlockOperand<RhsMatrixT> Rhs;

  typedef typename Lhs::PlainObjectType PlainLhsMatrixType;
  typedef typename Rhs::PlainObjectType PlainRhsMatrixType;

  const Lhs lhs;
  const Rhs rhs;
  enum { transposeLhs = Lhs::do_transpose, transposeRhs = Rhs::do_transpose };

  BinaryBlockOp(const LhsMatrixT& l, const RhsMatrixT& r,
                typename Lhs::Scalar lscaling = 1,
                typename Lhs::Scalar rscaling = 1)
      : lhs(l, lscaling), rhs(r, rscaling) {}
};

template <typename LhsMatrixT, typename RhsMatrixT>
struct Product : public BinaryBlockOp<Product, LhsMatrixT, RhsMatrixT> {
  typedef BinaryBlockOp<bogus::Product, LhsMatrixT, RhsMatrixT> Base;
  typedef typename Base::Scalar Scalar;
  typedef typename Base::Index Index;

  Product(const LhsMatrixT& l, const RhsMatrixT& r,
          typename Base::Lhs::Scalar lscaling = 1,
          typename Base::Lhs::Scalar rscaling = 1)
      : Base(l, r, lscaling, rscaling) {}

  typename Base::ConstTransposeReturnType transpose() const {
    return typename Base::ConstTransposeReturnType(
        Base::rhs.object.transpose(), Base::lhs.object.transpose(),
        Base::rhs.scaling, Base::lhs.scaling);
  }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT& rhs, ResT& res, Scalar alpha = 1,
                Scalar beta = 0) const;

  Index rows() const { return Base::lhs.object.rows(); }
  Index cols() const { return Base::rhs.object.cols(); }

  Index rowsOfBlocks() const { return Base::lhs.object.rowsOfBlocks(); }
  Index colsOfBlocks() const { return Base::rhs.object.colsOfBlocks(); }
  Index blockRows(Index row) const { return Base::lhs.object.blockRows(row); }
  Index blockCols(Index col) const { return Base::rhs.object.blockCols(col); }
  const Index* rowOffsets() const { return Base::lhs.object.rowOffsets(); }
  const Index* colOffsets() const { return Base::rhs.object.colOffsets(); }
};

template <typename LhsMatrixT, typename RhsMatrixT>
struct BlockMatrixTraits<Product<LhsMatrixT, RhsMatrixT> > {
  typedef BlockMatrixTraits<LhsMatrixT> LhsTraits;
  typedef BlockMatrixTraits<RhsMatrixT> RhsTraits;

  typedef typename LhsTraits::Index Index;
  typedef typename LhsTraits::Scalar Scalar;

  enum { is_transposed = 0, is_temporary = 1, is_symmetric = 0 };

  typedef Product<LhsMatrixT, RhsMatrixT> ProductType;

  typedef
      typename BlockMatrixTraits<typename LhsTraits::PlainObjectType>::BlockType
          LhsBlockType;
  typedef
      typename BlockMatrixTraits<typename RhsTraits::PlainObjectType>::BlockType
          RhsBlockType;

  typedef typename BlockBlockProductTraits<
      LhsBlockType, RhsBlockType, LhsTraits::is_transposed,
      RhsTraits::is_transposed>::ReturnType ResBlockType;

  typedef typename LhsTraits::PlainObjectType ::template MutableImpl<
      ResBlockType, false>::Type PlainObjectType;

  enum {
    RowsPerBlock = LhsTraits::RowsPerBlock,
    ColsPerBlock = RhsTraits::ColsPerBlock
  };

  typedef Product<
      typename BlockOperand<RhsMatrixT>::ObjectType::TransposeObjectType,
      typename BlockOperand<LhsMatrixT>::ObjectType::TransposeObjectType>
      ConstTransposeReturnType;
  typedef ConstTransposeReturnType TransposeObjectType;
};

template <typename LhsMatrixT, typename RhsMatrixT>
struct Addition : public BinaryBlockOp<Addition, LhsMatrixT, RhsMatrixT> {
  typedef BinaryBlockOp<bogus::Addition, LhsMatrixT, RhsMatrixT> Base;
  typedef typename Base::Scalar Scalar;
  typedef typename Base::Index Index;

  Addition(const LhsMatrixT& l, const RhsMatrixT& r,
           typename Base::Lhs::Scalar lscaling = 1,
           typename Base::Lhs::Scalar rscaling = 1)
      : Base(l, r, lscaling, rscaling) {}

  typename Base::ConstTransposeReturnType transpose() const {
    return typename Base::ConstTransposeReturnType(
        Base::lhs.object.transpose(), Base::rhs.object.transpose(),
        Base::lhs.scaling, Base::rhs.scaling);
  }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT& rhs, ResT& res, Scalar alpha = 1,
                Scalar beta = 0) const {
    Base::lhs.object.template multiply<DoTranspose>(
        rhs, res, alpha * Base::lhs.scaling, beta);
    Base::rhs.object.template multiply<DoTranspose>(
        rhs, res, alpha * Base::rhs.scaling, 1);
  }

  Index rows() const { return Base::lhs.object.rows(); }
  Index cols() const { return Base::lhs.object.cols(); }

  Index rowsOfBlocks() const { return Base::lhs.object.rowsOfBlocks(); }
  Index colsOfBlocks() const { return Base::lhs.object.colsOfBlocks(); }
  Index blockRows(Index row) const { return Base::lhs.object.blockRows(row); }
  Index blockCols(Index col) const { return Base::lhs.object.blockCols(col); }
  const Index* rowOffsets() const { return Base::lhs.object.rowOffsets(); }
  const Index* colOffsets() const { return Base::lhs.object.colOffsets(); }
};

template <typename LhsMatrixT, typename RhsMatrixT>
struct BlockMatrixTraits<Addition<LhsMatrixT, RhsMatrixT> > {
  typedef BlockMatrixTraits<LhsMatrixT> OrigTraits;
  typedef typename OrigTraits::Index Index;
  typedef typename OrigTraits::Scalar Scalar;

  typedef typename BlockMatrixTraits<
      typename OrigTraits::PlainObjectType>::BlockType ResBlockType;

  typedef typename OrigTraits::PlainObjectType ::template MutableImpl<
      ResBlockType, false>::Type PlainObjectType;

  enum {
    is_transposed = 0,
    is_temporary = 1,
    is_symmetric = (BlockMatrixTraits<LhsMatrixT>::is_symmetric &&
                    BlockMatrixTraits<RhsMatrixT>::is_symmetric)
  };
  enum {
    RowsPerBlock = OrigTraits::RowsPerBlock,
    ColsPerBlock = OrigTraits::ColsPerBlock
  };

  typedef Addition<
      typename BlockOperand<LhsMatrixT>::ObjectType::TransposeObjectType,
      typename BlockOperand<RhsMatrixT>::ObjectType::TransposeObjectType>
      ConstTransposeReturnType;
  typedef ConstTransposeReturnType TransposeObjectType;
};

template <typename MatrixT>
struct Scaling : public BlockObjectBase<Scaling<MatrixT> > {
  typedef BlockOperand<MatrixT> Operand;
  typedef typename Operand::PlainObjectType PlainOperandMatrixType;
  Operand operand;

  enum { transposeOperand = Operand::do_transpose };

  typedef BlockObjectBase<Scaling> Base;

  typedef typename Base::Scalar Scalar;
  typedef typename Base::Index Index;

  Scaling(const MatrixT& object, const typename MatrixT::Scalar scaling)
      : operand(object, scaling) {}

  typename Base::ConstTransposeReturnType transpose() const {
    return typename Base::ConstTransposeReturnType(operand.object.transpose(),
                                                   operand.scaling);
  }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT& rhs, ResT& res, Scalar alpha = 1,
                Scalar beta = 0) const {
    operand.object.template multiply<DoTranspose>(
        rhs, res, alpha * operand.scaling, beta);
  }

  typename Base::Index rows() const { return operand.object.rows(); }
  typename Base::Index cols() const { return operand.object.cols(); }

  Index rowsOfBlocks() const { return operand.object.rowsOfBlocks(); }
  Index colsOfBlocks() const { return operand.object.colsOfBlocks(); }
  Index blockRows(Index row) const { return operand.object.blockRows(row); }
  Index blockCols(Index col) const { return operand.object.blockCols(col); }
  const Index* rowOffsets() const { return operand.object.rowOffsets(); }
  const Index* colOffsets() const { return operand.object.colOffsets(); }
};

template <typename MatrixT>
struct BlockMatrixTraits<Scaling<MatrixT> > {
  typedef BlockMatrixTraits<MatrixT> OrigTraits;
  typedef typename OrigTraits::Index Index;
  typedef typename OrigTraits::Scalar Scalar;

  enum {
    is_symmetric = OrigTraits::is_symmetric,
    is_transposed = 0,
    is_temporary = 1
  };
  enum {
    RowsPerBlock = OrigTraits::RowsPerBlock,
    ColsPerBlock = OrigTraits::ColsPerBlock
  };

  typedef typename OrigTraits::PlainObjectType PlainObjectType;

  typedef Scaling<
      typename BlockOperand<MatrixT>::ObjectType::TransposeObjectType>
      ConstTransposeReturnType;
  typedef ConstTransposeReturnType TransposeObjectType;
};

template <typename ObjectT>
struct BlockOperand<Scaling<ObjectT> > : public BlockOperand<ObjectT> {
  typedef BlockOperand<ObjectT> Base;

  BlockOperand(const Scaling<ObjectT>& o, typename Base::Scalar s = 1)
      : Base(o.operand.object, s * o.operand.scaling) {}
  BlockOperand(const typename Base::ObjectType& o, typename Base::Scalar s = 1)
      : Base(o, s) {}
};

}  // namespace bogus

#endif  // EXPRESSIONS_HPP
