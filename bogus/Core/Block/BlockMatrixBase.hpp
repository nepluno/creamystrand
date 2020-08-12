/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCKMATRIX_HPP
#define BOGUS_BLOCKMATRIX_HPP

#include "IterableBlockObject.hpp"

namespace bogus {

//! Base class for dense and sparse block matrices, thought dense don't exist
//! yet
template <typename Derived>
class BlockMatrixBase : public IterableBlockObject<Derived> {
 public:
  typedef BlockMatrixTraits<Derived> Traits;
  typedef typename Traits::Index Index;
  typedef typename Traits::Scalar Scalar;

  // Supplemental Traits interface for BlockMatrixBase
  typedef typename Traits::BlockType BlockType;
  typedef typename Traits::BlockRef BlockRef;
  typedef typename Traits::ConstBlockRef ConstBlockRef;
  typedef typename Traits::BlockPtr BlockPtr;
  typedef typename Traits::BlocksArrayType BlocksArrayType;

  typedef IterableBlockObject<Derived> Base;
  using Base::derived;

  //! Return value of blockPtr( Index, Index ) for non-existing block
  static const BlockPtr InvalidBlockPtr;

  BlockMatrixBase() : m_rows(0), m_cols(0) {}

  virtual ~BlockMatrixBase() {}

  //! Multiplies a given block-row of the matrix with \p rhs, omitting the
  //! diagonal block
  /*! I.e. res = [ M( row, 0 ) ... M( row, row-1 ) 0 M( row, row+1 ) ... M( row,
     colsOfBlocks()-1 ) ] * rhs

          Especially useful inside a Gauss-Seidel algorithm.
          \warning Does not work on sparse, non-symmetric column major matrices
  */
  template <typename RhsT, typename ResT>
  void splitRowMultiply(const Index row, const RhsT& rhs, ResT& res) const {
    derived().splitRowMultiply(row, rhs, res);
  }

  //! Return a BlockPtr to the block a (row, col) or InvalidBlockPtr if it does
  //! not exist
  BlockPtr blockPtr(Index row, Index col) const {
    return derived().blockPtr(row, col);
  }
  //! Return a BlockPtr to the block a (row, row) or InvalidBlockPtr if it does
  //! not exist
  BlockPtr diagonalBlockPtr(Index row) const {
    return derived().diagonalBlockPtr(row);
  }

  //! Returns a reference to a block using the result from blockPtr() or
  //! diagonalBlockPtr()
  ConstBlockRef block(BlockPtr ptr) const { return derived().block(ptr); }
  //! Returns a reference to a block using the result from blockPtr() or
  //! diagonalBlockPtr()
  BlockRef block(BlockPtr ptr) { return derived().block(ptr); }

  //! Returns the total number of rows of the matrix ( expanding blocks )
  Index rows() const { return m_rows; }
  //! Returns the total number of columns of the matrix ( expanding blocks )
  Index cols() const { return m_cols; }

  //! Returns the total number of blocks of the matrix
  Index size() const;

  //! Access to blocks data
  const typename Traits::BlocksArrayType& blocks() const { return m_blocks; }
  //! Access to blocks data as a raw pointer
  const BlockType* data() const { return data_pointer(m_blocks); }
  //! Access to blocks data as a raw pointer
  BlockType* data() { return data_pointer(m_blocks); }

  //! \warning block has to exist
  BlockRef diagonal(const Index row) { return block(diagonalBlockPtr(row)); }
  //! \warning block has to exist
  ConstBlockRef diagonal(const Index row) const {
    return block(diagonalBlockPtr(row));
  }

  //! \warning block has to exist
  BlockRef block(Index row, Index col) { return block(blockPtr(row, col)); }
  //! \warning block has to exist
  ConstBlockRef block(Index row, Index col) const {
    return block(blockPtr(row, col));
  }

  //! Iterates over each block of a given row. Calls func( col, block )
  template <typename Func>
  void eachBlockOfRow(const Index row, Func func) const {
    derived().template eachBlockOf<false, Func>(row, func);
  }

  //! Iterates over each block of a given col. Calls func( row, block )
  template <typename Func>
  void eachBlockOfCol(const Index col, Func func) const {
    derived().template eachBlockOf<true, Func>(col, func);
  }

  //! Compile-time block properties
  enum CompileTimeProperties {
    RowsPerBlock = Traits::RowsPerBlock,
    ColsPerBlock = Traits::ColsPerBlock,
    has_row_major_blocks = BlockTraits<BlockType>::is_row_major,
    has_square_or_dynamic_blocks = ColsPerBlock == RowsPerBlock,
    has_fixed_rows_blocks = ((int)RowsPerBlock != internal::DYNAMIC),
    has_fixed_cols_blocks = ((int)ColsPerBlock != internal::DYNAMIC),
    has_fixed_size_blocks = has_fixed_cols_blocks && has_fixed_rows_blocks
  };

 protected:
  Index m_rows;
  Index m_cols;

  BlocksArrayType m_blocks;
};

}  // namespace bogus

#endif  // BLOCKMATRIX_HPP
