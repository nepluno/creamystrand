/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCKOBJECTBASE_HPP
#define BOGUS_BLOCKOBJECTBASE_HPP

#include "../Block.fwd.hpp"

namespace bogus {

//! Base class for anything block
template <typename Derived>
struct BlockObjectBase {
  //! Returns a const reference to the implementation
  const Derived& derived() const { return static_cast<const Derived&>(*this); }
  //! Returns a reference to the implementation
  Derived& derived() { return static_cast<Derived&>(*this); }

  typedef BlockMatrixTraits<Derived> Traits;

  typedef typename Traits::Index Index;
  typedef typename Traits::Scalar Scalar;
  typedef typename Traits::ConstTransposeReturnType ConstTransposeReturnType;
  typedef typename Traits::TransposeObjectType TransposeObjectType;
  enum { is_transposed = Traits::is_transposed };

  typedef typename Traits::PlainObjectType PlainObjectType;

  //! Returns the total number of rows of the matrix ( expanding blocks )
  Index rows() const { return derived().rows(); }
  //! Returns the total number of columns of the matrix ( expanding blocks )
  Index cols() const { return derived().cols(); }

  //! Returns the number of rows of a given block row
  Index blockRows(Index row) const { return derived().blockRows(row); }
  //! Returns the number of columns of a given block columns
  Index blockCols(Index col) const { return derived().blockCols(col); }

  //! Returns the number of block rows of the matrix
  Index rowsOfBlocks() const { return derived().rowsOfBlocks(); }
  //! Returns the number of block columns of the matrix
  Index colsOfBlocks() const { return derived().colsOfBlocks(); }

  //! Returns an array containing the first index of each row
  const Index* rowOffsets() const { return derived().rowOffsets(); }
  //! Returns an array containing the first index of each column
  const Index* colOffsets() const { return derived().colOffsets(); }

  //! Returns an array containing the first index of a given row
  Index rowOffset(Index row) const { return rowOffsets()[row]; }
  //! Returns an array containing the first index of a given columns
  Index colOffset(Index col) const { return colOffsets()[col]; }

  //! Return a const transposed view of this object
  ConstTransposeReturnType transpose() const { return derived().transpose(); }

  //! Performs a matrix vector multiplication
  /*! \tparam DoTranspose If true, performs \c res = \c alpha * \c M' * \c rhs +
     beta * res, otherwise \c res = \c alpha * M * \c rhs + beta * res
    */
  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT& rhs, ResT& res, Scalar alpha = 1,
                Scalar beta = 0) const {
    derived().template multiply<DoTranspose>(rhs, res, alpha, beta);
  }

  // To use as block type
  enum {
    RowsAtCompileTime = internal::DYNAMIC,
    ColsAtCompileTime = internal::DYNAMIC,
    is_self_transpose = Traits::is_symmetric
  };
};

//! Default specialization of traits for BlockMatrices
/*! Re-specialized for derived classes, see e.g. BlockMatrixTraits<
 * SparseBlockMatrix > */
template <typename Derived>
struct BlockMatrixTraits<BlockObjectBase<Derived> > {
  //! Index type -- for accessing elements, defining offsets, etc
  typedef BOGUS_DEFAULT_INDEX_TYPE Index;

  //! Type representing the transpose of this object
  typedef Transpose<Derived> TransposeObjectType;
  //! Type returned by the transpose() method
  /*! Generally a TransposeObjectType or a const reference to it */
  typedef TransposeObjectType ConstTransposeReturnType;

  enum {
    is_symmetric = 0,   //!< Whether the object is self-transpose
    is_transposed = 0,  //!< Whether this object represents a transposition
    is_temporary = 0    //!< Whether this object should be copied or ref'd to
  };

  enum {
    //! Number of rows per block at compile time, or DYNAMIC if unknown
    RowsPerBlock = internal::DYNAMIC,
    //! Number of cols per block at compile time, or DYNAMIC if unknown
    ColsPerBlock = internal::DYNAMIC
  };

  //! Type in which the expression may be evaluated into
  typedef Derived PlainObjectType;
};

}  // namespace bogus

#endif  // BLOCKMATRIX_HPP
