/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_ITERABLE_BLOCK_OBJECT
#define BOGUS_ITERABLE_BLOCK_OBJECT

#include "../Utils/CppTools.hpp"
#include "BlockObjectBase.hpp"

namespace bogus {

//! Base class for matrix-like objects that define a block structure, but not a
//! block type
template <typename Derived>
class IterableBlockObject : public BlockObjectBase<Derived> {
 public:
  typedef BlockMatrixTraits<Derived> Traits;
  typedef typename Traits::Index Index;
  typedef typename Traits::Scalar Scalar;

  typedef BlockObjectBase<Derived> Base;
  using Base::derived;

  //! Returns the total number of blocks of the matrix
  Index size() const { return derived().size(); }

  //! Multiplication with a single row
  template <bool DoTranspose, typename RhsT, typename ResT>
  void rowMultiply(const Index row, const RhsT& rhs, ResT& res) const {
    rowMultiplyPrecompose<DoTranspose>(row, rhs, res, make_constant_array(1));
  }
  template <bool DoTranspose, typename RhsT, typename ResT, typename PreOp>
  void rowMultiplyPrecompose(const Index row, const RhsT& rhs, ResT& res,
                             const PreOp& op) const {
    derived().template rowMultiplyPrecompose<DoTranspose>(row, rhs, res, op);
  }

  //! Multiplication with a single column
  template <bool DoTranspose, typename RhsT, typename ResT>
  void colMultiply(const Index col, const RhsT& rhs, ResT& res) const {
    colMultiplyPostcompose<DoTranspose>(col, rhs, res, make_constant_array(1));
  }
  template <bool DoTranspose, typename RhsT, typename ResT, typename PostOp>
  void colMultiplyPostcompose(const Index col, const RhsT& rhs, ResT& res,
                              const PostOp& op) const {
    derived().template colMultiplyPostcompose<DoTranspose>(col, rhs, res, op);
  }

  //! Iterates over each block of a given row. Calls func( col, block )
  template <typename Func>
  void eachBlockOfRow(const Index row, Func func) const {
    derived().template eachBlockOfRow<Func>(row, func);
  }

  //! Iterates over each block of a given col. Calls func( row, block )
  template <typename Func>
  void eachBlockOfCol(const Index col, Func func) const {
    derived().template eachBlockOfCol<Func>(col, func);
  }

  //! Should be overidden by mutable types
  template <typename OtherBlockType, bool PreserveSymmetry = true,
            bool SwitchDirection = false>
  struct MutableImpl {
    typedef typename Derived::PlainObjectType Type;
  };
};

}  // namespace bogus

#endif
