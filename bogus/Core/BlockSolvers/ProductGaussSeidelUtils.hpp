
/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_PRODUCT_GAUSS_SEIDEL_UTILS_HPP
#define BOGUS_PRODUCT_GAUSS_SEIDEL_UTILS_HPP

#include <vector>

#include "../Utils/CppTools.hpp"

namespace bogus {

namespace block_solvers_impl {

//! Wrapper for the block-diagonal matrix of the ProductGaussSeidel solver
/*! Specializations allows to:
 *   - Use a constant times the Identity matrix by default without requiring the
 * user to call any function or specify template parameters
 *   - Specify a pointer to an arbitrary block-diagonal matrix when the second
 * template parameter of ProductGaussSeidel is explicitely set, and precompute
 *     the indices of its diagonal blocks
 */
template <typename Type, bool IsScalar>
struct DiagonalMatrixWrapper {
  typedef Type Scalar;
  ConstantArray<Scalar> array;
  DiagonalMatrixWrapper(const Scalar& s = 1) : array(s) {}

  bool valid() const { return true; }
  Scalar get() const { return array[0]; }
  const ConstantArray<Scalar>& asArray() const { return array; }
};
template <typename Type>
struct DiagonalMatrixWrapper<Type, false> {
  typedef typename Type::Index Index;
  typedef typename Type::BlockPtr BlockPtr;

  DiagonalMatrixWrapper() : m_matrixPtr(BOGUS_NULL_PTR(const Type)) {}
  DiagonalMatrixWrapper(const Type& diag) : m_matrixPtr(&diag) {
    computeBlockIndices();
  }

  bool valid() const { return m_matrixPtr; }
  const Type& get() const { return *m_matrixPtr; }
  const DiagonalMatrixWrapper& asArray() const { return *this; }

  void computeBlockIndices();

  inline bool has_element(const Index i) const {
    return m_blockIndices[i] != Type::InvalidBlockPtr;
  }

  // \warning has_block(i) must be true
  inline const typename Type::BlockType& operator[](const Index i) const {
    return m_matrixPtr->block(m_blockIndices[i]);
  }

 private:
  const Type* m_matrixPtr;
  std::vector<BlockPtr> m_blockIndices;
};

//! Utility struct to precompute or reference the (D M') part of the product
template <typename MType, typename DType, bool Precompute>
struct DMtStorage {
  DMtStorage()
      : m_M(BOGUS_NULL_PTR(const MType)), m_D(BOGUS_NULL_PTR(const DType)) {}

  void compute(const MType& M, const DType& D) {
    m_M = &M;
    m_D = &D;
  }

  template <typename Rhs, typename Intermediate, typename Res>
  void multiply(const Rhs& rhs, Intermediate& itm, Res& res) const;

  template <typename Rhs, typename Res>
  void colMultiply(typename MType::Index col, const Rhs& rhs, Res& res) const;
  template <typename Rhs, typename Res>
  void rowMultiply(typename MType::Index row, const Rhs& rhs, Res& res) const;

  const MType* m_M;
  const DType* m_D;
};

template <typename MType, typename DType>
struct DMtStorage<MType, DType, true> {
  DMtStorage() : m_M(BOGUS_NULL_PTR(const MType)) {}

  void compute(const MType& M, const DType& D);

  template <typename Rhs, typename Intermediate, typename Res>
  void multiply(const Rhs& rhs, Intermediate& itm, Res& res) const;

  template <typename Rhs, typename Res>
  void colMultiply(typename MType::Index col, const Rhs& rhs, Res& res) const;
  template <typename Rhs, typename Res>
  void rowMultiply(typename MType::Index row, const Rhs& rhs, Res& res) const;

 private:
  const MType* m_M;
  typedef
      typename MType::template MutableImpl<typename MType::TransposeBlockType,
                                           false, true>::Type DMtType;
  DMtType m_DMt;
};

}  // namespace block_solvers_impl

}  // namespace bogus

#endif
