/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_ACCESS_HPP
#define BOGUS_BLOCK_ACCESS_HPP

#include "../Utils/CppTools.hpp"
#include "Constants.hpp"
#include "Traits.hpp"

namespace bogus {

//! Specialization of transpose_block() for self-adjoint types
template <typename SelfTransposeT>
inline typename EnableIf<BlockTraits<SelfTransposeT>::is_self_transpose,
                         const SelfTransposeT&>::ReturnType
transpose_block(const SelfTransposeT& block) {
  return block;
}

//! Specialization of transpose_block() for types that define a
//! ConstTransposeReturnType
template <typename BlockType>
inline
    typename EnableIf<!BlockTraits<BlockType>::is_self_transpose,
                      typename BlockType::ConstTransposeReturnType>::ReturnType
    transpose_block(const BlockType& block) {
  return block.transpose();
}

//! Defines the transpose type of a \p BlockType using self-introspection
/*! Process af follow:
                - If the BlockType is self transpose, the transpose type is the
   BlockType itself
                - If the BlockType defines a ConstTransposeReturnType, use it
                - If BlockTransposeTraits< BlockType > defines a ReturnType, use
   it
                - If the BlockType defines a Base, retry with this Base
  */
template <typename BlockType,
          bool IsSelfTranspose = BlockTraits<BlockType>::is_self_transpose,
          bool DefinesConstTranspose =
              HasConstTransposeReturnType<BlockType>::Value,
          bool DefinesTransposeTraits =
              HasReturnType<BlockTransposeTraits<BlockType> >::Value>
struct BlockTranspose {
  enum { is_defined = 0 };
};

// Self-transpose
template <typename BlockType, bool DCT, bool DTT>
struct BlockTranspose<BlockType, true, DCT, DTT> {
  typedef const BlockType& ReturnType;
  enum { is_defined = 1 };
};
// ConstTransposeReturnType
template <typename BlockType, bool DTT>
struct BlockTranspose<BlockType, false, true, DTT> {
  typedef typename BlockType::ConstTransposeReturnType ReturnType;
  enum { is_defined = 1 };
};
// BlockTransposeTraits
template <typename BlockType>
struct BlockTranspose<BlockType, false, false, true> {
  typedef typename BlockTransposeTraits<BlockType>::ReturnType ReturnType;
  enum { is_defined = 1 };
};

template <typename BlockType>
struct IsTransposable {
  enum { Value = BlockTranspose<BlockType>::is_defined };
};

//! Utility struct for expressing a compile-time conditional transpose of a
//! block
template <bool DoTranspose>
struct TransposeIf {
  template <typename BlockT>
  inline static const BlockT& get(const BlockT& src) {
    return src;
  }
};
template <>
struct TransposeIf<true> {
  template <typename BlockT>
  inline static typename BlockTranspose<BlockT>::ReturnType get(
      const BlockT& src) {
    return transpose_block(src);
  }
};

//! Utility struct to handle both compile-time and runtime optionally transposed
//! ops
// Case with runtime test
template <bool RuntimeTest, bool>
struct BlockTransposeOption {
  enum { runtime_test = 1, runtime_transpose = 0 };

  template <typename RhsT, typename ResT, typename Scalar>
  inline static void assign(const RhsT& rhs, ResT& res, const Scalar scale,
                            bool doTranspose) {
    if (doTranspose)
      res = scale * transpose_block(rhs);
    else
      res = scale * rhs;
  }
};

// Case with compile-time check
template <bool CompileTimeTranspose>
struct BlockTransposeOption<false, CompileTimeTranspose> {
  enum { runtime_test = 0, runtime_transpose = 1 };

  template <typename RhsT, typename ResT, typename Scalar>
  inline static void assign(const RhsT& rhs, ResT& res, const Scalar scale,
                            bool) {
    res = scale * TransposeIf<CompileTimeTranspose>::get(rhs);
  }
};

//! Access to the dimensions of a block
template <typename BlockT, bool Transpose_ = false>
struct BlockDims {
  typedef BlockTraits<BlockT> Traits;
  typedef SwapIf<Transpose_, Traits::RowsAtCompileTime,
                 Traits::ColsAtCompileTime>
      Dims;

  enum { Rows = Dims::First, Cols = Dims::Second };
};

namespace internal {

template <int dimension, typename VectorType>
struct SegmenterTraits {
  typedef
      typename VectorType::template NRowsBlockXpr<dimension>::Type ReturnType;
  typedef typename VectorType::template ConstNRowsBlockXpr<dimension>::Type
      ConstReturnType;
};

template <int dimension, typename VectorType>
struct SegmenterTraits<dimension, const VectorType> {
  typedef typename VectorType::template ConstNRowsBlockXpr<dimension>::Type
      ConstReturnType;
  typedef ConstReturnType ReturnType;
};

template <typename VectorType>
struct SegmenterTraits<DYNAMIC, VectorType> {
  typedef typename VectorType::RowsBlockXpr ReturnType;
  typedef typename VectorType::ConstRowsBlockXpr ConstReturnType;
};

template <typename VectorType>
struct SegmenterTraits<DYNAMIC, const VectorType> {
  typedef typename VectorType::ConstRowsBlockXpr ConstReturnType;
  typedef ConstReturnType ReturnType;
};
}  // namespace internal

//! Access to segment of a vector corresponding to a given block-row
template <int DimensionAtCompileTime, typename VectorType, typename Index>
struct Segmenter {
  enum { dimension = DimensionAtCompileTime };

  typedef internal::SegmenterTraits<dimension, VectorType> Traits;
  typedef typename Traits::ReturnType ReturnType;
  typedef typename Traits::ConstReturnType ConstReturnType;

  Segmenter(VectorType& vec, const Index*) : m_vec(vec) {}

  inline ReturnType operator[](const Index inner) {
    return m_vec.template middleRows<dimension>(dimension * inner);
  }

  inline ConstReturnType operator[](const Index inner) const {
    return m_vec.template middleRows<dimension>(dimension * inner);
  }

 private:
  VectorType& m_vec;
};

template <typename VectorType, typename Index>
struct Segmenter<internal::DYNAMIC, VectorType, Index> {
  typedef internal::SegmenterTraits<internal::DYNAMIC, VectorType> Traits;
  typedef typename Traits::ReturnType ReturnType;
  typedef typename Traits::ConstReturnType ConstReturnType;

  Segmenter(VectorType& vec, const Index* offsets)
      : m_vec(vec), m_offsets(offsets) {}

  inline ReturnType operator[](const Index inner) {
    return m_vec.middleRows(m_offsets[inner],
                            m_offsets[inner + 1] - m_offsets[inner]);
  }

  inline ConstReturnType operator[](const Index inner) const {
    return m_vec.middleRows(m_offsets[inner],
                            m_offsets[inner + 1] - m_offsets[inner]);
  }

 private:
  VectorType& m_vec;
  const Index* m_offsets;
};

}  // namespace bogus

#endif
