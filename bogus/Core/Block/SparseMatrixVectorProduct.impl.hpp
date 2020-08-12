/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_SPARSE_MATRIXVECTOR_PRODUCT_HPP
#define BOGUS_SPARSE_MATRIXVECTOR_PRODUCT_HPP

#include "Access.hpp"
#include "Expressions.hpp"
#include "SparseBlockIndexComputer.hpp"
#include "SparseBlockMatrixBase.hpp"

namespace bogus {

namespace mv_impl {

template <typename ResT, typename Scalar>
void scale_vector(ResT& res, const Scalar beta) {
  if ((Scalar)0 == beta) res.setZero();
  if ((Scalar)1 != beta) res *= beta;
}

template <bool DoTranspose, typename Matrix, typename RhsT, typename ResT,
          typename Scalar>
inline void mv_add_pre(const Matrix& matrix, const RhsT& rhs, ResT& res,
                       Scalar alpha) {
  res.noalias() += TransposeIf<DoTranspose>::get(matrix) * (alpha * rhs);
}

template <bool DoTranspose, typename Matrix, typename RhsT, typename ResT,
          typename Scalar>
inline void mv_add_post(const Matrix& matrix, const RhsT& rhs, ResT& res,
                        Scalar alpha) {
  res.noalias() += alpha * (TransposeIf<DoTranspose>::get(matrix) * rhs);
}

template <bool Transpose, typename BlockType, typename BlocksT, typename IndexT,
          typename RhsT, typename ResT, typename ScalarT>
static inline void innerRowMultiply(const BlocksT& blocks, const IndexT& index,
                                    const typename IndexT::Index outerIdx,
                                    const RhsT& rhs, ResT& res, ScalarT alpha) {
  const Segmenter<BlockDims<BlockType, Transpose>::Cols, const RhsT,
                  typename IndexT::Index>
      segmenter(rhs, index.innerOffsetsData());

  for (typename IndexT::InnerIterator it(index, outerIdx); it; ++it) {
    if (alpha.has_element(it.inner())) {
      mv_add_pre<Transpose>(blocks[it.ptr()], segmenter[it.inner()], res,
                            alpha[it.inner()]);
    }
  }
}

template <bool Transpose, typename BlockType, typename BlocksT, typename IndexT,
          typename RhsT, typename ResT, typename ScalarT>
static inline void innerColMultiply(const BlocksT& blocks, const IndexT& index,
                                    const typename IndexT::Index outerIdx,
                                    const RhsT& rhs, ResT& res, ScalarT alpha) {
  typedef Segmenter<BlockDims<BlockType, Transpose>::Rows, ResT,
                    typename IndexT::Index>
      ResSegmenter;
  ResSegmenter segmenter(res, index.innerOffsetsData());

  for (typename IndexT::InnerIterator it(index, outerIdx); it; ++it) {
    typename ResSegmenter::ReturnType res_seg(segmenter[it.inner()]);
    if (alpha.has_element(it.inner())) {
      mv_add_post<Transpose>(blocks[it.ptr()], rhs, res_seg, alpha[it.inner()]);
    }
  }
}

//! Implementation for non-symmetric, in order matrix/vector products
template <bool Symmetric, bool NativeOrder, bool Transpose>
struct SparseBlockMatrixVectorMultiplier {
  template <typename Derived, typename RhsT, typename ResT, typename ScalarT>
  static void multiply(const SparseBlockMatrixBase<Derived>& matrix,
                       const RhsT& rhs, ResT& res, ScalarT alpha) {
    const int ResSegDim =
        BlockDims<typename Derived::BlockType, Transpose>::Rows;
    typedef Segmenter<ResSegDim, ResT, typename Derived::Index> ResSegmenter;
    ResSegmenter resSegmenter(res, matrix.minorIndex().innerOffsetsData());

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
    for (typename Derived::Index i = 0;
         i < (typename Derived::Index)matrix.majorIndex().outerSize(); ++i) {
      typename ResSegmenter::ReturnType seg = resSegmenter[i];
      innerRowMultiply<Transpose, typename Derived::BlockType>(
          matrix.blocks(), matrix.majorIndex(), i, rhs, seg,
          make_constant_array(alpha));
    }
  }
};

//! Implementation for symmetric products
template <bool NativeOrder, bool Transpose>
struct SparseBlockMatrixVectorMultiplier<true, NativeOrder, Transpose> {
#ifndef BOGUS_DONT_PARALLELIZE
  template <typename Derived, typename RhsT, typename ResT, typename LocalResT,
            typename ScalarT>
  static void multiplyAndReduct(const SparseBlockMatrixBase<Derived>& matrix,
                                const RhsT& rhs, ResT& res, const LocalResT&,
                                ScalarT alpha) {
    const int SegDim = BlockDims<typename Derived::BlockType, Transpose>::Rows;

    typedef Segmenter<SegDim, LocalResT, typename Derived::Index> ResSegmenter;
    typedef Segmenter<SegDim, const RhsT, typename Derived::Index> RhsSegmenter;
    const RhsSegmenter rhsSegmenter(rhs,
                                    matrix.majorIndex().innerOffsetsData());

    const Lock& lock = matrix.lock();

    typedef
        typename SparseBlockMatrixBase<Derived>::MajorIndexType MajorIndexType;
#pragma omp parallel
    {
      LocalResT locRes(res.rows(), res.cols());
      locRes.setZero();

      ResSegmenter resSegmenter(locRes, matrix.minorIndex().innerOffsetsData());

#pragma omp for
      for (typename Derived::Index i = 0; i < matrix.majorIndex().outerSize();
           ++i) {
        typename RhsSegmenter::ConstReturnType rhs_seg(rhsSegmenter[i]);
        typename ResSegmenter::ReturnType res_seg(resSegmenter[i]);
        for (typename MajorIndexType::InnerIterator it(matrix.majorIndex(), i);
             it; ++it) {
          const typename Derived::BlockType& b = matrix.block(it.ptr());
          mv_add_pre<false>(b, rhsSegmenter[it.inner()], res_seg, alpha);
          if (it.inner() != i) {
            typename ResSegmenter::ReturnType inner_res_seg(
                resSegmenter[it.inner()]);
            mv_add_pre<true>(b, rhs_seg, inner_res_seg, alpha);
          }
        }
      }

      {
        Lock::Guard<> guard(lock);
        res += locRes;
      }
    }
  }
#endif

  template <typename Derived, typename RhsT, typename ResT, typename ScalarT>
  static void multiply(const SparseBlockMatrixBase<Derived>& matrix,
                       const RhsT& rhs, ResT& res, ScalarT alpha) {
    const int SegDim = BlockDims<typename Derived::BlockType, Transpose>::Rows;

    typedef Segmenter<SegDim, ResT, typename Derived::Index> ResSegmenter;
    ResSegmenter resSegmenter(res, matrix.minorIndex().innerOffsetsData());
    typedef Segmenter<SegDim, const RhsT, typename Derived::Index> RhsSegmenter;
    const RhsSegmenter rhsSegmenter(rhs,
                                    matrix.majorIndex().innerOffsetsData());

    typedef typename Derived::Index Index;
    if (matrix.transposeIndex().valid) {
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
      for (Index i = 0; i < matrix.majorIndex().outerSize(); ++i) {
        typename ResSegmenter::ReturnType seg(resSegmenter[i]);
        innerRowMultiply<false, typename Derived::BlockType>(
            matrix.blocks(), matrix.majorIndex(), i, rhs, seg,
            make_constant_array(alpha));
        innerRowMultiply<false, typename Derived::TransposeBlockType>(
            matrix.transposeBlocks(), matrix.transposeIndex(), i, rhs, seg,
            make_constant_array(alpha));
      }
    } else {
#ifdef BOGUS_DONT_PARALLELIZE
      for (Index i = 0; i < matrix.majorIndex().outerSize(); ++i) {
        typename RhsSegmenter::ConstReturnType rhs_seg(rhsSegmenter[i]);
        typename ResSegmenter::ReturnType res_seg(resSegmenter[i]);
        for (typename SparseBlockMatrixBase<
                 Derived>::MajorIndexType::InnerIterator it(matrix.majorIndex(),
                                                            i);
             it; ++it) {
          const typename Derived::BlockType& b = matrix.block(it.ptr());
          mv_add_pre<false>(b, rhsSegmenter[it.inner()], res_seg, alpha);
          if (it.inner() != i) {
            typename ResSegmenter::ReturnType inner_res_seg(
                resSegmenter[it.inner()]);
            mv_add_pre<true>(b, rhs_seg, inner_res_seg, alpha);
          }
        }
      }
#else
      multiplyAndReduct(matrix, rhs, res, get_mutable_vector(res), alpha);
#endif
    }
  }
};

//! Implemntation for non-symmetric, out-of-order products
template <bool Transpose>
struct OutOfOrderSparseBlockMatrixVectorMultiplier {
#ifndef BOGUS_DONT_PARALLELIZE
  template <typename Derived, typename RhsT, typename ResT, typename LocalResT,
            typename ScalarT>
  static void multiplyAndReduct(const SparseBlockMatrixBase<Derived>& matrix,
                                const RhsT& rhs, ResT& res, const LocalResT&,
                                ScalarT alpha) {
    typedef typename Derived::Index Index;

    const int RhsSegDim =
        BlockDims<typename Derived::BlockType, Transpose>::Cols;
    typedef Segmenter<RhsSegDim, const RhsT, Index> RhsSegmenter;
    const RhsSegmenter rhsSegmenter(rhs,
                                    matrix.minorIndex().innerOffsetsData());

    const Lock& lock = matrix.lock();

#pragma omp parallel
    {
      LocalResT locRes(res.rows(), res.cols());
      locRes.setZero();

#pragma omp for
      for (Index i = 0; i < matrix.majorIndex().outerSize(); ++i) {
        innerColMultiply<Transpose, typename Derived::BlockType>(
            matrix.blocks(), matrix.majorIndex(), i, rhsSegmenter[i], locRes,
            make_constant_array(alpha));
      }

      {
        Lock::Guard<> guard(lock);
        res += locRes;
      }
    }
  }
#endif

  template <typename Derived, typename RhsT, typename ResT, typename ScalarT>
  static void multiply(const SparseBlockMatrixBase<Derived>& matrix,
                       const RhsT& rhs, ResT& res, ScalarT alpha) {
    typedef typename SparseBlockMatrixBase<Derived>::Index Index;

    const int RhsSegDim =
        BlockDims<typename Derived::BlockType, Transpose>::Cols;
    typedef Segmenter<RhsSegDim, const RhsT, Index> RhsSegmenter;
    const RhsSegmenter rhsSegmenter(rhs,
                                    matrix.minorIndex().innerOffsetsData());

#ifdef BOGUS_DONT_PARALLELIZE
    for (Index i = 0; i < matrix.majorIndex().outerSize(); ++i) {
      innerColMultiply<Transpose, typename Derived::BlockType>(
          matrix.blocks(), matrix.majorIndex(), i, rhsSegmenter[i], res,
          make_constant_array(alpha));
    }
#else
    multiplyAndReduct(matrix, rhs, res, get_mutable_vector(res), alpha);
#endif
  }
};

//! Implemntation for non-symmetric, out-of-order products
template <bool Transpose>
struct SparseBlockMatrixVectorMultiplier<false, false, Transpose>
    : public OutOfOrderSparseBlockMatrixVectorMultiplier<Transpose> {};

//! Implemntation for non-symmetric, out-of-order transposed products using
//! cached transpose
template <>
struct SparseBlockMatrixVectorMultiplier<false, false, true>
    : public OutOfOrderSparseBlockMatrixVectorMultiplier<true> {
  typedef OutOfOrderSparseBlockMatrixVectorMultiplier<true> Base;

  template <typename Derived, typename RhsT, typename ResT, typename ScalarT>
  static void multiply(const SparseBlockMatrixBase<Derived>& matrix,
                       const RhsT& rhs, ResT& res, ScalarT alpha) {
    if (matrix.transposeIndex().valid) {
      const int ResSegDim = BlockDims<typename Derived::BlockType, true>::Rows;
      typedef Segmenter<ResSegDim, ResT, typename Derived::Index> ResSegmenter;
      ResSegmenter resSegmenter(res, matrix.majorIndex().innerOffsetsData());

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
      for (typename Derived::Index i = 0;
           i < matrix.transposeIndex().outerSize(); ++i) {
        typename ResSegmenter::ReturnType seg(resSegmenter[i]);
        innerRowMultiply<false, typename Derived::TransposeBlockType>(
            matrix.transposeBlocks(), matrix.transposeIndex(), i, rhs, seg,
            make_constant_array(alpha));
      }
    } else {
      Base::multiply(matrix, rhs, res, alpha);
    }
  }
};

// Split-row-multiply

//! Implementation for for non-symmetric matrices
template <bool Symmetric, bool NativeOrder>
struct SparseBlockRowMultiplier {
  template <typename Derived, typename RhsT, typename ResT>
  static void splitRowMultiply(const SparseBlockMatrixBase<Derived>& matrix,
                               typename Derived::Index row, const RhsT& rhs,
                               ResT& res) {
    typedef SparseBlockIndexComputer<Derived, false, false> IndexCptr;
    typedef typename IndexCptr::ReturnType IndexType;

    const Segmenter<BlockDims<typename Derived::BlockType, false>::Cols,
                    const RhsT, typename Derived::Index>
        segmenter(rhs, matrix.colOffsets());

    IndexCptr cptr(matrix);

    for (typename IndexType::InnerIterator it(cptr.get(), row); it; ++it) {
      if (it.inner() != row)
        mv_add_pre<false>(matrix.block(it.ptr()), segmenter[it.inner()], res,
                          1);
    }
  }

  template <bool DoTranspose, typename Derived, typename RhsT, typename ResT,
            typename PreOp>
  static void rowMultiply(const SparseBlockMatrixBase<Derived>& matrix,
                          typename Derived::Index row, const RhsT& rhs,
                          ResT& res, const PreOp& preop) {
    SparseBlockIndexComputer<Derived, false, DoTranspose> cptr(matrix);
    innerRowMultiply<DoTranspose, typename Derived::BlockType>(
        matrix.blocks(), cptr.get(), row, rhs, res, preop);
  }

  template <bool DoTranspose, typename Derived, typename RhsT, typename ResT,
            typename PostOp>
  static void colMultiply(const SparseBlockMatrixBase<Derived>& matrix,
                          typename Derived::Index col, const RhsT& rhs,
                          ResT& res, const PostOp& postop) {
    SparseBlockIndexComputer<Derived, true, DoTranspose> cptr(matrix);
    innerColMultiply<DoTranspose, typename Derived::BlockType>(
        matrix.blocks(), cptr.get(), col, rhs, res, postop);
  }
};

//! Implementation for symmetric matrices
template <bool NativeOrder>
struct SparseBlockRowMultiplier<true, NativeOrder> {
  template <typename Derived, typename RhsT, typename ResT>
  static void splitRowMultiply(const SparseBlockMatrixBase<Derived>& matrix,
                               typename Derived::Index row, const RhsT& rhs,
                               ResT& res) {
    const Segmenter<BlockDims<typename Derived::BlockType, !NativeOrder>::Cols,
                    const RhsT, typename Derived::Index>
        segmenter(rhs, matrix.colOffsets());

    for (typename SparseBlockMatrixBase<Derived>::MajorIndexType::InnerIterator
             it(matrix.majorIndex(), row);
         it && it.inner() < row; ++it) {
      mv_add_pre<!NativeOrder>(matrix.block(it.ptr()), segmenter[it.inner()],
                               res, 1);
    }

    if (matrix.transposeIndex().valid) {
      innerRowMultiply<!NativeOrder, typename Derived::TransposeBlockType>(
          matrix.transposeBlocks(), matrix.transposeIndex(), row, rhs, res,
          make_constant_array(1));
    } else {
      assert(matrix.minorIndex().valid);
      innerRowMultiply<NativeOrder, typename Derived::BlockType>(
          matrix.blocks(), matrix.minorIndex(), row, rhs, res,
          make_constant_array(1));
    }
  }

  template <bool DoTranspose, typename Derived, typename RhsT, typename ResT,
            typename PreOp>
  static void rowMultiply(const SparseBlockMatrixBase<Derived>& matrix,
                          typename Derived::Index row, const RhsT& rhs,
                          ResT& res, const PreOp& preop) {
    innerRowMultiply<!NativeOrder, typename Derived::BlockType>(
        matrix.blocks(), matrix.majorIndex(), row, rhs, res, preop);

    if (matrix.transposeIndex().valid) {
      innerRowMultiply<!NativeOrder, typename Derived::TransposeBlockType>(
          matrix.transposeBlocks(), matrix.transposeIndex(), row, rhs, res,
          preop);
    } else {
      assert(matrix.minorIndex().valid);
      innerRowMultiply<NativeOrder, typename Derived::BlockType>(
          matrix.blocks(), matrix.minorIndex(), row, rhs, res, preop);
    }
  }

  template <bool DoTranspose, typename Derived, typename RhsT, typename ResT,
            typename PostOp>
  static void colMultiply(const SparseBlockMatrixBase<Derived>& matrix,
                          typename Derived::Index col, const RhsT& rhs,
                          ResT& res, const PostOp& postop) {
    innerColMultiply<!NativeOrder, typename Derived::BlockType>(
        matrix.blocks(), matrix.majorIndex(), col, rhs, res, postop);

    if (matrix.transposeIndex().valid) {
      innerColMultiply<!NativeOrder, typename Derived::TransposeBlockType>(
          matrix.transposeBlocks(), matrix.transposeIndex(), col, rhs, res,
          postop);
    } else {
      assert(matrix.minorIndex().valid);
      innerColMultiply<NativeOrder, typename Derived::BlockType>(
          matrix.blocks(), matrix.minorIndex(), col, rhs, res, postop);
    }
  }
};

template <bool Transposed>
struct ProductVector {
  template <typename LhsDerived, typename RhsDerived, typename VecT,
            typename ResT, typename LocalResT, typename ScalarT>
  static void multiply(const BlockObjectBase<LhsDerived>& prodLhs,
                       const BlockObjectBase<RhsDerived>& prodRhs,
                       const VecT& rhs, ResT& res, const LocalResT&,
                       const ScalarT alpha, const ScalarT beta) {
    LocalResT buf(prodRhs.rows(), rhs.cols());
    prodRhs.template multiply<false>(rhs, buf, alpha, 0);
    prodLhs.template multiply<false>(buf, res, 1, beta);
  }
};
template <>
struct ProductVector<true> {
  template <typename LhsDerived, typename RhsDerived, typename VecT,
            typename ResT, typename LocalResT, typename ScalarT>
  static void multiply(const BlockObjectBase<LhsDerived>& prodLhs,
                       const BlockObjectBase<RhsDerived>& prodRhs,
                       const VecT& rhs, ResT& res, const LocalResT&,
                       const ScalarT alpha, const ScalarT beta) {
    LocalResT buf(prodLhs.cols(), rhs.cols());
    prodLhs.template multiply<true>(rhs, buf, alpha, 0);
    prodRhs.template multiply<true>(buf, res, 1, beta);
  }
};

}  // namespace mv_impl

// Proxy -- can be specialized to use external libraries

template <bool BSRFormat, bool BlocksAreRowMajor, typename Scalar,
          typename Index>
struct SparseBlockMatrixOpProxy {
  template <bool Transpose, typename Derived, typename RhsT, typename ResT>
  static void multiply(const SparseBlockMatrixBase<Derived>& matrix,
                       const RhsT& rhs, ResT& res, Scalar alpha, Scalar beta) {
    mv_impl::scale_vector(res, beta);

    typedef BlockMatrixTraits<Derived> Traits;
    mv_impl::SparseBlockMatrixVectorMultiplier<
        Traits::is_symmetric, Transpose == bool(Traits::is_col_major),
        Transpose>::multiply(matrix, rhs, res, alpha);
  }

  template <typename Derived, typename RhsT, typename ResT>
  static void splitRowMultiply(const SparseBlockMatrixBase<Derived>& matrix,
                               typename Derived::Index row, const RhsT& rhs,
                               ResT& res) {
    typedef BlockMatrixTraits<Derived> Traits;
    mv_impl::SparseBlockRowMultiplier<
        Traits::is_symmetric, !Traits::is_col_major>::splitRowMultiply(matrix,
                                                                       row, rhs,
                                                                       res);
  }

  template <bool DoTranspose, typename Derived, typename RhsT, typename ResT,
            typename PreOp>
  static void rowMultiply(const SparseBlockMatrixBase<Derived>& matrix,
                          typename Derived::Index row, const RhsT& rhs,
                          ResT& res, const PreOp& preop) {
    typedef BlockMatrixTraits<Derived> Traits;
    mv_impl::SparseBlockRowMultiplier<
        Traits::is_symmetric, DoTranspose == bool(Traits::is_col_major)>::
        template rowMultiply<DoTranspose>(matrix, row, rhs, res, preop);
  }

  template <bool DoTranspose, typename Derived, typename RhsT, typename ResT,
            typename PostOp>
  static void colMultiply(const SparseBlockMatrixBase<Derived>& matrix,
                          typename Derived::Index col, const RhsT& rhs,
                          ResT& res, const PostOp& postop) {
    typedef BlockMatrixTraits<Derived> Traits;
    mv_impl::SparseBlockRowMultiplier<
        Traits::is_symmetric, DoTranspose != bool(Traits::is_col_major)>::
        template colMultiply<DoTranspose>(matrix, col, rhs, res, postop);
  }
};

template <typename Derived>
template <bool Transpose, typename RhsT, typename ResT>
void SparseBlockMatrixBase<Derived>::multiply(
    const RhsT& rhs, ResT& res,
    typename SparseBlockMatrixBase<Derived>::Scalar alpha,
    typename SparseBlockMatrixBase<Derived>::Scalar beta) const {
  BOGUS_STATIC_ASSERT(!Transpose || IsTransposable<BlockType>::Value,
                      TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE);

  SparseBlockMatrixOpProxy<is_bsr_compatible, Base::has_row_major_blocks,
                           Scalar, Index>::template multiply<Transpose>(*this,
                                                                        rhs,
                                                                        res,
                                                                        alpha,
                                                                        beta);
}

template <typename Derived>
template <typename RhsT, typename ResT>
void SparseBlockMatrixBase<Derived>::splitRowMultiply(const Index row,
                                                      const RhsT& rhs,
                                                      ResT& res) const {
  SparseBlockMatrixOpProxy<is_bsr_compatible, Base::has_row_major_blocks,
                           Scalar, Index>::splitRowMultiply(*this, row, rhs,
                                                            res);
}

template <typename Derived>
template <bool Transpose, typename RhsT, typename ResT, typename PreOp>
void SparseBlockMatrixBase<Derived>::rowMultiplyPrecompose(
    const Index row, const RhsT& rhs, ResT& res, const PreOp& preop) const {
  SparseBlockMatrixOpProxy<is_bsr_compatible, Base::has_row_major_blocks,
                           Scalar,
                           Index>::template rowMultiply<Transpose>(*this, row,
                                                                   rhs, res,
                                                                   preop);
}

template <typename Derived>
template <bool Transpose, typename RhsT, typename ResT, typename PostOp>
void SparseBlockMatrixBase<Derived>::colMultiplyPostcompose(
    const Index col, const RhsT& rhs, ResT& res, const PostOp& postop) const {
  SparseBlockMatrixOpProxy<is_bsr_compatible, Base::has_row_major_blocks,
                           Scalar,
                           Index>::template colMultiply<Transpose>(*this, col,
                                                                   rhs, res,
                                                                   postop);
}

template <typename LhsMatrixT, typename RhsMatrixT>
template <bool DoTranspose, typename RhsT, typename ResT>
void Product<LhsMatrixT, RhsMatrixT>::multiply(const RhsT& rhs, ResT& res,
                                               Scalar alpha,
                                               Scalar beta) const {
  mv_impl::ProductVector<DoTranspose>::multiply(
      Base::lhs.object, Base::rhs.object, rhs, res, get_mutable_vector(res),
      alpha * Base::lhs.scaling * Base::rhs.scaling, beta);
}

template <typename Expression>
template <bool DoTranspose, typename RhsT, typename ResT>
void NarySum<Expression>::multiply(const RhsT& rhs, ResT& res, Scalar alpha,
                                   Scalar beta) const {
  mv_impl::scale_vector(res, beta);

  for (typename Sum::const_iterator it = members.begin(); it != members.end();
       ++it) {
    it->template multiply<DoTranspose>(rhs, res, alpha, 1);
  }
}

template <typename Scalar>
template <bool DoTranspose, typename RhsT, typename ResT>
void Zero<Scalar>::multiply(const RhsT&, ResT& res, Scalar, Scalar beta) const {
  mv_impl::scale_vector(res, beta);
}

}  // namespace bogus

#endif
