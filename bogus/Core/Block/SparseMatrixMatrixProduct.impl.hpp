/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_SPARSE_MATRIXMATRIX_PRODUCT_IMPL_HPP
#define BOGUS_SPARSE_MATRIXMATRIX_PRODUCT_IMPL_HPP

#include "Access.hpp"
#include "Expressions.hpp"
#include "SparseBlockIndexComputer.hpp"
#include "SparseBlockMatrixBase.hpp"

#ifndef BOGUS_DONT_PARALLELIZE
#include <omp.h>
#endif

namespace bogus {

template <typename Derived>
template <typename LhsT, typename RhsT>
Derived &SparseBlockMatrixBase<Derived>::operator=(
    const Product<LhsT, RhsT> &prod) {
  // Default to column-major startegy
  // TODO use row-major strategy without inner indices comparison for  ( H H^T
  // )-type products
  setFromProduct<true>(prod);
  return derived();
}

namespace mm_impl {

// A block-block product
template <typename Index_, typename BlockPtr>
struct SparseBlockProductTerm {
  typedef Index_ Index;
  Index index;  // Inner index in product matrix

  BlockPtr lhsPtr;
  BlockPtr rhsPtr;
  bool lhsIsAfterDiag;
  bool rhsIsAfterDiag;

  SparseBlockProductTerm(Index idx, BlockPtr lPtr, BlockPtr rPtr, bool lIAD,
                         bool rIAD)
      : index(idx),
        lhsPtr(lPtr),
        rhsPtr(rPtr),
        lhsIsAfterDiag(lIAD),
        rhsIsAfterDiag(rIAD) {}

  // Ordering so all terms for a single final block are consecutive
  bool operator<(const SparseBlockProductTerm &rhs) const {
    return index < rhs.index;
  }
};

// Block structure computation using row-major startegy (deprecated)
template <bool ColWise, typename Index, typename BlockPtr, bool is_symmetric,
          bool is_col_major>
struct SparseBlockProductIndex {
  typedef SparseBlockProductTerm<Index, BlockPtr> Term;

  typedef std::vector<Term> InnerType;
  std::vector<InnerType> to_compute;

  typedef SparseBlockIndex<true, Index, BlockPtr> CompressedIndexType;
  CompressedIndexType compressed;

  template <typename LhsIndex, typename RhsIndex>
  void compute(const Index outerSize, const Index innerSize,
               const LhsIndex &lhsIdx, const RhsIndex &rhsIdx) {
    assert(lhsIdx.innerSize() == rhsIdx.innerSize());

    to_compute.resize(outerSize);

#ifndef BOGUS_DONT_PARALLELIZE
    SparseBlockIndex<false, Index, BlockPtr> uncompressed;
    uncompressed.resizeOuter(outerSize);
#else
    compressed.resizeOuter(outerSize);
#endif

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
    for (Index i = 0; i < outerSize; ++i) {
      const Index last = is_symmetric ? i + 1 : innerSize;
      for (Index j = 0; j != last; ++j) {
        bool nonZero = false;

        const Index lhsOuter = is_col_major ? j : i;
        const Index rhsOuter = is_col_major ? i : j;
        typename LhsIndex::InnerIterator lhs_it(lhsIdx, lhsOuter);
        typename RhsIndex::InnerIterator rhs_it(rhsIdx, rhsOuter);

        while (lhs_it && rhs_it) {
          if (lhs_it.inner() > rhs_it.inner())
            ++rhs_it;
          else if (lhs_it.inner() < rhs_it.inner())
            ++lhs_it;
          else {
            to_compute[i].push_back(Term(j, lhs_it.ptr(), rhs_it.ptr(),
                                         lhs_it.after(lhsOuter),
                                         rhs_it.after(rhsOuter)));
            nonZero = true;
            ++lhs_it;
            ++rhs_it;
          }
        }

        if (nonZero) {
#ifndef BOGUS_DONT_PARALLELIZE
          uncompressed.insertBack(i, j, 0);
#else
          compressed.insertBack(i, j, 0);
#endif
        }
      }
    }

#ifndef BOGUS_DONT_PARALLELIZE
    compressed = uncompressed;
#else
    compressed.finalize();
#endif
  }
};

// Block structure computation using col-major strategy
// (Using less conditionals, as it does not have to compare inner indices in the
// tighter loop)
template <typename Index, typename BlockPtr, bool is_symmetric,
          bool is_col_major>
struct SparseBlockProductIndex<true, Index, BlockPtr, is_symmetric,
                               is_col_major> {
  typedef SparseBlockProductTerm<Index, BlockPtr> Term;

  typedef std::vector<Term> InnerType;
  std::vector<InnerType> to_compute;

  typedef SparseBlockIndex<true, Index, BlockPtr> CompressedIndexType;
  CompressedIndexType compressed;

  template <typename LhsIndex, typename RhsIndex>
  void compute(const Index outerSize, const Index innerSize,
               const LhsIndex &lhsIdx, const RhsIndex &rhsIdx) {
    assert(lhsIdx.outerSize() == rhsIdx.outerSize());
    (void)innerSize;

    const Index productSize = lhsIdx.outerSize();

    to_compute.resize(outerSize);
    compressed.resizeOuter(outerSize);

#ifdef BOGUS_DONT_PARALLELIZE
    std::vector<InnerType> &loc_compute = to_compute;
#else

    std::vector<std::vector<InnerType> > temp_compute;

#pragma omp parallel
    {

#pragma omp master
      { temp_compute.resize(omp_get_num_threads()); }
#pragma omp barrier

      std::vector<InnerType> &loc_compute = temp_compute[omp_get_thread_num()];
      loc_compute.resize(outerSize);

#pragma omp for
#endif
    for (Index i = 0; i < productSize; ++i) {
      if (is_col_major) {
        for (typename RhsIndex::InnerIterator rhs_it(rhsIdx, i); rhs_it;
             ++rhs_it) {
          for (typename LhsIndex::InnerIterator lhs_it(lhsIdx, i);
               lhs_it && (!is_symmetric || lhs_it.inner() <= rhs_it.inner());
               ++lhs_it) {
            loc_compute[rhs_it.inner()].push_back(
                Term((Index)lhs_it.inner(), lhs_it.ptr(), rhs_it.ptr(),
                     lhs_it.after(i), rhs_it.after(i)));
          }
        }
      } else {
        for (typename LhsIndex::InnerIterator lhs_it(lhsIdx, i); lhs_it;
             ++lhs_it) {
          for (typename RhsIndex::InnerIterator rhs_it(rhsIdx, i);
               rhs_it && (!is_symmetric || rhs_it.inner() <= lhs_it.inner());
               ++rhs_it) {
            loc_compute[lhs_it.inner()].push_back(
                Term((Index)rhs_it.inner(), lhs_it.ptr(), rhs_it.ptr(),
                     lhs_it.after(i), rhs_it.after(i)));
          }
        }
      }
    }

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp for
    for (Index i = 0; i < outerSize; ++i) {
      for (int t = 0; t < (int)temp_compute.size(); ++t) {
        const std::vector<InnerType> &loc_compute = temp_compute[t];
        to_compute[i].insert(to_compute[i].end(), loc_compute[i].begin(),
                             loc_compute[i].end());
      }
      std::sort(to_compute[i].begin(), to_compute[i].end());
    }
  }
#endif

  for (Index i = 0; i < outerSize; ++i) {
#ifdef BOGUS_DONT_PARALLELIZE
    std::sort(to_compute[i].begin(), to_compute[i].end());
#endif
    Index prevIndex = -1;
    for (std::size_t j = 0; j != to_compute[i].size(); ++j) {
      if (to_compute[i][j].index != prevIndex) {
        prevIndex = to_compute[i][j].index;
        compressed.insertBack(i, prevIndex, 0);
      }
    }
  }

  compressed.finalize();
}

};  // namespace mm_impl

template <bool LhsRuntimeTest, bool RhsRunTimeTest,
          bool LhsCompileTimeTranspose, bool RhsCompileTimeTranspose>
struct BinaryTransposeOption {
  template <typename LhsT, typename RhsT, typename ResT>
  static inline void mm_assign(const LhsT &lhs, const RhsT &rhs, ResT &res,
                               bool transposeLhs, bool transposeRhs) {
    if (transposeLhs)
      if (transposeRhs)
        res = transpose_block(lhs) * transpose_block(rhs);
      else
        res = transpose_block(lhs) * rhs;
    else if (transposeRhs)
      res = lhs * transpose_block(rhs);
    else
      res = lhs * rhs;
  }

  template <typename LhsT, typename RhsT, typename ResT>
  static inline void mm_add(const LhsT &lhs, const RhsT &rhs, ResT &res,
                            bool transposeLhs, bool transposeRhs) {
    if (transposeLhs)
      if (transposeRhs)
        res += transpose_block(lhs) * transpose_block(rhs);
      else
        res += transpose_block(lhs) * rhs;
    else if (transposeRhs)
      res += lhs * transpose_block(rhs);
    else
      res += lhs * rhs;
  }
};

// Lhs transpose known at compile time
template <bool LhsTranspose, bool RhsTranspose>
struct BinaryTransposeOption<false, true, LhsTranspose, RhsTranspose> {
  template <typename LhsT, typename RhsT, typename ResT>
  static inline void mm_assign(const LhsT &lhs, const RhsT &rhs, ResT &res,
                               bool, bool transposeRhs) {
    if (transposeRhs)
      res = TransposeIf<LhsTranspose>::get(lhs) * transpose_block(rhs);
    else
      res = TransposeIf<LhsTranspose>::get(lhs) * rhs;
  }
  template <typename LhsT, typename RhsT, typename ResT>
  static inline void mm_add(const LhsT &lhs, const RhsT &rhs, ResT &res, bool,
                            bool transposeRhs) {
    if (transposeRhs)
      res += TransposeIf<LhsTranspose>::get(lhs) * transpose_block(rhs);
    else
      res += TransposeIf<LhsTranspose>::get(lhs) * rhs;
  }
};

// Rhs transpose known at compile time
template <bool LhsTranspose, bool RhsTranspose>
struct BinaryTransposeOption<true, false, LhsTranspose, RhsTranspose> {
  template <typename LhsT, typename RhsT, typename ResT>
  static inline void mm_assign(const LhsT &lhs, const RhsT &rhs, ResT &res,
                               bool transposeLhs, bool) {
    if (transposeLhs)
      res = transpose_block(lhs) * TransposeIf<RhsTranspose>::get(rhs);
    else
      res = lhs * TransposeIf<RhsTranspose>::get(rhs);
  }
  template <typename LhsT, typename RhsT, typename ResT>
  static inline void mm_add(const LhsT &lhs, const RhsT &rhs, ResT &res,
                            bool transposeLhs, bool) {
    if (transposeLhs)
      res += transpose_block(lhs) * TransposeIf<RhsTranspose>::get(rhs);
    else
      res += lhs * TransposeIf<RhsTranspose>::get(rhs);
  }
};

// Both transpose known at compile time
template <bool LhsTranspose, bool RhsTranspose>
struct BinaryTransposeOption<false, false, LhsTranspose, RhsTranspose> {
  template <typename LhsT, typename RhsT, typename ResT>
  static inline void mm_assign(const LhsT &lhs, const RhsT &rhs, ResT &res,
                               bool, bool) {
    res = TransposeIf<LhsTranspose>::get(lhs) *
          TransposeIf<RhsTranspose>::get(rhs);
  }

  template <typename LhsT, typename RhsT, typename ResT>
  static inline void mm_add(const LhsT &lhs, const RhsT &rhs, ResT &res, bool,
                            bool) {
    res = res + TransposeIf<LhsTranspose>::get(lhs) *
                    TransposeIf<RhsTranspose>::get(rhs);
  }
};

template <typename TransposeOption, typename ResBlockRef, typename ProductIndex,
          typename LhsBlocks, typename RhsBlocks, typename ResBlocks,
          typename Scalar>
static void compute_blocks(const ProductIndex &productIndex,
                           const std::size_t nBlocks,
                           const LhsBlocks &lhsBlocks,
                           const RhsBlocks &rhsBlocks, ResBlocks &resBlocks,
                           Scalar scaling) {
  typedef typename ProductIndex::Term Term;
  typedef typename Term::Index Index;
  typedef std::pair<const Term *, unsigned> BlockComputation;

  std::vector<BlockComputation> flat_compute(nBlocks);

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
  for (Index i = 0; i < (Index)productIndex.to_compute.size(); ++i) {
    const typename ProductIndex::InnerType &innerVec =
        productIndex.to_compute[i];

    const Index n = static_cast<Index>(innerVec.size());
    if (n != 0) {
      typename ProductIndex::CompressedIndexType::InnerIterator c_it(
          productIndex.compressed, i);

      for (Index j = 0, curStart = 0;; ++j) {
        if (j == n || innerVec[j].index != c_it.inner()) {
          flat_compute[c_it.ptr()].first = &innerVec[curStart];
          flat_compute[c_it.ptr()].second = j - curStart;

          ++c_it;
          if (!c_it) break;

          curStart = j;
        }
      }

      assert(!c_it);
    }
  }

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
  for (std::ptrdiff_t i = 0; i < (std::ptrdiff_t)nBlocks; ++i) {
    ResBlockRef b = resBlocks[i];

    const Term *bc = flat_compute[i].first;
    const unsigned n = flat_compute[i].second;

    TransposeOption::mm_assign(lhsBlocks[bc->lhsPtr], rhsBlocks[bc->rhsPtr], b,
                               bc->lhsIsAfterDiag, bc->rhsIsAfterDiag);

    for (unsigned j = 1; j != n; ++j) {
      ++bc;
      TransposeOption::mm_add(lhsBlocks[bc->lhsPtr], rhsBlocks[bc->rhsPtr], b,
                              bc->lhsIsAfterDiag, bc->rhsIsAfterDiag);
    }

    if (scaling != 1) b = b * scaling;
  }
}

}  // namespace bogus

template <typename Derived>
template <bool ColWise, typename LhsT, typename RhsT>
void SparseBlockMatrixBase<Derived>::setFromProduct(
    const Product<LhsT, RhsT> &prod) {
  typedef Product<LhsT, RhsT> Prod;

  Evaluator<typename Prod::Lhs::ObjectType> lhs(prod.lhs.object);
  Evaluator<typename Prod::Rhs::ObjectType> rhs(prod.rhs.object);
  typedef BlockMatrixTraits<typename Prod::PlainLhsMatrixType> LhsTraits;
  typedef BlockMatrixTraits<typename Prod::PlainRhsMatrixType> RhsTraits;

  BOGUS_STATIC_ASSERT(!Prod::transposeLhs ||
                          IsTransposable<typename LhsTraits::BlockType>::Value,
                      TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE);
  BOGUS_STATIC_ASSERT(!Prod::transposeRhs ||
                          IsTransposable<typename RhsTraits::BlockType>::Value,
                      TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE);

  clear();
  if (Prod::transposeLhs) {
    m_rows = lhs->cols();
    colMajorIndex().innerOffsets = lhs->rowMajorIndex().innerOffsets;
  } else {
    m_rows = lhs->rows();
    colMajorIndex().innerOffsets = lhs->colMajorIndex().innerOffsets;
  }
  if (Prod::transposeRhs) {
    m_cols = rhs->rows();
    rowMajorIndex().innerOffsets = rhs->colMajorIndex().innerOffsets;
  } else {
    m_cols = rhs->cols();
    rowMajorIndex().innerOffsets = rhs->rowMajorIndex().innerOffsets;
  }

  rowMajorIndex().resizeOuter(colMajorIndex().innerSize());
  colMajorIndex().resizeOuter(rowMajorIndex().innerSize());

  {
    typedef mm_impl::SparseBlockProductIndex<
        ColWise, Index, BlockPtr, Traits::is_symmetric, Traits::is_col_major>
        ProductIndex;
    ProductIndex productIndex;

    {
      SparseBlockIndexComputer<typename Prod::PlainLhsMatrixType, ColWise,
                               Prod::transposeLhs>
          lhsIndexComputer(*lhs);
      SparseBlockIndexComputer<typename Prod::PlainRhsMatrixType, !ColWise,
                               Prod::transposeRhs>
          rhsIndexComputer(*rhs);

      productIndex.compute(majorIndex().outerSize(), minorIndex().outerSize(),
                           lhsIndexComputer.get(), rhsIndexComputer.get());
    }

    const unsigned outerSize = majorIndex().outerSize();
    createBlockShapes(productIndex.compressed.outer[outerSize],
                      productIndex.compressed, m_blocks);

    typedef mm_impl::BinaryTransposeOption<
        LhsTraits::is_symmetric &&
            !(BlockTraits<typename LhsTraits::BlockType>::is_self_transpose),
        RhsTraits::is_symmetric &&
            !(BlockTraits<typename RhsTraits::BlockType>::is_self_transpose),
        Prod::transposeLhs, Prod::transposeRhs>
        TransposeOption;

    mm_impl::template compute_blocks<TransposeOption, BlockRef>(
        productIndex, nBlocks(), lhs->blocks(), rhs->blocks(), m_blocks,
        prod.lhs.scaling * prod.rhs.scaling);

    productIndex.compressed.valid = true;
    m_majorIndex.move(productIndex.compressed);
    m_minorIndex.valid = false;
  }

  assert(m_majorIndex.valid);

  Finalizer::finalize(*this);
}

}  // namespace mm_impl

#endif  // SPARSEPRODUCT_IMPL_HPP
