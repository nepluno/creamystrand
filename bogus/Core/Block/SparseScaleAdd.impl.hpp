/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_SPARSE_SCALEADD_HPP
#define BOGUS_SPARSE_SCALEADD_HPP

#include "Access.hpp"
#include "SparseBlockIndexComputer.hpp"
#include "SparseBlockMatrixBase.hpp"

namespace bogus {

template <typename Derived>
Derived &SparseBlockMatrixBase<Derived>::scale(const Scalar alpha) {
  if (alpha != 1) {
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
    for (std::ptrdiff_t i = 0; i < (std::ptrdiff_t)nBlocks(); ++i) {
      block(i) *= alpha;
    }
  }

  return derived();
}

template <typename Derived>
template <bool Transpose, typename OtherDerived>
Derived &SparseBlockMatrixBase<Derived>::add(
    const SparseBlockMatrixBase<OtherDerived> &rhs, Scalar alpha) {
  BOGUS_STATIC_ASSERT(
      !Transpose || IsTransposable<typename OtherDerived::BlockType>::Value,
      TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE);

  typedef typename SparseBlockMatrixBase<OtherDerived>::Traits OtherTraits;
  typedef std::pair<BlockPtr, typename OtherTraits::BlockPtr> PtrPair;
  typedef std::pair<Index, PtrPair> NonZero;

  if (rhs.empty()) return derived();

  std::vector<std::vector<NonZero> > nonZeros;

  // I - Compute non-zeros
  {
    SparseBlockIndexComputer<OtherDerived, Traits::is_col_major, Transpose>
        indexComputer(rhs);
    typedef
        typename SparseBlockIndexComputer<OtherDerived, Traits::is_col_major,
                                          Transpose>::ReturnType
            SourceIndexType;
    const SourceIndexType &rhsIndex = indexComputer.get();

    const MajorIndexType &lhsIndex = majorIndex();

    assert(rhsIndex.outerSize() == lhsIndex.outerSize());
    assert(rhsIndex.innerSize() == lhsIndex.innerSize());

    nonZeros.resize(lhsIndex.outerSize());

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
    for (Index i = 0; i < lhsIndex.outerSize(); ++i) {
      typename MajorIndexType::InnerIterator lhs_it(lhsIndex, i);
      typename SourceIndexType::InnerIterator rhs_it(rhsIndex, i);

      NonZero nz;
      while (lhs_it || rhs_it) {
        if (lhs_it && (!rhs_it || lhs_it.inner() < rhs_it.inner())) {
          nz.first = lhs_it.inner();
          nz.second.first = lhs_it.ptr();
          nz.second.second = OtherDerived::InvalidBlockPtr;
          ++lhs_it;
        } else if (rhs_it && (!lhs_it || rhs_it.inner() < lhs_it.inner())) {
          nz.first = rhs_it.inner();
          nz.second.first = InvalidBlockPtr;
          nz.second.second = rhs_it.ptr();
          ++rhs_it;
        } else {
          nz.first = lhs_it.inner();
          nz.second.first = lhs_it.ptr();
          nz.second.second = rhs_it.ptr();
          ++lhs_it;
          ++rhs_it;
        }

        if (Traits::is_symmetric && nz.first > i) break;

        nonZeros[i].push_back(nz);
      }
    }
  }

  std::vector<BlockPtr> offsets(nonZeros.size() + 1);

  MajorIndexType resIndex;
  resIndex.resizeOuter(nonZeros.size());
  offsets[0] = 0;
  for (unsigned i = 0; i < nonZeros.size(); ++i) {
    offsets[i + 1] = offsets[i] + nonZeros[i].size();
  }

  resIndex.reserve(offsets.back());
  for (unsigned i = 0; i < nonZeros.size(); ++i) {
    for (unsigned j = 0; j < nonZeros[i].size(); ++j) {
      const BlockPtr ptr = offsets[i] + j;
      const NonZero &nz = nonZeros[i][j];
      resIndex.insertBack(i, nz.first, ptr);
    }
  }
  resIndex.finalize();

  typename Traits::BlocksArrayType resBlocks;
  createBlockShapes(offsets.back(), resIndex, resBlocks);

  typedef BlockTransposeOption<
      OtherTraits::is_symmetric &&
          !(BlockTraits<typename OtherTraits::BlockType>::is_self_transpose),
      Transpose>
      RhsGetter;

  // II - Proper addition

#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for
#endif
  for (std::ptrdiff_t i = 0; i < (std::ptrdiff_t)nonZeros.size(); ++i) {
    for (unsigned j = 0; j < nonZeros[i].size(); ++j) {
      const BlockPtr ptr = offsets[i] + j;
      const NonZero &nz = nonZeros[i][j];

      BlockRef res = resBlocks[ptr];
      const bool afterDiag =
          ((bool)Traits::is_col_major) == ((bool)OtherTraits::is_col_major)
              ? (nz.first > i)
              : (nz.first < i);
      if (nz.second.first == InvalidBlockPtr) {
        RhsGetter::assign(rhs.block(nz.second.second), res, alpha, afterDiag);
      } else if (nz.second.second == OtherDerived::InvalidBlockPtr) {
        res = block(nz.second.first);
      } else {
        RhsGetter::assign(rhs.block(nz.second.second), res, alpha, afterDiag);
        res += block(nz.second.first);
      }
    }
  }

  clear();
  m_majorIndex.move(resIndex);
  resBlocks.swap(m_blocks);
  m_minorIndex.valid = false;

  Finalizer::finalize(*this);

  return derived();
}

template <typename Derived>
template <typename LhsT, typename RhsT>
Derived &SparseBlockMatrixBase<Derived>::operator=(
    const Addition<LhsT, RhsT> &addition) {
  typedef Addition<LhsT, RhsT> Add;

  // WARNING -- Not safe w.r.t aliasing

  Scaling<typename Add::Lhs::ObjectType> lhs(addition.lhs.object,
                                             addition.lhs.scaling);
  *this = lhs;

  Evaluator<typename Add::Rhs::ObjectType> rhs(addition.rhs.object);
  return add<Add::transposeRhs>(*rhs, addition.rhs.scaling);
}

template <typename Derived>
template <typename Expression>
Derived &SparseBlockMatrixBase<Derived>::operator=(
    const NarySum<Expression> &sum) {
  if (sum.empty()) {
    setZero();
  } else {
    // WARNING -- Not safe w.r.t aliasing
    typedef typename NarySum<Expression>::Sum::const_iterator Iterator;
    Iterator it = sum.members.begin();

    *this = *it;
    for (; it != sum.members.end(); ++it) {
      *this += *it;
    }
  }

  return derived();
}

}  // namespace bogus

#endif
