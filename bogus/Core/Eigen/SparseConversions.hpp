/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_EIGEN_SPARSE_CONVERSIONS_HPP
#define BOGUS_EIGEN_SPARSE_CONVERSIONS_HPP

#include <map>

#include "../Block/Access.hpp"
#include "../Block/SparseBlockIndexComputer.hpp"
#include "../Block/SparseBlockMatrix.hpp"
#include "SparseHeader.hpp"

// Conversions to Sparse Matrix
namespace bogus {

template <typename EigenDerived, typename BogusDerived>
void convert(const Eigen::SparseMatrixBase<EigenDerived>& source,
             SparseBlockMatrixBase<BogusDerived>& dest,
             int destRowsPerBlock = 0, int destColsPerBlock = 0) {
  typedef BlockMatrixTraits<BogusDerived> Traits;
  typedef typename Traits::Index Index;
  typedef typename Traits::BlockPtr BlockPtr;

  const Index RowsPerBlock =
      destRowsPerBlock ? (Index)destRowsPerBlock : (Index)Traits::RowsPerBlock;
  const Index ColsPerBlock =
      destColsPerBlock ? (Index)destColsPerBlock : (Index)Traits::ColsPerBlock;

  assert(RowsPerBlock != (Index)-1);
  assert(ColsPerBlock != (Index)-1);

  BOGUS_STATIC_ASSERT(
      (((bool)Eigen::SparseMatrixBase<EigenDerived>::IsRowMajor) ^
       ((bool)Traits::is_col_major)),
      MATRICES_ORDERING_IS_INCONSISTENT);

  assert(0 == (source.rows() % RowsPerBlock));
  assert(0 == (source.cols() % ColsPerBlock));

  dest.clear();
  dest.setRows(source.rows() / RowsPerBlock, RowsPerBlock);
  dest.setCols(source.cols() / ColsPerBlock, ColsPerBlock);
  dest.reserve(
      1 + (source.derived().nonZeros() / (RowsPerBlock * ColsPerBlock)),
      source.derived().nonZeros());

  const Index blockSize = Traits::is_col_major ? ColsPerBlock : RowsPerBlock;
  const Index innerBlockSize =
      Traits::is_col_major ? RowsPerBlock : ColsPerBlock;
  for (Index outer = 0; outer < dest.majorIndex().outerSize(); ++outer) {
    // I - compute non-zero blocks
    std::map<Index, BlockPtr> nzBlocks;

    for (Index i = 0; i < blockSize; ++i) {
      for (typename EigenDerived::InnerIterator innerIt(source.derived(),
                                                        outer * blockSize + i);
           innerIt; ++innerIt) {
        const Index blockId = (Index)(innerIt.index()) / innerBlockSize;
        if (Traits::is_symmetric && blockId > outer) break;
        nzBlocks[blockId] = 0;
      }
    }

    // II - Insert them in block mat
    for (typename std::map<Index, BlockPtr>::iterator bIt = nzBlocks.begin();
         bIt != nzBlocks.end(); ++bIt) {
      bIt->second = (BlockPtr)dest.nBlocks();
      typename BogusDerived::BlockRef block =
          dest.template insertByOuterInner<true>(outer, bIt->first);
      resize(block, RowsPerBlock, ColsPerBlock);
      block.setZero();
    }

    // III - copy values
    for (Index i = 0; i < blockSize; ++i) {
      for (typename EigenDerived::InnerIterator innerIt(source.derived(),
                                                        outer * blockSize + i);
           innerIt; ++innerIt) {
        const Index blockId = (Index)(innerIt.index()) / innerBlockSize;
        if (Traits::is_symmetric && blockId > outer) break;
        const Index binn = innerIt.index() - blockId * innerBlockSize;
        const Index brow = Traits::is_col_major ? binn : i;
        const Index bcol = Traits::is_col_major ? i : binn;

        dest.block(nzBlocks[blockId])(brow, bcol) = innerIt.value();
      }
    }

    // IV - Symmetrify diagonal block if required
    if (Traits::is_symmetric) {
      typename std::map<Index, BlockPtr>::const_iterator diagPtr =
          nzBlocks.find(outer);
      if (diagPtr != nzBlocks.end()) {
        const typename Traits::BlockType diagBlock =
            dest.block(diagPtr->second);
        dest.block(diagPtr->second) =
            .5 *
            (diagBlock + TransposeIf<Traits::is_symmetric>::get(diagBlock));
      }
    }
  }

  dest.finalize();
}

template <typename BogusDerived, typename EigenScalar, int EigenOptions,
          typename EigenIndex>
void convert(const SparseBlockMatrixBase<BogusDerived>& source,
             Eigen::SparseMatrix<EigenScalar, EigenOptions, EigenIndex>& dest) {
  typedef BlockMatrixTraits<BogusDerived> Traits;
  typedef typename Traits::Index Index;

  typedef Eigen::SparseMatrix<EigenScalar, EigenOptions, EigenIndex>
      EigenMatrixType;

  typedef SparseBlockIndexGetter<BogusDerived,
                                 Traits::is_symmetric ||
                                     ((bool)EigenMatrixType::IsRowMajor) !=
                                         ((bool)Traits::is_col_major)>
      IndexGetter;

  dest.setZero();
  dest.resize(source.rows(), source.cols());

  typename BogusDerived::UncompressedIndexType auxIndex;
  const typename IndexGetter::ReturnType& index =
      IndexGetter::getOrCompute(source, auxIndex);

  const std::vector<Index>& outerOffsets =
      ((bool)EigenMatrixType::IsRowMajor) == ((bool)Traits::is_col_major)
          ? source.majorIndex().innerOffsets
          : source.minorIndex().innerOffsets;

  for (Index outerBlock = 0; outerBlock < index.outerSize(); ++outerBlock) {
    for (Index outer = outerOffsets[outerBlock];
         outer < outerOffsets[outerBlock + 1]; ++outer) {
      dest.startVec(outer);
      for (typename IndexGetter::ReturnType::InnerIterator it(index,
                                                              outerBlock);
           it; ++it) {
        for (Index inner = index.innerOffsets[it.inner()];
             inner < index.innerOffsets[it.inner() + 1]; ++inner) {
          if (Traits::is_symmetric && inner > outer) break;

          const Index bout = outer - outerOffsets[outerBlock];
          const Index binn = inner - index.innerOffsets[it.inner()];
          const Index bcol = EigenMatrixType::IsRowMajor ? binn : bout;
          const Index brow = EigenMatrixType::IsRowMajor ? bout : binn;
          const EigenScalar val = source.block(it.ptr())(brow, bcol);

          dest.insertBackByOuterInner(outer, inner) = val;
        }
      }
    }
  }
}

}  // namespace bogus

#endif
