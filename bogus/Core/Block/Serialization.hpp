/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_SERIALIZATION_HPP
#define BOGUS_BLOCK_SERIALIZATION_HPP

#include "SparseBlockMatrix.hpp"

namespace boost {
namespace serialization {

template <typename Index, typename BlockPtr>
struct version<bogus::SparseBlockIndex<true, Index, BlockPtr> > {
  enum { value = 1 };
};

template <typename BlockT, int Flags>
struct version<bogus::SparseBlockMatrix<BlockT, Flags> > {
  enum { value = 1 };
};

template <typename Archive, typename Index, typename BlockPtr>
inline void serialize(Archive& ar,
                      bogus::SparseBlockIndex<false, Index, BlockPtr>& index,
                      const unsigned int file_version) {
  (void)file_version;
  ar& index.valid;
  ar& index.innerOffsets;
  ar& index.outer;
}

template <typename Archive, typename Index, typename BlockPtr>
inline void serialize(Archive& ar,
                      bogus::SparseBlockIndex<true, Index, BlockPtr>& index,
                      const unsigned int file_version) {
  BlockPtr dummy_base;

  ar& index.valid;
  ar& index.innerOffsets;
  if (file_version == 0) ar& dummy_base;
  ar& index.inner;
  ar& index.outer;
}

}  // namespace serialization
}  // namespace boost

namespace bogus {

template <typename Derived>
template <typename Archive>
void SparseBlockMatrixBase<Derived>::serialize(
    Archive& ar, const unsigned int file_version) {
  std::size_t dummyNBlocks(0);

  ar& m_rows;
  ar& m_cols;
  ar& m_blocks;
  if (file_version == 0) ar& dummyNBlocks;
  ar& m_majorIndex;
  ar& m_minorIndex;
  ar& m_transposeIndex;

  if (file_version == 0) {
    m_transposeBlocks.clear();
    m_transposeIndex.clear();
    m_transposeIndex.valid = 0;
  } else
    ar& m_transposeBlocks;
}

}  // namespace bogus

#endif
