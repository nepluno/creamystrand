/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_FWD_HPP
#define BOGUS_BLOCK_FWD_HPP

#include "Block/Constants.hpp"
#include "Block/Traits.hpp"

namespace bogus {

template <typename Derived>
struct BlockObjectBase;

template <typename Derived>
class IterableBlockObject;

template <typename Derived>
struct Transpose;

template <typename Derived>
class BlockMatrixBase;

template <typename Derived>
class SparseBlockMatrixBase;

template <typename BlockT, int Flags = flags::NONE>
class SparseBlockMatrix;

template <typename BlockT, int Flags = flags::NONE>
class FlatSparseBlockMatrix;

template <typename BlockT, int Flags = flags::NONE,
          typename Index = BOGUS_DEFAULT_INDEX_TYPE>
class MappedSparseBlockMatrix;

}  // namespace bogus

#endif
