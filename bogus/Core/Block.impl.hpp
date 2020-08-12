/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_IMPL_HPP
#define BOGUS_BLOCK_IMPL_HPP

#include "Block.hpp"
#include "Block/CompoundMatrix.impl.hpp"
#include "Block/SparseAssign.impl.hpp"
#include "Block/SparseBlockIndex.impl.hpp"
#include "Block/SparseBlockMatrixBase.impl.hpp"
#include "Block/SparseMatrixMatrixProduct.impl.hpp"
#include "Block/SparseMatrixVectorProduct.impl.hpp"
#include "Block/SparseScaleAdd.impl.hpp"

#ifdef BOGUS_WITH_MKL
#include "Block/MklBindings.hpp"
#endif

#endif
