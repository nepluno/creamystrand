/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_HPP
#define BOGUS_BLOCK_HPP

#include "Block.fwd.hpp"
#include "Block/BlockMatrixBase.hpp"
#include "Block/BlockObjectBase.hpp"
#include "Block/ScalarBindings.hpp"

#ifndef BOGUS_WITHOUT_EIGEN
#include "Eigen/BlockBindings.hpp"
#include "Eigen/EigenBlockContainers.hpp"
#endif

#include "Block/CompoundMatrix.hpp"
#include "Block/Evaluators.hpp"
#include "Block/FlatSparseBlockMatrix.hpp"
#include "Block/MappedSparseBlockMatrix.hpp"
#include "Block/Operators.hpp"
#include "Block/SparseBlockMatrix.hpp"
#include "Block/Zero.hpp"

#endif
