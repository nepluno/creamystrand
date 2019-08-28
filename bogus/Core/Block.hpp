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

#include "Block/ScalarBindings.hpp"
#include "Block/BlockObjectBase.hpp"
#include "Block/BlockMatrixBase.hpp"

#ifndef BOGUS_WITHOUT_EIGEN
#include "Eigen/EigenBlockContainers.hpp"
#include "Eigen/BlockBindings.hpp"
#endif

#include "Block/Evaluators.hpp"

#include "Block/SparseBlockMatrix.hpp"
#include "Block/MappedSparseBlockMatrix.hpp"
#include "Block/FlatSparseBlockMatrix.hpp"
#include "Block/CompoundMatrix.hpp"
#include "Block/Zero.hpp"
#include "Block/Operators.hpp"

#endif
