/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "Block.hpp"
#include "Block/Streams.hpp"

#ifdef BOGUS_WITH_BOOST_SERIALIZATION
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#ifndef BOGUS_WITHOUT_EIGEN
#include "Eigen/EigenSerialization.hpp"
#endif

#include "Block/Serialization.hpp"
#endif

#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
#include "Eigen/SparseConversions.hpp"
#endif
