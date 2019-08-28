/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_EIGEN_BLOCK_CONTAINERS_HPP
#define BOGUS_EIGEN_BLOCK_CONTAINERS_HPP

#include <Eigen/StdVector>
#include "../Utils/CppTools.hpp"

namespace bogus {

// We do not have to use the specialized allocator if the size is dynamic or not a multiple of 16 bytes
template< typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols >
struct ResizableSequenceContainer< Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> >
{
	typedef Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> BlockType ;
	typedef typename TypeSwapIf< 
		Rows == Eigen::Dynamic || Cols == Eigen::Dynamic || 0 != ( static_cast< std::size_t >( Rows * Cols * sizeof( Scalar ) ) & 0xf ),
		std::vector< BlockType, Eigen::aligned_allocator<BlockType> >,
		std::vector< BlockType > >::First Type ;
} ;


} //namespace bogus

#endif
