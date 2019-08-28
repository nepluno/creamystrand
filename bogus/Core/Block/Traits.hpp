/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_BLOCK_TRAITS_HPP
#define BOGUS_BLOCK_TRAITS_HPP

#include <vector>

namespace bogus {

template< typename Derived >
struct BlockMatrixTraits { } ;

//! Default container type, that should resizable and use contiguous storage
template< typename ElementType >
struct ResizableSequenceContainer {
	typedef std::vector< ElementType > Type ;
	enum { is_mutable = 1 } ;
} ;

template< typename BlockType >
struct BlockTraits
{
   typedef typename BlockType::Scalar Scalar ;
	//! Type for storing the result of transpose_block( BlockType ), useful for cacheTranspose()
   typedef BlockType TransposeStorageType ;

	enum {
	   //! Number of rows spanned by a block at compile time ; useful for efficient segmentation
	   RowsAtCompileTime = BlockType::RowsAtCompileTime,
	   //! Number of cols spanned by a block at compile time ; useful for efficient segmentation
	   ColsAtCompileTime = BlockType::ColsAtCompileTime,

	   //! Can be set to true if "const Scalar* data_pointer( const BlockType& )" exist.
	   uses_plain_array_storage = 0,
	   //! Ordering inside the block ; only useful uses_plain_array_storage is true
	   is_row_major = false,
	   //! Whether this block is equal to its transpose
	   is_self_transpose = 0
	}  ;
} ;

// Transpose, blobk-matrix/vector and block/block product return types
// Specialization of these structures should define a ReturnType if the operation is allowed

//! Defines the return type of an associated transpose_block( const BlockType& ) function
template< typename BlockType >
struct BlockTransposeTraits {} ;

//! Defines the type of the vectors resulting from the multiplication of a BlockMatrix and an instance of VectorTypea
/*! Should be the return type of an associated get_mutable_vector( const VectorType & ) function */
template< typename VectorType >
struct BlockVectorProductTraits {} ;

//! Defines the return type of the product of two blocks potentially transposed
template< typename LhsBlockType, typename RhsBlockType, bool TransposeLhs, bool TransposeRhs >
struct BlockBlockProductTraits {
	typedef LhsBlockType ReturnType ;
} ;

} // namespace bogus

#endif
