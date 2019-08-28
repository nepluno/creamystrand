/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_SPARSE_BLOCK_INDEX_IMPL_HPP
#define BOGUS_SPARSE_BLOCK_INDEX_IMPL_HPP

#include "SparseBlockIndex.hpp"

#include "CompoundSparseBlockIndex.hpp"

namespace bogus {

template< typename Derived >
Derived & SparseBlockIndexBase< Derived >::derived()
{
	return static_cast< Derived& >( *this ) ;
}

template< typename Derived >
const Derived & SparseBlockIndexBase< Derived >::derived() const
{
	return static_cast< const Derived& >( *this ) ;
}

template< typename Derived >
typename SparseBlockIndexBase< Derived >::Index SparseBlockIndexBase< Derived >::innerSize() const
{
	return innerOffsetsArray().size() - 1  ;
}

template< typename Derived >
typename SparseBlockIndexBase< Derived >::Index SparseBlockIndexBase< Derived >::outerSize() const
{
	return derived().outerSize() ;
}

template< typename Derived >
const typename SparseBlockIndexBase< Derived >::InnerOffsetsType & SparseBlockIndexBase< Derived >::innerOffsetsArray() const
{
	return derived().innerOffsetsArray() ;
}

template< typename Derived >
bool SparseBlockIndexBase< Derived >::hasInnerOffsets() const
{
	return !innerOffsetsArray().empty();
}

template< typename Derived >
typename SparseBlockIndexBase< Derived >::Index SparseBlockIndexBase< Derived >::size( Index outerIdx ) const
{
	return derived().size( outerIdx ) ;
}

template< typename Derived >
typename SparseBlockIndexBase< Derived >::InnerIterator SparseBlockIndexBase< Derived >::begin( Index outerIdx ) const
{
	return InnerIterator( *this, outerIdx ) ;
}

template< typename Derived >
typename SparseBlockIndexBase< Derived >::InnerIterator SparseBlockIndexBase< Derived >::end( Index outerIdx ) const
{
	return InnerIterator( *this, outerIdx ).toEnd() ;
}

template< typename Derived >
typename SparseBlockIndexBase< Derived >::InnerIterator SparseBlockIndexBase< Derived >::last( Index outerIdx ) const
{
	InnerIterator it ( derived(), outerIdx ) ;
	return it ? --(it.toEnd()) : it ;
}



template < bool Compressed, typename Index, typename BlockPtr, template <typename> class ArrayType  >
template < typename SourceDerived >
SparseBlockIndex< Compressed, Index, BlockPtr, ArrayType > & SparseBlockIndex< Compressed, Index, BlockPtr, ArrayType >::operator=(
		const SparseBlockIndexBase< SourceDerived > &source )
{
	clear() ;
	resizeOuter( source.outerSize() ) ;

	for( typename SourceDerived::Index i = 0 ; i < source.outerSize() ; ++i )
	{
		for( typename SourceDerived::InnerIterator it( source.derived(), i ) ;
			 it ; ++ it )
		{
			insertBack( i, it.inner(), it.ptr() ) ;
		}
	}

	finalize() ;
	valid = source.valid ;
	if( source.hasInnerOffsets() ) {
		innerOffsets.resize( source.innerOffsetsArray().size() ) ;
		std::copy( source.innerOffsetsArray().begin(), source.innerOffsetsArray().end(), innerOffsets.begin() ) ;
	}

	return *this ;
}


} // namespace bogus

#endif

