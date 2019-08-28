/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_COLORING_IMPL_HPP
#define BOGUS_COLORING_IMPL_HPP

#include "Coloring.hpp"

#include "../Block/SparseBlockMatrixBase.hpp"

namespace bogus {

template < typename Derived >
void Coloring::compute( const SparseBlockMatrixBase< Derived >& matrix )
{

	// Affects each contact to the first color with which it does not have any interaction
	// This is not optimimal, but optimality would be costly

	const std::ptrdiff_t n = static_cast< std::ptrdiff_t >( matrix.rowsOfBlocks() ) ;

	std::vector< std::vector< std::ptrdiff_t > > tmp_colors ;
	for( std::ptrdiff_t i = 0 ; i < n ; ++i )
	{
		bool ok = false ;

		const unsigned nColors = tmp_colors.size() ;

		for( unsigned c = 0 ; c < nColors ; ++c )
		{
			std::vector< std::ptrdiff_t > &color = tmp_colors[ c ] ;
			ok = true ;
			for( typename Derived::MajorIndexType::InnerIterator it( matrix.majorIndex(), i ) ;
				 it && it.inner() < i ; ++ it )
			{
				if( color.end() !=
						std::lower_bound( color.begin(), color.end(), it.inner() ) )
				{
					ok = false ;
					break ;
				}
			}
			if( ok )
			{
				color.push_back( i ) ;
				break ;
			}
		}
		if( !ok )
		{
			tmp_colors.resize( nColors+1 ) ;
			tmp_colors.back().push_back( i ) ;
		}
	}

	permutation.clear() ;

	colors.resize( tmp_colors.size()+1 ) ;
	colors[ 0 ] = 0 ;

	for( unsigned c = 0 ; c < tmp_colors.size() ; ++c )
	{
		colors[ c+1 ] = colors[c] + tmp_colors[c].size() ;
		permutation.insert( permutation.end(), tmp_colors[c].begin(), tmp_colors[c].end() ) ;
	}

}

template < typename Derived >
void Coloring::compute( const BlockMatrixBase< Derived >& matrix )
{ reset( matrix.rowsOfBlocks() ) ; }

template < typename Derived >
void Coloring::update( const bool enable, const BlockMatrixBase< Derived >& matrix )
{
	if( enable )
	{
		compute( matrix.derived() ) ;
	} else {
		reset( matrix.rowsOfBlocks() ) ;
	}
}


} //namespace bogus

#endif
