/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_BLOCK_STREAMS_HPP
#define BOGUS_BLOCK_STREAMS_HPP

#include <iostream>
#include "SparseBlockMatrix.hpp"
#include "SparseBlockIndexComputer.hpp"

template < typename Derived >
std::ostream& operator<<( std::ostream &out, const bogus::SparseBlockMatrixBase< Derived > &sbm )
{

	typedef bogus::SparseBlockIndexComputer< Derived, false, false > IndexComputerType ;
	IndexComputerType indexComputer( sbm ) ;
	typedef typename IndexComputerType::ReturnType SourceIndexType ;
	const SourceIndexType &sourceIndex = indexComputer.get() ;

	out << " Total rows: " << sbm.rows() << " / cols: " << sbm.cols() << std::endl ;
	for ( unsigned i = 0 ; i < (unsigned) sourceIndex.outerSize() ; ++ i )
	{
		out << "Row " << i << ": " ;
		for( typename SourceIndexType::InnerIterator it( sourceIndex, i ) ;
			 it ; ++ it )
		{
			out << " " << it.inner() << "@" << it.ptr() << "; " ;
		}
		out << std::endl ;
	}
	out << " Blocks (" << sbm.nBlocks() << ")" << std::endl ;
	for ( unsigned i = 0 ; i < sbm.nBlocks() ; ++ i )
	{
		out << sbm.block(i) << std::endl ;
		out << "^-- " << i << std::endl ;
	}
	return out ;
}

#endif
