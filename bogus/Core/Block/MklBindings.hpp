/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_BLOCK_MKL_BINDINGS_HPP
#define BOGUS_BLOCK_MKL_BINDINGS_HPP

#include "SparseBlockMatrix.hpp"
#include "CompressedSparseBlockIndex.hpp"

#include <mkl.h>
#include <mkl_spblas.h>

// Creates a compile error with boost
#ifdef P4
#undef P4
#endif

namespace bogus
{

//! Specializations using mkl
namespace mkl
{

//! Wrapper over scalar-specific mkl calls
template< typename Scalar >
struct bindings {} ;

template< >
struct bindings< double >
{
	typedef double Scalar ;
	static void bsrmv (char *transa, MKL_INT *m, MKL_INT *k, MKL_INT *lb, Scalar *alpha, char *matdescra, Scalar  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, Scalar *x, Scalar *beta, Scalar *y)
	{
		mkl_dbsrmv( transa, m, k, lb, alpha,
					matdescra, val, indx, pntrb, pntre,
					x, beta, y ) ;
	}
} ;

/* FIXME sbsrmv no supported on older MKL versions ; do version detection
template< >
struct bindings< float >
{
	typedef float Scalar ;
	static void bsrmv (char *transa, MKL_INT *m, MKL_INT *k, MKL_INT *lb, Scalar *alpha, char *matdescra, Scalar  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, Scalar *x, Scalar *beta, Scalar *y)
	{
		mkl_sbsrmv( transa, m, k, lb, alpha,
					matdescra, val, indx, pntrb, pntre,
					x, beta, y ) ;
	}
} ; */


template< typename Scalar >
void bsrmv(
		bool Symmetric, bool Transpose, MKL_INT Dimension,
		MKL_INT m, MKL_INT k, const std::size_t offset,
		const MKL_INT *rowIndex, const MKL_INT *columns, const Scalar *data,
		const Scalar *rhs, int rhsCols, Scalar *res, Scalar alpha, Scalar beta )
{

	char matdescra[4] = { Symmetric ? 'S' : 'G', 'L', 'N', 'C'} ;

	char transa = Transpose ? 'T' : 'N' ;
	MKL_INT lb = Dimension ;
	Scalar *x = const_cast< Scalar* > ( rhs ) ;
	Scalar *y = res ;

	MKL_INT* pntrb = const_cast< MKL_INT* >( rowIndex+offset ) ;
	MKL_INT* pntre = const_cast< MKL_INT* >( rowIndex+offset+1 ) ;
	MKL_INT* indx  = const_cast< MKL_INT* >( columns )  ;

	Scalar *a = const_cast< Scalar* > ( data ) + pntrb[0] * lb * lb ;

	for( int i = 0 ; i < rhsCols ; ++i )
	{
		bindings< Scalar >::bsrmv( &transa, &m, &k, &lb, &alpha,
								   matdescra, a, indx, pntrb, pntre,
								   x + i*Dimension*k, &beta, y + i*Dimension*k ) ;
	}

}

template < typename BlockPtr, typename BlockType >
static void rowmv( const SparseBlockIndex< true, MKL_INT, BlockPtr >& index,
				   const BlockType* data, MKL_INT row,
				   const typename BlockTraits< BlockType >::Scalar *rhs, int rhsCols,
				   typename BlockTraits< BlockType >::Scalar *res )
{
	const MKL_INT *rowIndexOrig  = index.rowIndex() ;
	const MKL_INT  offset        = rowIndexOrig[ row ] ;
	const MKL_INT  rowIndex[2]   = {0, rowIndexOrig[ row+1 ] - offset} ;
	const MKL_INT *columns       = index.columns() + offset ;

	mkl::bsrmv< typename BlockTraits< BlockType >::Scalar >
			( false, false, BlockTraits< BlockType >::RowsAtCompileTime,
			  1, index.innerSize(), 0,
			  rowIndex, columns,
			  data_pointer( data[offset + index.base] ),
			  rhs, rhsCols, res, 1, 1 ) ;
}

} //namespace mkl

template <>
struct SparseBlockMatrixOpProxy< true, true, double, MKL_INT >
{
	typedef double Scalar ;

	template < bool Transpose, typename Derived, typename RhsT, typename ResT >
	static void multiply( const SparseBlockMatrixBase< Derived >& matrix,  const RhsT& rhs, ResT& res,
						  Scalar alpha, Scalar beta )
	{
		typedef BlockMatrixTraits< Derived > Traits ;

		mkl::bsrmv< Scalar >
				( Traits::is_symmetric, Transpose, Derived::RowsPerBlock,
				  matrix.rowsOfBlocks(), matrix.colsOfBlocks(), 0,
				  matrix.majorIndex().rowIndex(), matrix.majorIndex().columns(),
				  data_pointer( matrix.data()[0] ),
				  rhs.data(), rhs.cols(), res.data(), alpha, beta ) ;
	}

	template < typename Derived, typename RhsT, typename ResT >
	static void splitRowMultiply( const SparseBlockMatrixBase< Derived >& matrix, typename Derived::Index row, const RhsT& rhs, ResT& res  )
	{
		typedef BlockMatrixTraits< Derived > Traits ;

		if( Traits::is_symmetric && !matrix.transposeIndex().valid )
		{
			SparseBlockSplitRowMultiplier< Traits::is_symmetric, !Traits::is_col_major >
					::splitRowMultiply( matrix, row, rhs, res )  ;

			return ;
		}

		mkl::rowmv( matrix.majorIndex(), matrix.data(), row,
					rhs.data(), rhs.cols(), res.data() ) ;

		if( Traits::is_symmetric )
		{
			mkl::rowmv( matrix.transposeIndex(), matrix.data(), row,
						rhs.data(), rhs.cols(), res.data() ) ;
		}

		// Remove diagonal block if it exist
		const typename Traits::BlockPtr diagPtr = matrix.diagonalBlockPtr( row ) ;
		if( diagPtr != matrix.InvalidBlockPtr )
		{

			const Segmenter< BlockDims< typename Derived::BlockType, false >::Cols, const RhsT, typename Derived::Index >
					segmenter( rhs, matrix.colOffsets() ) ;
			res -= matrix.block( diagPtr ) * segmenter[row] ;
		}

	}
} ;


} //namespace bogus

#endif
