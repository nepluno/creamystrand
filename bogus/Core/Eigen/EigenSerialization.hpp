/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_EIGEN_SERIALIZATION_HPP
#define BOGUS_EIGEN_SERIALIZATION_HPP

#include <Eigen/Core>

#ifndef BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE
#include "SparseHeader.hpp"
#endif

namespace boost
{
namespace serialization
{

template<typename Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void load(
	   Archive & ar,
	   Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	ar & make_array( matrix.data(), matrix.size() ) ;
}

template<typename Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void save(
	   Archive & ar,
	   const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	ar & make_array( matrix.data(), matrix.size() ) ;
}


template<typename Archive, typename _Scalar, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void load(
	   Archive & ar,
	   Eigen::Matrix<_Scalar, Eigen::Dynamic, _Cols, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	int rows ;
	ar & rows ;
	matrix.resize( rows, _Cols ) ;
	ar & make_array( matrix.data(), matrix.size() ) ;
}

template<typename Archive, typename _Scalar, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void save(
	   Archive & ar,
	   const Eigen::Matrix<_Scalar, Eigen::Dynamic, _Cols, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	int rows = matrix.rows()  ;
	ar & rows ;
	ar & make_array( matrix.data(), matrix.size() ) ;
}


template<typename Archive, typename _Scalar, int _Rows, int _Options, int _MaxRows, int _MaxCols>
inline void load(
	   Archive & ar,
	   Eigen::Matrix<_Scalar, _Rows, Eigen::Dynamic, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	int cols ;
	ar & cols ;
	matrix.resize( _Rows, cols ) ;
	ar & make_array( matrix.data(), matrix.size() ) ;
}

template<typename Archive, typename _Scalar, int _Rows, int _Options, int _MaxRows, int _MaxCols>
inline void save(
	   Archive & ar,
	   const Eigen::Matrix<_Scalar, _Rows, Eigen::Dynamic, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	int cols = matrix.cols() ;
	ar & cols ;
	ar & make_array( matrix.data(), matrix.size() ) ;
}

template<typename Archive, typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
inline void load(
	   Archive & ar,
	   Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	int rows, cols ;
	ar & rows ;
	ar & cols ;
	matrix.resize( rows, cols ) ;
	ar & make_array( matrix.data(), matrix.size() ) ;
}

template<typename Archive, typename _Scalar, int _Options, int _MaxRows, int _MaxCols>
inline void save(
	   Archive & ar,
	   const Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	int rows = matrix.rows(), cols = matrix.cols() ;
	ar & rows ;
	ar & cols ;
	ar & make_array( matrix.data(), matrix.size() ) ;
}

template<typename Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void serialize(
	   Archive & ar,
	   Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	split_free( ar, matrix, file_version ) ;
}

#ifdef BOGUS_WITH_EIGEN_STABLE_SPARSE_API

template<typename Archive, typename _Scalar, int _Options, typename _Index >
inline void load(
	   Archive & ar,
	   Eigen::SparseMatrix<_Scalar, _Options, _Index> & matrix,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	_Index rows, cols ;
	unsigned nnz ;
	ar & rows ;
	ar & cols ;
	ar & nnz ;
	matrix.resize( rows, cols ) ;
	matrix.resizeNonZeros( nnz ) ;
	ar & make_array( matrix.outerIndexPtr(), matrix.outerSize()+1 ) ;
	ar & make_array( matrix.innerIndexPtr(), nnz ) ;
	ar & make_array( matrix.valuePtr(), nnz ) ;
}

template<typename Archive, typename _Scalar, int _Options, typename _Index >
inline void save(
	   Archive & ar,
	   const Eigen::SparseMatrix<_Scalar, _Options, _Index> & matrix,
	   const unsigned int file_version
   )
{
	(void) file_version ;
	assert( matrix.isCompressed() ) ;

	_Index rows = matrix.rows(), cols = matrix.cols() ;
	unsigned nnz = matrix.data().size() ;
	ar & rows ;
	ar & cols ;
	ar & nnz ;
	ar & make_array( matrix.outerIndexPtr(), matrix.outerSize()+1 ) ;
	ar & make_array( matrix.innerIndexPtr(), nnz ) ;
	ar & make_array( matrix.valuePtr(), nnz ) ;
}

template<typename Archive, typename _Scalar, int _Options, typename _Index >
inline void serialize(
	   Archive & ar,
	   Eigen::SparseMatrix<_Scalar, _Options, _Index> & matrix,
	   const unsigned int file_version
   )
{
	split_free( ar, matrix, file_version ) ;
}

#endif

} // serialization
} // boost

#endif
