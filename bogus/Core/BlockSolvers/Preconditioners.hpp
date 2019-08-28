/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_PRECONDITIONERS_HPP
#define BOGUS_PRECONDITIONERS_HPP

namespace bogus {

//! Trivial ( identity ) preconditioner. Does nothing.
template < typename MatrixType >
class TrivialPreconditioner
{
public:
	void setMatrix( const MatrixType & )
	{}

	template < bool transpose, typename ResT, typename RhsT >
	void apply( const RhsT& rhs, ResT &res ) const
	{
		res = rhs ;
	}
} ;

//! Diagonal preconditioner
/*! Defines the preconditioner matrix \f$ P^{-1} \f$ as a diagonal matrix
	where each of the coefficient is the scalar inverse of the corresponding
	diagonal coefficient in the system matrix
	Works well for diagonally dominant matrices.
	*/
template < typename MatrixType >
class DiagonalPreconditioner
{
} ;

//! Diagonal Block-LU preconditioner
/*! Defines the preconditioner matrix \f$ P^{-1} \f$ as a block-diagonal matrix
	where each diagonal block is the LU factorization of the corresponding
	diagonal block in the system matrix
	*/
template < typename MatrixType >
class DiagonalLUPreconditioner
{
} ;

//! Diagonal Block-LDLT preconditioner
/*! Defines the preconditioner matrix \f$ P^{-1} \f$ as a block-diagonal matrix
	where each diagonal block is the LDLT factorization of the corresponding
	diagonal block in the system matrix
	*/
template < typename MatrixType >
class DiagonalLDLTPreconditioner
{
} ;

//! Matrix preconditioner
/*! Explicitely define the preconditioner with an arbitray matrix*/

template < typename PreconditionerMatrixType >
struct MatrixPreconditioner {

	template < typename MatrixType >
	class Type
	{
		const PreconditionerMatrixType* m_preconditionerMatrix ;

	public:

		Type() : m_preconditionerMatrix( 0 )
		{}

		void setMatrix( const MatrixType & )
		{}

		void setPreconditionerMatrix( const PreconditionerMatrixType & preconditioner )
		{
			m_preconditionerMatrix = &preconditioner ;
		}

		template < bool transpose, typename ResT, typename RhsT >
		void apply( const RhsT& rhs, ResT &res ) const
		{
			m_preconditionerMatrix->template multiply< transpose >( rhs, res ) ;
		}
	} ;
} ;

}

#endif

