/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_KRYLOV_IMPL_HPP
#define BOGUS_KRYLOV_IMPL_HPP

#include "Krylov.hpp"

#include "Preconditioners.impl.hpp"
#include "KrylovMethods.impl.hpp"

namespace bogus {



template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
Krylov< BlockMatrixType, PreconditionerType >::Krylov(
		const BlockObjectBase< BlockMatrixType > & matrix )
	: Base( &matrix, matrix.rows(), NumTraits< Scalar >::epsilon() )
{
	setMatrix( matrix ) ;
}

template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
Krylov< BlockMatrixType, PreconditionerType >::Krylov()
	: Base( BOGUS_NULL_PTR(const BlockObjectBase< BlockMatrixType >),
			100, NumTraits< Scalar >::epsilon() )
{
}

template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
Krylov< BlockMatrixType, PreconditionerType >&
Krylov< BlockMatrixType, PreconditionerType >::setMatrix(
		const BlockObjectBase< BlockMatrixType > & matrix )
{
	Base::m_matrix = &matrix ;
	m_preconditioner.setMatrix( matrix.derived() ) ;

	return *this ;
}

template < typename BlockMatrixType, template< typename BlockMatrixT > class PreconditionerType >
template < typename RhsT, typename ResT >
typename Krylov< BlockMatrixType, PreconditionerType >::Scalar
Krylov< BlockMatrixType, PreconditionerType >::solve( const RhsT &b, ResT &x,
																	 krylov::Method method ) const
{
	switch(method)
	{
#define BOGUS_PROCESS_KRYLOV_METHOD( MethodName ) \
	case krylov::MethodName : return solve_##MethodName( b, x ) ;
BOGUS_KRYLOV_METHODS
#undef BOGUS_PROCESS_KRYLOV_METHOD
	}
	return -1 ;
}

} // namespace bogus

#endif
