/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_BLOCK_OPERATORS_HPP
#define BOGUS_BLOCK_OPERATORS_HPP

#include "Expressions.hpp"

// operators +, -, *, /

namespace bogus {

template < typename LhsT, typename RhsT >
Addition< LhsT, RhsT > operator+ ( const BlockObjectBase< LhsT >& lhs,
								   const BlockObjectBase< RhsT > &rhs )
{
	return Addition< LhsT, RhsT >( lhs.derived(), rhs.derived() ) ;
}

template < typename LhsT, typename RhsT >
Addition< LhsT, RhsT > operator- ( const BlockObjectBase< LhsT >& lhs,
								   const BlockObjectBase< RhsT > &rhs )
{
	return Addition<  LhsT, RhsT >( lhs.derived(), rhs.derived(), 1, -1 ) ;
}

template < typename Derived >
Scaling< Derived > operator* ( const BlockObjectBase< Derived >& lhs,
							   typename Derived::Scalar rhs )
{
	return Scaling< Derived >( lhs.derived(), rhs ) ;
}

template < typename Derived >
Scaling< Derived > operator* ( typename Derived::Scalar lhs ,
							   const BlockObjectBase< Derived >& rhs)
{
	return Scaling< Derived >( rhs.derived(), lhs ) ;
}

template < typename LhsT, typename RhsT  >
Product< LhsT, RhsT > operator* ( const BlockObjectBase< LhsT >& lhs,
								  const BlockObjectBase< RhsT > &rhs )
{
	return Product< LhsT, RhsT >( lhs.derived(), rhs.derived() ) ;
}

template < typename Derived >
Scaling< Derived > operator/ ( const BlockObjectBase< Derived >& lhs,
							   typename Derived::Scalar rhs )
{
	return Scaling< Derived >( lhs.derived(), 1/rhs ) ;
}

} //namespace bogus

#endif
