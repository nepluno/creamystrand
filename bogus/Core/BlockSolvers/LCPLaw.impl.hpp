/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_LCPLAW_IMPL_HPP
#define BOGUS_LCPLAW_IMPL_HPP

#include "LCPLaw.hpp"
#include "../Utils/NumTraits.hpp"

namespace bogus {

template< typename Scalar >
void LCPLaw< Scalar >::projectOnConstraint( const unsigned , typename Traits::Vector &x ) const
{
	x[0] = std::max( x[0], 0. ) ;
}

template< typename Scalar >
Scalar LCPLaw< Scalar >::eval( const unsigned ,
							   const typename Traits::Vector &x,
							   const typename Traits::Vector &y ) const
{
	const Scalar fb = x[0] + y[0] - std::sqrt( x[0]*x[0] + y[0]*y[0] ) ;
	return fb * fb ;
}

template< typename Scalar >
bool LCPLaw<Scalar>::solveLocal(
		const unsigned ,
		const typename Traits::Matrix &A,
		const typename Traits::Vector &b,
		typename Traits::Vector &x,
		const Scalar
		) const
{
	if( b[0] >=0. ) {
		x[0] = 0. ;
		return true ;
	}

	if( A(0,0) > NumTraits< Scalar >::epsilon() ) {
		x[0] = -b[0] / A(0,0) ;
		return true ;
	}

	x[0] = std::max( x[0], 0. ) ;
	return false;

}

} // bogus


#endif
