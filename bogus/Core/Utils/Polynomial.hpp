/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_POLYNOMIAL_HPP
#define BOGUS_POLYNOMIAL_HPP

namespace bogus
{

namespace polynomial {

//! Filter values for getRealRoots() functions
enum RealRootsFilter
{
	StrictlyPositiveRoots,
	StrictlyNegativeRoots,
	AllRoots
} ;

template< unsigned Dimension, typename Scalar >
struct RootsFinder
{
	static unsigned getRealRoots( const Scalar* coeffs,
								  Scalar* realRoots,
								  RealRootsFilter filter = AllRoots ) ;
} ;

template< unsigned Dimension, typename Scalar >
struct PossiblyDegenerateRootsFinder
{
	static unsigned getRealRoots( Scalar* coeffs,
								  Scalar* realRoots,
								  RealRootsFilter filter = AllRoots ) ;
} ;

//! Finds the real roots of a polynomial of dimension \p Dimension whose leading coeffcient is 1
/*! \param coeffs the polynomial's coefficients, ordered by increasing degree [ x^0, ... , x^(Dimension-1) ]
	\param the real roots, in arbitrary order
	\param filter An optional filter
	\return the number of real roots found
*/
template< unsigned Dimension, typename Scalar >
unsigned getRealRoots( const Scalar (&coeffs)[Dimension],
					   Scalar (&realRoots)[Dimension],
					   RealRootsFilter filter = AllRoots )
{
	return RootsFinder< Dimension, Scalar >::getRealRoots( coeffs, realRoots, filter ) ;
}

//! Finds the real roots of a polynomial of dimension \p Dimension with arbitrary leading coefficient
/*! \param coeffs the polynomial's coefficients, ordered by increasing degree [ x^0, ... , x^(Dimension) ]
	\param the real roots, in arbitrary order
	\param filter An optional filter
	\return the number of real roots found
*/
template< unsigned Dimension, typename Scalar >
unsigned getRealRoots( Scalar (&coeffs)[Dimension+1],
					   Scalar (&realRoots)[Dimension],
					   RealRootsFilter filter = AllRoots )
{
	return PossiblyDegenerateRootsFinder< Dimension, Scalar >::getRealRoots( coeffs, realRoots, filter ) ;
}

} //namespace polynomial

} //namespace bogus

#endif
