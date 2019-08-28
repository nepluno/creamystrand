/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_NUMTRAITS_HPP
#define BOGUS_NUMTRAITS_HPP

#include <limits>
#include <cmath>

namespace bogus
{

template <typename Scalar>
struct NumTraits
{
	static Scalar epsilon()
	{ return std::numeric_limits< Scalar >::epsilon() ; }
	static bool isZero( Scalar s )
	{ return std::fabs(s) < std::numeric_limits< Scalar >::epsilon() ; }
	static bool isSquareZero( Scalar s )
	{ return s*s < std::numeric_limits< Scalar >::epsilon() ; }
} ;

template <typename MatrixType>
struct MatrixTraits ;

}

#endif
