/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2015 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_BLOCK_EVALUATORS_HPP
#define BOGUS_BLOCK_EVALUATORS_HPP

#include "Traits.hpp"

namespace bogus {

//!  Evaluates an expression inside a temporary if necessary, otherwise returns directly a matrix reference
template< typename Src, typename Dest = typename BlockMatrixTraits<Src>::PlainObjectType >
struct Evaluator ;

template< typename Src, typename Dest >
struct Evaluator {
	Dest val ;

	Evaluator( const Src& src )
		: val(src)
	{}

	const Dest& operator * ( ) const {
		return val;
	}
	const Dest* operator -> ( ) const {
		return &val;
	}
};


template< typename Src >
struct Evaluator< Src, Src > {
	const Src &val ;

	Evaluator( const Src& src )
		: val(src)
	{}

	const Src& operator * ( ) const {
		return val;
	}
	const Src* operator -> ( ) const {
		return &val;
	}
};

template< typename Src, typename Dest >
struct Evaluator< Transpose<Src>, Dest >
: public Evaluator< Src, Dest >
{
	Evaluator( const Transpose< Src >& src )
		: Evaluator< Src, Dest >( src.matrix )
	{}
} ;

}

#endif
