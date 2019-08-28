/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_BLOCK_SCALAR_BINDGINS
#define BOGUS_BLOCK_SCALAR_BINDGINS

#include "Expressions.hpp"
#include <cmath>
#include <cstdlib>

#define BOGUS_BLOCK_SCALAR_TYPES \
	BOGUS_PROCESS_SCALAR( double   ) \
	BOGUS_PROCESS_SCALAR( float    ) \
	BOGUS_PROCESS_SCALAR( int      ) \
	BOGUS_PROCESS_SCALAR( char     ) \

namespace bogus {

#define BOGUS_PROCESS_SCALAR( Scalar ) \
	inline bool is_zero( Scalar s, Scalar precision ) { return std::abs( s ) <= precision ; } \
	inline void set_zero( Scalar &s ) { s = 0 ; } \
	inline void set_identity( Scalar &s ) { s = 1 ; } \
	inline void resize( Scalar &, int ,int ) { } \
	inline const Scalar* data_pointer( const Scalar &s ) { return &s ; }
BOGUS_BLOCK_SCALAR_TYPES
#undef BOGUS_PROCESS_SCALAR

#define BOGUS_PROCESS_SCALAR( Scalar ) \
BOGUS_BLOCK_SCALAR_TYPES
#undef BOGUS_PROCESS_SCALAR

#define BOGUS_PROCESS_SCALAR( Scalar_ ) \
	template< > struct BlockTraits< Scalar_ > { \
		typedef Scalar_ Scalar ;  \
		typedef Scalar_ TransposeStorageType ;  \
		enum { RowsAtCompileTime = 1, ColsAtCompileTime = 1, \
				 is_row_major = 0, uses_plain_array_storage = 1, \
				 is_self_transpose = 1 } ; \
	} ;
BOGUS_BLOCK_SCALAR_TYPES
#undef BOGUS_PROCESS_SCALAR

}


#endif
