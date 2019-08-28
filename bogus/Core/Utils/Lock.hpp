/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_LOCK_HPP
#define BOGUS_LOCK_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

namespace bogus {

namespace lock_impl {

#ifdef _OPENMP
typedef omp_lock_t RawData ;
#else
typedef void*      RawData ;
#endif

struct AbstractLock
{
	mutable union {
		RawData api_type ;
		unsigned char padding[24] ;
	} data;

	void (*set_f  )(RawData*) ;
	void (*unset_f)(RawData*) ;

	void set() const {
		if( set_f ) set_f( &data.api_type ) ;
	}

	void unset() const {
		if( unset_f ) unset_f( &data.api_type ) ;
	}

} ;

} //namespace lock_impl

#ifndef _OPENMP
class Lock {
public:
	lock_impl::AbstractLock  for_abi_compat ;

	Lock()
	{
		for_abi_compat.  set_f = 0 ;
		for_abi_compat.unset_f = 0 ;
	}

	template< bool DoLock = true >
	struct Guard {
		explicit Guard( Lock& ) {}
		~Guard() {} 
	} ;
};
#else

class Lock {

public:
	template< bool DoLock = true >
	struct Guard {
		explicit Guard( const Lock& lock )
			: m_lock( lock )
		{
			if(DoLock) m_lock.set() ;
		}

		~Guard()
		{
			if(DoLock) m_lock.unset() ;
		}

	private:
		Guard(const Guard &guard) ;
		Guard& operator=(const Guard &guard) ;

		const Lock &m_lock ;
	} ;

	Lock()
	{
		init() ;
	}

	Lock( const Lock& )
	{
		init() ;
	}

	Lock& operator=( const Lock& )
	{
		return *this ;
	}

	~Lock()
	{
		omp_destroy_lock( data() ) ;
	}

	void set() const {
		m_impl.set() ;
	}

	void unset() const {
		m_impl.unset() ;
	}

	lock_impl::RawData* data() const {
		return &m_impl.data.api_type ;
	}

private:

	void init()
	{
		m_impl.  set_f = &omp_set_lock ;
		m_impl.unset_f = &omp_unset_lock ;
		omp_init_lock( data() ) ;
	}

	lock_impl::AbstractLock m_impl ;
};

#endif

} //namespace bogus

#endif
