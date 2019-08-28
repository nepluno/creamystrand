/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "CppTools.hpp"

#if BOGUS_HAS_CPP11 && !defined(BOGUS_SHARED_PTR_NS)
#define BOGUS_SHARED_PTR_NS std
#endif

// Standard shared_ptr -- BOGUS_SHARED_PTR_NS should be set to std, std::tr1 or boost
#ifdef BOGUS_SHARED_PTR_NS
#ifndef BOGUS_SHARED_PTR
#include <memory>
#define BOGUS_SHARED_PTR( Type, Name ) BOGUS_SHARED_PTR_NS::shared_ptr< Type > Name
#endif
#endif

#ifndef BOGUS_NAIVE_SHARED_PTR_HPP
#define BOGUS_NAIVE_SHARED_PTR_HPP

#ifndef BOGUS_SHARED_PTR
#define BOGUS_SHARED_PTR( Type, Name ) NaiveSharedPtr< Type > Name
#endif

#if defined( _MSC_VER ) && !defined( BOGUS_DONT_USE_BUILTIN_ATOMICS )
	#define BOGUS_DONT_USE_BUILTIN_ATOMICS
#endif

namespace bogus {

//! Naive reference-counting Shared Pointer
/*! There's no reason to use this implementation over any other one from boost, tr1, or c++11,
  except for c++98 compability. Fortunately it is only used for an edge case ( see EigenSparseLinearSolvers.hpp )

  \warning This class is *NOT* thread-safe. However, multiple NaiveSmartPointers pointing
  to the same shared instance should be able to be safely reset in parallel.

*/
template < typename T >
class NaiveSharedPtr
{
private:

	T* m_instance ;
	int* m_refCount ;

public:

	explicit NaiveSharedPtr( T * instance = BOGUS_NULL_PTR(T) )
	{
		acquire( instance ) ;
	}

	NaiveSharedPtr( const NaiveSharedPtr< T >& rhs )
	{
		add_ref( rhs ) ;
	}

	const NaiveSharedPtr& operator=( const NaiveSharedPtr< T >& rhs )
	{
		if( this != &rhs )
		{
			release() ;
			add_ref( rhs ) ;
		}
		return *this ;
	}

	~NaiveSharedPtr()
	{
		release() ;
	}

	void reset( T * instance = BOGUS_NULL_PTR(T) )
	{
		release() ;
		acquire( instance ) ;
	}

	void release()
	{
		T* instance = BOGUS_NULL_PTR(T) ;
		std::swap( m_instance, instance ) ;


		if( instance && 0 == sync_add< -1 > () )
		{
			delete instance ;
			delete m_refCount ;
		}
	}

	T& operator * ( ) {
		return *m_instance;
	}
	T const& operator * ( ) const {
		return *m_instance;
	}
	T* operator -> ( ) {
		return m_instance;
	}
	T const* operator -> ( ) const {
		return m_instance;
	}
	T* get() {
		return m_instance ;
	}
	const T* get() const {
		return m_instance ;
	}

	operator bool() const {
		return BOGUS_NULL_PTR(T) != m_instance ;
	}

private:

	T* add_ref() const
	{
		T* instance = m_instance ;

		if( instance && sync_add< 1 >() <= 1 )
		{
			// Here m_refCount may already been deleted if another thread is simultaneously releasing
			// this object. We hope that the memory location is still accessible and check for destruction
			// But as I said previously, this class is *NOT* thread-safe
			instance = BOGUS_NULL_PTR(T) ;
		}

		return instance ;
	}

	void add_ref( const NaiveSharedPtr< T >& rhs )
	{
		m_refCount = rhs.m_refCount ;
		m_instance = rhs.add_ref() ;
	}

	void acquire( T* instance )
	{
		if( instance ) {
			m_refCount = new int ;
			*m_refCount = 1 ;
		}
		m_instance = instance ;
	}


	template< int val >
	inline int sync_add() const
	{

#ifndef BOGUS_DONT_USE_BUILTIN_ATOMICS
		return __sync_add_and_fetch( m_refCount, val );
#else
		int t;

#ifdef OPENMP_3_1
#pragma omp atomic capture
#else
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp critical (BogusNaiveSharedPtr)
#endif
#endif
		{  *m_refCount += val ; t = *m_refCount ; }
		return t;
#endif
	}

} ;

} //naemspace bogus

#endif
