/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_KRYLOV_METHODS_IMPL_HPP
#define BOGUS_KRYLOV_METHODS_IMPL_HPP

#include "KrylovMethods.hpp"

#include "../Utils/NumTraits.hpp"
#include "../Block/Access.hpp"

namespace bogus {

namespace krylov {

namespace solvers {

//! Check init guess, reset it to zero if that would give a lower residual
template < typename Matrix, typename Vector, typename RhsT, typename ResT, typename Scalar >
Scalar init(
		const Matrix& A, const RhsT &b, ResT x, Vector &r0,
		Scalar &scale )
{
	r0 = b ;
	r0.noalias() -= A*x ;

	Scalar res = r0.squaredNorm() ;
	const Scalar resAt0 = b.squaredNorm()  ;

	if( res > resAt0 ) {
		r0 = b;
		x.setZero() ;
		res = resAt0 ;
	}

	scale =  1. / ( 1 + b.size() ) ;
	return res * scale ;
}

template < bool DoTranspose, typename Preconditioner, typename RhsT, typename ResT >
void apply( const Preconditioner *P,  const RhsT &b, ResT& x )
{
	if( P ) P->template apply< DoTranspose >( b, x ) ;
	else    x = b ;
}

// Conjugate Gradient

template < typename Mat, typename Prec, typename Traits >
template < typename RhsT, typename ResT >
typename CG< Mat, Prec, Traits >::Scalar
CG< Mat, Prec, Traits >::vectorSolve( const RhsT &b, ResT x ) const
{
	typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
	Vector r ;

	Scalar scale ;
	Scalar res = init( *m_A, b, x, r, scale ) ;
	if( res < m_tol ) return res ;

	Vector z( m_A->rows() ) ;
	apply< false >( m_P, r, z ) ;
	Vector p = z;

	Scalar zr0 = r.dot( z ) ;
	Scalar zr1 ;

	Vector Mp( m_A->rows() ) ;

	for( unsigned k = 0 ; k < m_maxIters ; ++k )
	{
		Mp = (*m_A) * p ;
		const Scalar alpha = zr0 / ( p.dot( Mp ) ) ;
		x += alpha * p ;
		r -= alpha * Mp ;

		res = r.squaredNorm() * scale ;
		if( m_callback ) m_callback->trigger( k, res ) ;
		if( res < m_tol ) break ;

		apply< false >( m_P, r, z ) ;
		zr1 = z.dot( r ) ;

		p = z + ( zr1 / zr0 ) * p ;

		zr0 = zr1 ;
	}

	return res ;
}

// BiCG
template < typename Mat, typename Prec, typename Traits >
template < typename RhsT, typename ResT >
typename BiCG< Mat, Prec, Traits >::Scalar
BiCG< Mat, Prec, Traits >::vectorSolve( const RhsT &b, ResT x ) const
{
	typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
	Vector r ;

	Scalar scale ;
	Scalar res = init( *m_A, b, x, r, scale ) ;
	if( res < m_tol ) return res ;

	Vector x_ = x ;

	Vector r_ = b ;
	r_.noalias() -= m_A->transpose() * x ;

	Vector p ( m_A->rows() ),
			p_( m_A->rows() ) ;

	apply< false >( m_P, r, p ) ;
	apply< true  >( m_P, r_, p_ ) ;

	Vector z = p, z_ = p_ ;
	Scalar zr0 = r_.dot( z ) ;
	Scalar zr1 ;

	Vector Mp( m_A->rows() ) ;

	for( unsigned k = 0 ; k < m_maxIters ; ++k )
	{
		Mp = (*m_A) * p ;
		const Scalar alpha = zr0 / ( p_.dot( Mp ) ) ;
		x  += alpha * p  ;
		x_ += alpha * p_ ;

		r  -= alpha * Mp ;
		r_.noalias() -= alpha * m_A->transpose() * p_ ;

		res = r.squaredNorm() * scale ;
		if( m_callback ) m_callback->trigger( k, res ) ;
		if( res < m_tol ) break ;

		apply< false >( m_P, r, z ) ;
		zr1 = z.dot( r_ ) ;

		const Scalar beta = ( zr1 / zr0 ) ;

		p = z + beta * p ;
		apply< true  >( m_P, r_, z_ ) ;
		p_ = z_ + beta * p_ ;

		zr0 = zr1 ;
	}


	return res ;
}

// BicGSTAB
template < typename Mat, typename Prec, typename Traits >
template < typename RhsT, typename ResT >
typename BiCGSTAB< Mat, Prec, Traits >::Scalar
BiCGSTAB< Mat, Prec, Traits >::vectorSolve( const RhsT &b, ResT x ) const
{
	typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
	Vector r ;

	Scalar scale ;
	Scalar res = init( *m_A, b, x, r, scale ) ;
	if( res < m_tol ) return res ;

	const Vector r0h = r ;

	Scalar rho0 = 1, rho1 ;
	Scalar alpha = 1, w = 1 ;

	Vector nu;
	nu.setZero( r.rows() );
	Vector p = nu ;
	Vector s ;
	Vector t( m_A->rows() ),
		   y( m_A->rows() ),
		   z( m_A->rows() ) ;

	for( unsigned k = 0 ; k < m_maxIters ; ++k )
	{
		rho1 = r0h.dot( r ) ;

		const Scalar beta = ( rho1 / rho0 ) * ( alpha / w ) ;
		p = r + beta * ( p - w * nu ) ;
		apply< false >( m_P, p, y ) ;
		nu = ( *m_A ) * y ;

		alpha = rho1 / r0h.dot( nu ) ;
		s = r - alpha * nu ;
		apply< false >( m_P, s, z ) ;
		t = ( *m_A ) * z ;

		const Scalar nt2 = t.squaredNorm() ;
		if ( nt2 < NumTraits< Scalar >::epsilon( ) )
		{
			w = 1 ;
		} else {
			w = t.dot( s ) / nt2 ;
		}

		x += alpha*y + w*z ;
		r = s - w*t ;

		res = r.squaredNorm() * scale;
		if( m_callback ) m_callback->trigger( k, res ) ;
		if( res < m_tol ) break ;

	}

	return res ;
}

//CGS
template < typename Mat, typename Prec, typename Traits >
template < typename RhsT, typename ResT >
typename CGS< Mat, Prec, Traits >::Scalar
CGS< Mat, Prec, Traits >::vectorSolve( const RhsT &b, ResT x ) const
{
	typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
	Vector r ;

	Scalar scale ;
	Scalar res = init( *m_A, b, x, r, scale ) ;
	if( res < m_tol ) return res ;

	const Vector r0h = r ;

	Vector u = r, p = r, q ;
	Scalar rho1, rho0 = 1, 
	       alpha, beta ;

	Vector y  ( m_A->rows() ),
			nu ( m_A->rows() ) ;

	for( unsigned k = 0 ; k < m_maxIters ; ++k )
	{
		rho1 = r0h.dot( r ) ;

		if( NumTraits< Scalar >::isZero( rho1 ) )
			break ;

		if( k > 0 )
		{
			beta = rho1 / rho0 ;
			u = r + beta *  q ;
			p = u + beta * (q + beta*p) ;
		}


		apply< false >( m_P, p, y ) ;
		nu = ( *m_A ) * y ;
		alpha = rho1 / r0h.dot( nu ) ;
		q = u - alpha*nu ;

		apply< false >( m_P, u+q, y ) ;
		nu = ( *m_A ) * y ;

		x += alpha * y ;
		r -= alpha * nu ;

		res = r.squaredNorm() * scale;
		if( m_callback ) m_callback->trigger( k, res ) ;
		if( res < m_tol ) break ;

		rho0 = rho1 ;
	}

	return res ;
}

//GMRES
template < typename Mat, typename Prec, typename Traits >
template < typename RhsT, typename ResT >
typename GMRES< Mat, Prec, Traits >::Scalar
GMRES< Mat, Prec, Traits >::vectorSolve( const RhsT &b, ResT x ) const
{
	typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
	typedef typename Traits::DynMatrix WorkMatrix ;
	typedef typename Traits::DynVector WorkVector ;
	typedef typename LocalProblemTraits< 2, Scalar >::Matrix Matrix22 ;

	Vector r ;

	Scalar scale ;
	Scalar res = init( *m_A, b, x, r, scale ) ;
	if( res < m_tol ) return res ;

	const unsigned n = b.rows() ;

	const unsigned restart = ( m_restart == 0 ) ? n : m_restart ;
	const unsigned m = std::min( restart, m_maxIters ) ;


	// Allocate working memory

	WorkMatrix H ( std::max( n, m+1 ), m + 1 ) ;		// Hessenberg
	WorkMatrix V ( n, m + 1 ) ; // Krylov subspace basis

	WorkMatrix U ( m+1, m ) ;   // Upper triangular matrix
	WorkMatrix O ( m+1, m+1 ) ; // Orthogonal matrix such that U = O*H
	Matrix22 G ;			// Givens rotation

	WorkVector g(m+1),     // O * beta * e1
			y(m+1) ;    // Res of least-square


	// Restart loop
	unsigned globalIter = 0 ;
	do
	{
		typename WorkMatrix::ColXpr v0 ( V.col(0) ) ;
		apply< false >( m_P, r, v0 ) ;
		Scalar beta = v0.norm() ;
		v0 /= beta ; // ~ r.normalized()

		O(0,0) = 1 ; // Initialize O to identity
		g(0,0) = beta ;

		unsigned k ;
		for( k = 0 ; k < m && res >= m_tol ; ++k )
		{

			// 1 - Arnoldi iteration
			typename WorkMatrix::ColXpr v ( V.col(k+1) ) ;
			r = ( *m_A ) * V.col(k) ;
			apply< false >( m_P, r, v ) ;

			H.col(k  ).head( k+1 ) = V.leftCols(k+1).transpose() * v ;
			H.row(k+1).head( k   ).setZero() ;

			v -= V.leftCols( k+1 ) * H.col( k ).head( k+1 );

			const Scalar vhn = v.norm() ;
			H(k+1, k) = vhn ;

			if( vhn > NumTraits< Scalar >::epsilon( ) ) //If vhn is zero, the algorithm shall stop at this step
				V.col(k+1) /= vhn ;

			// 2 - Least squares
			// a. Grow orthogonal matrix O and vector g = O * res0 * (1,0,0, ... )'

			O.row( k+1 ).head( k+1 ).setZero() ;     // Set last row to zero
			O.col( k+1 ).head( k+1 ).setZero() ;     // Set last col to zero
			O ( k+1, k+1 ) = 1 ;
			g ( k+1 )      = 0 ;

			// a' Store temporary, before-rotation, not yet upper-triangular new U
			U.col(k).head( k+1 ) = O.topLeftCorner( k+1, k+1 ) * H.col(k).head( k+1 ) ;
			U( k+1, k ) = vhn ;

			// b. Apply givens rotation
			G.row(0) = U.col(k).template segment< 2 >( k ).transpose() ;

			const Scalar l = G.row(0).norm() ;
			G.row(0) /= l ;

			G(1, 0) = -G( 0, 1 ) ;
			G(1, 1) =  G( 0, 0 ) ;

			O.block( k, 0, 2, k+2 ).applyOnTheLeft( G );
			g.template segment< 2 >( k ).applyOnTheLeft( G ) ;

			// c. Update upper-triagular matrix U = O * H
			U.row( k+1 ).head( k+1 ).setZero() ; // Last U line to zero
			U(k,k) = l ;

			// d. Solve triangular system
			y = U.topLeftCorner( k+1, k+1 ).template triangularView< Eigen::Upper >().solve( g.head( k+1 ) ) ;

			// 3 - Update residual
			res = g( k+1 ) * g( k+1 ) * scale ;

			//			std::cout << " ==== Iteration " << globalIter << " + " << k <<std::endl
			//			          << "H" << std::endl
			//			          << H.topLeftCorner( k+2, k+1 )<< std::endl
			//			          << "O" << std::endl
			//			          << O.topLeftCorner( k+2, k+2 )<< std::endl
			//			          << "U" << std::endl
			//			          << U.topLeftCorner( k+2, k+1 )<< std::endl
			//			          << "V" << std::endl
			//			          << V.leftCols( k+2 ) << std::endl
			//			          << "Orthogonality" << std::endl
			//			          << ( O.topLeftCorner(k+2,k+2) * O.topLeftCorner(k+2,k+2).transpose() - Matrix::Identity( k+2, k+2 ) ).squaredNorm() << std::endl
			//			          << "Equality" << std::endl
			//			          << ( O.topLeftCorner(k+2,k+2) * H.topLeftCorner(k+2,k+1) - U.topLeftCorner(k+2,k+1) ).squaredNorm()<< std::endl
			//			          << "Solve" << std::endl
			//			          << ( U.topLeftCorner( k+1, k+1 )*y - g.segment(0,k+1) ).transpose()<< std::endl
			//			          << "res" << std::endl
			//			          << g(k+1) << std::endl
			//			          << ( m_A*x - b ).norm()
			//			          << std::endl ;


			if( m_callback ) m_callback->trigger( k + globalIter, res ) ;
		}

		x += V.leftCols( k ) * y.head( k ) ;
		globalIter += m ;

		if( res < m_tol || globalIter >= m_maxIters  )
			break ;

		// Restart

		r = b ;
		r.noalias() -= ( *m_A ) * x ;
		res = r.squaredNorm() * scale ;

	} while( res >= m_tol ) ;

	return res ;
}

// TFQMR

template < typename Mat, typename Prec, typename Traits >
template < typename RhsT, typename ResT >
typename TFQMR< Mat, Prec, Traits >::Scalar
TFQMR< Mat, Prec, Traits >::vectorSolve( const RhsT &b, ResT x ) const
{
	typedef typename Traits::template MutableClone< RhsT >::Type Vector ;
	Vector r ;

	Scalar scale ;
	Scalar res = init( *m_A, b, x, r, scale ) ;
	if( res < m_tol ) return res ;

	Vector u = r, w = r ;
	const Vector& r0h = r ;

	Vector d ( m_A->rows() ),
		   v ( m_A->rows() ),
		  nu ( m_A->rows() ),
		   y ( m_A->rows() ) ;

	d.setZero() ;

	apply< false >( m_P, u, y ) ;
	v = ( *m_A ) * y ;

	Scalar tau2 = res/scale ;
	Scalar rho = tau2, rho1 ; // rho = r.dot( r0h )
	Scalar theta2 = 0, eta = 0, alpha = 0,  beta, c2 ;

	for( unsigned k = 0 ; k < m_maxIters ; ++k )
	{


		const bool odd = k&1u;
		if( !odd )
		{
			alpha = rho / r0h.dot( v ) ;
		}

		apply< false >( m_P, u, y ) ;
		nu = ( *m_A ) * y ;
		w -= alpha * nu ;
		d = y + ( theta2 / alpha ) * eta * d ;

		theta2 = w.squaredNorm() ;
		theta2 /= tau2 ;

		c2 = 1. / ( 1 + theta2 ) ;
		tau2 *= theta2 * c2 ;
		eta = c2 * alpha ;

		x += eta * d ;
		// m_A->template multiply< false >( d, r, -eta, 1 ) ;
		res = tau2 * scale ; // Approx

		if( m_callback ) m_callback->trigger( k, res ) ;
		if( res < m_tol ) break ;

		if( odd )
		{
			rho1 = w.dot( r0h ) ;
			beta = rho1 / rho ;

			u = w + beta * u ;

			v = nu + beta * v ;
			apply< false >( m_P, u, y ) ;
			v *= beta ;
			v.noalias() += ( *m_A ) * y ;

			rho = rho1 ;
		} else {
			u -= alpha * v ;
		}

		//		rho0 = rho1 ;
	}

	return res ;
}

} //namespace solvers

} //namespace krylov

} //namespace bogus

#endif
