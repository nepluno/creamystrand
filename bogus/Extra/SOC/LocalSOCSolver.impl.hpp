/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * So-bogus is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * So-bogus is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with So-bogus.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef BOGUS_LOCAL_SOC_SOLVER_IMPL_HPP
#define BOGUS_LOCAL_SOC_SOLVER_IMPL_HPP

#include "LocalSOCSolver.hpp"

#include "FischerBurmeister.impl.hpp"

#include "../../Core/Utils/LinearSolverBase.hpp"
#include "../../Core/Utils/NonSmoothNewton.impl.hpp"
#include "../../Core/Utils/Polynomial.impl.hpp"

namespace bogus {

// No analytic solution in the general case
template < DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV >
struct AnalyticLocalSOCSolver
{
	typedef LocalProblemTraits< Dimension, Scalar > Traits ;
	typedef typename Traits::Vector Vector ;
	typedef typename Traits::Matrix Matrix ;

	static bool solveOrthogonality(
			const typename Traits::Matrix &,
			const typename Traits::Vector &,
			typename Traits::Vector &,
			const Scalar
			)
	{
		return false ;
	}
} ;

// Specialization for Coulomb 3D Friction
template< typename Scalar >
struct AnalyticLocalSOCSolver< 3u, Scalar, true >
{
	enum { Dimension = 3 } ;
	typedef LocalProblemTraits< Dimension, Scalar > Traits ;
	typedef typename Traits::Vector Vector ;
	typedef typename Traits::Matrix Matrix ;

	static bool solveOrthogonality(
			const typename Traits::Matrix &W,
			const typename Traits::Vector &b,
			typename Traits::Vector &r,
			const Scalar mu
			)
	{
		// see [Daviet et al 2011], Appendix B.1

		typedef typename LocalProblemTraits< 2, Scalar >::Vector Vec2 ;
		typedef typename LocalProblemTraits< 2, Scalar >::Matrix Mat2 ;

		const Scalar wN = W(0,0) ;
		if( wN < NumTraits< Scalar >::epsilon() )
			return false ; // Could we do something better ?

		const Vec2 wT = W.template block< 2, 1 >( 1, 0 ) ;
		const Mat2 WT = W.template block< 2, 2 >( 1, 1 ) ;

		const Scalar bN = Traits::np( b );
		const Vec2 bT = Traits::tp( b ) ;

		const Vec2 wT_wN = wT/wN ;
		const Mat2 Wbar = WT - wT_wN * wT.transpose() ;
		const Vec2 bbar = bT/bN - wT_wN ;

		const Scalar A = Wbar.trace() - wT.dot( bbar ) ;
		const Vec2   B ( Wbar(1,1)*bbar[0] - Wbar(1,0)*bbar[1],
				Wbar(0,0)*bbar[1] - Wbar(0,1)*bbar[0] ) ;
		const Scalar C = Wbar.determinant() - wT.dot( B ) ;
		const Scalar D = wN*wN / ( mu*mu ) ;

		const Scalar coeffs[4] = {
			C*C - D * B.squaredNorm(),
			2*( C*A - D * bbar.dot( B ) ),
			2*C + A*A - D * bbar.squaredNorm(),
			2*A
		} ;

		Scalar roots[4] ;
		const unsigned nRoots =
				polynomial::getRealRoots( coeffs, roots, polynomial::StrictlyPositiveRoots ) ;

		if( 0 == nRoots ) return false ;

		Scalar alpha = roots[0] ;
		//Get the minimal one, this is as good an heuristic as any
		for ( unsigned i = 1 ; i != nRoots ; ++ i )
		{
			if( roots[i] < alpha ) alpha = roots[i] ;
		}

		//	 std::cout << "Found " << alpha << std::endl ;

		const Mat2 M = Wbar + alpha * Mat2::Identity() ;
		const typename MatrixTraits< Mat2 >::LUType lu (M) ;
		Traits::tp( r ) = lu.solve( -bN * bbar ) ;
		Traits::np( r ) = Traits::tp( r ).norm() / mu ;

		return true ;
	}

} ;


// Specialization for Coulomb 2D Friction
template< typename Scalar >
struct AnalyticLocalSOCSolver< 2u, Scalar, true >
{
	enum { Dimension = 2 } ;
	typedef LocalProblemTraits< Dimension, Scalar > Traits ;
	typedef typename Traits::Vector Vector ;
	typedef typename Traits::Matrix Matrix ;

	static bool solveOrthogonality(
			const typename Traits::Matrix &W,
			const typename Traits::Vector &b,
			typename Traits::Vector &r,
			const Scalar mu
			)
	{
		bool found = false ;

		// Case rT = mu * rN
		const Scalar w_p = W(0,0) + mu * W(0,1) ;
		if( w_p > NumTraits< Scalar >::epsilon( ) )
		{
			r[ 0 ] = -b[0] / w_p ;
			r[ 1 ] = mu * r[0] ;
			const Scalar minus_alpha = W(1,0)/mu + W(1,1) + b(1)/(mu*r[0]) ;
			found = 0 > minus_alpha ;
		}

		if( !found )
		{
			// Case rT = - mu * rN
			const Scalar w_m = W(0,0) - mu * W(0,1) ;
			if( w_m > NumTraits< Scalar >::epsilon( ) )
			{
				r[ 0 ] = -b[0] / w_m ;
				r[ 1 ] = -mu * r[0] ;
				const Scalar alpha = W(1,0)/mu - W(1,1) + b(1)/(mu*r[0]) ;
				found = 0 < alpha ;
			}
		}

		return found ;
	}
} ;

// Specialization for 3D SOC complementarity
template< typename Scalar >
struct AnalyticLocalSOCSolver< 3u, Scalar, false >
{
	enum { Dimension = 3 } ;
	typedef LocalProblemTraits< Dimension, Scalar > Traits ;
	typedef typename Traits::Vector Vector ;
	typedef typename Traits::Matrix Matrix ;

	static bool solveOrthogonality(
			const typename Traits::Matrix &W,
			const typename Traits::Vector &b,
			typename Traits::Vector &r,
			const Scalar mu
			)
	{
		// see doc/sage/poySOC.sage

		const Scalar A = (b[1]*W(1,2) - b[2]*W(1,1))*mu*mu + b[0]*W(0,2) - b[2]*W(0,0) ;
		const Scalar B = (b[0]*W(1,2) + b[1]*W(0,2) - 2*b[2]*W(0,1))*mu ;
		const Scalar C = 2*( (b[1]*W(2,2) - b[2]*W(1,2))*mu*mu - b[0]*W(0,1) + b[1]*W(0,0) );
		const Scalar D = 2*( b[0]*(W(1,1) - W(2,2)) - b[1]*W(0,1) + b[2]*W(0,2))*mu ;
		const Scalar E = -6*(b[0]*W(1,2) - b[1]*W(0,2))*mu ;

		Scalar coeffs[5] ;
		coeffs[0] =  A + B ;
		coeffs[1] =  C - D ;
		coeffs[2] =  E ;
		coeffs[3] =  C + D ;
		coeffs[4] =  B - A ;

		Scalar roots[4] ;
		const unsigned nRoots =
				bogus::polynomial::getRealRoots( coeffs, roots, bogus::polynomial::AllRoots ) ;

		bool found = false ;

		for ( unsigned i = 0 ; i != nRoots ; ++ i )
		{
			Scalar t = roots[i] ;

			const Scalar CT = ( 1 - t*t ) / ( 1 + t*t ) ;
			const Scalar ST = 2*t / ( 1 + t*t ) ;

			const typename Traits::Vector dir ( 1, mu*CT, mu*ST ) ;

			Scalar rN ;

			// Chose the better conditionned way to go back to r
			const Scalar den0 = ( mu * W.col(1) + CT * W.col( 0 )).dot( dir ) ;
			const Scalar den1 = ( mu * W.col(2) + ST * W.col( 0 )).dot( dir ) ;

			if( std::fabs( den0 ) > std::fabs( den1 ) )
			{
				if( bogus::NumTraits< Scalar >::isZero( den0 ) )
				{
					continue ;
				}
				rN = -(CT*b[0] + b[1]*mu)/den0 ;
			} else {
				if( bogus::NumTraits< Scalar >::isZero( den1 ) )
				{
					continue ;
				}
				rN = -(ST*b[0] + b[2]*mu)/den1 ;
			}

			if( rN <= 0 )
				continue ;

			r = rN * dir ;
			const Scalar uN = W.col(0).dot(r) + b[0] ;

			if( uN > 0 )
			{
				found = true ;
				break ;
			}
		}

		return found ;
	}

} ;

// Specialization for 2D SOC complementarity
template< typename Scalar >
struct AnalyticLocalSOCSolver< 2u, Scalar, false >
{
	enum { Dimension = 2 } ;
	typedef LocalProblemTraits< Dimension, Scalar > Traits ;
	typedef typename Traits::Vector Vector ;
	typedef typename Traits::Matrix Matrix ;

	static bool solveOrthogonality(
			const typename Traits::Matrix &W,
			const typename Traits::Vector &b,
			typename Traits::Vector &r,
			const Scalar mu
			)
	{
		bool found = false ;

		// Case rT = mu * rN

		for( int s = 1 ; !found && s > -2 ; s -= 2 )
		{
			const Scalar w = W(0,0) + 2 * mu * s * W(0,1) + mu * mu * W( 1, 1 ) ;
			if( !NumTraits< Scalar >::isZero( w ) )
			{
				r[ 0 ] = - ( b[0] + mu * s * b[1] ) / w ;
				if( r[0] > 0 )
				{
					r[ 1 ] = mu * s * r[0] ;

					const Scalar uN = W(0,0)*r[0] + W(1,0)*r[1] + b[0] ;
					found = 0 < uN ;
				}
			}
		}

		return found ;
	}
} ;


template< DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV, local_soc_solver::Strategy Strat >
Scalar LocalSOCSolver< Dimension, Scalar, DeSaxceCOV, Strat >::solve(
		const typename Traits::Matrix &A,
		const typename Traits::Vector &b,
		typename Traits::Vector &x,
		const Scalar mu, const Scalar tol, const Scalar scaling
		)
{
        if(mu < 0 )
        {
          x = A.fullPivHouseholderQr().solve( -b ) ;
          return (A*x + b).squaredNorm() ;
        }

	// see [Daviet et al 2011], Appendix B.2

	// Newton solver
	typedef FischerBurmeister< Dimension, Scalar, DeSaxceCOV > FBFunc ;
	typedef typename Traits::LUType LUType ;

	FBFunc fb( mu, A, b, scaling ) ;
	NonSmoothNewton< FBFunc > nsNewton( fb, tol )  ;

	if( Strat == local_soc_solver::PureNewton )
	{
		return nsNewton.solve( x ) ;
	}

	if( Traits::np(b) >= ( DeSaxceCOV ? 0 : mu * Traits::tp(b).norm() ) )
	{
		// Take-off case ( -b in normal cone of constraint )
		x.setZero() ;
		return 0. ;
	}
	if( NumTraits< Scalar >::isZero( mu ) )
	{
		//Frictionless case
		if( A(0,0) < NumTraits< Scalar >::epsilon() )
		{
			// Degenerate problem
			x.setZero() ;
			return b[0]*b[0] ;
		} else {
			Traits::tp( x ).setZero() ;
			Traits::np( x ) = - Traits::np( b ) / A(0,0);
			return 0. ;
		}
	}

	double res = 0. ;

	if( Strat == local_soc_solver::Hybrid )
	{
		res = nsNewton.solve( x ) ;
		if( res < tol ) return res ;
	}

	// Continuing enumerative fallback

	Vector x0 = x ;
	LUType( A ).solve( -b, x ) ;
	if( mu * Traits::np( x ) >= Traits::tp( x ).norm() )
	{
		// Sticking case
		return 0. ;
	}

	// Sliding case
	if( ! AnalyticLocalSOCSolver< Dimension, Scalar, DeSaxceCOV >::solveOrthogonality( A, b, x, mu ) )
	{
		x = x0 ;
	}

	// Refinement of final solution
	if( Strat == local_soc_solver::RevHybrid ) {
		res = nsNewton.solve( x ) ;
	} else if( Strat == local_soc_solver::Hybrid  ) {
		const double refinedRes = nsNewton.solve( x ) ;
		if( refinedRes <= res )
			return refinedRes ;

		//This can happen if the quartic solver returned a very bad value, like an
		// unreastically huge alpha
		x = x0 ;
	}

	return res ;
}


}

#endif
