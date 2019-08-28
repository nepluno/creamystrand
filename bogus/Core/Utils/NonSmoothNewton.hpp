/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BOGUS_NS_NETWON_HPP
#define BOGUS_NS_NETWON_HPP

namespace bogus {

//! Dense, naive Newton implementation
/*!
  \tparam NSFunction The type of the function for which to find the zeros.
  Should define \c compute() and \c computeJacobian() methods.
  Does not have to be smooth or convex, but that helps.
  \sa FischerBurmeister
  */
template < typename NSFunction >
class NonSmoothNewton
{
public:
  typedef typename NSFunction::Traits Traits ;
  typedef typename Traits::Scalar Scalar ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  NonSmoothNewton( const NSFunction &func, Scalar tol )
  : m_func( func ), m_tol( tol ), m_maxIters( 20 )
  {}

  //! Tries to find \p x such that m_func( x ) = 0
  /*! Or, more precisely, such that \f[ \Phi(x) = \frac 1 2 \vert f(x) \vert_2^2 < tol \f]
	\param x Used both as an initial guess and to return the approximate solution
	\return \f$ \Phi(x) \f$
  */
  Scalar solve ( Vector &x ) const ;

  //! Sets the maximum number of iterations
  void setMaxIters( unsigned maxIters ) { m_maxIters = maxIters ; }

  //! Sets the solver tolerance
  void setTol( Scalar tol ) { m_tol = tol ; }

private:
  const NSFunction& m_func ;
  Scalar m_tol ;
  unsigned m_maxIters ;

} ;


}

#endif
