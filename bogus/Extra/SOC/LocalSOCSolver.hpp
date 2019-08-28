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


#ifndef BOGUS_LOCAL_SOC_SOLVER_HPP
#define BOGUS_LOCAL_SOC_SOLVER_HPP

#include "../SecondOrder.fwd.hpp"

namespace bogus {


template< DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV,
		   local_soc_solver::Strategy Strat = local_soc_solver::Hybrid  >
struct LocalSOCSolver
{

  typedef LocalProblemTraits< Dimension, Scalar > Traits ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  static Scalar solve(
		  const typename Traits::Matrix &A,
		  const typename Traits::Vector &b,
		  typename Traits::Vector &x,
		  const Scalar mu, const Scalar tol,
		  const Scalar scaling = 1
		  ) ;

} ;

}

#endif
