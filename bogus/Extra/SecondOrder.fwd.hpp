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

#ifndef BOGUS_SECOND_ORDER_FWD_HPP
#define BOGUS_SECOND_ORDER_FWD_HPP

#include "../Core/Block.fwd.hpp"

namespace bogus {

//! Configuration properties of local Second Order Cone solver
namespace local_soc_solver {
//! Strategy to be used by the local SOC solver.
/*! Note that some strategies may be unavailable for some loval problem types,
 in which case the solver will revert to the PureNewton strategy */
enum Strategy {
  PureNewton  //!< Newton algorithm on the SOC FischerBurmeister function. \sa
              //!< NonSmoothNewton
#ifndef BOGUS_WITHOUT_EIGEN
  ,
  PureEnumerative  //!< Enumerative algorithm, such as describer in Appendix B
                   //!< of \cite DBB11
  ,
  Hybrid  //!< Newton algorithm, then Enumerative as failsafe
  ,
  RevHybrid  //!< Enumerative algorithm, then Newton to refine the solution
#endif
};
}  // namespace local_soc_solver

template <DenseIndexType Dimension, typename Scalar>
struct LocalProblemTraits;

template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV,
#ifndef BOGUS_WITHOUT_EIGEN
          local_soc_solver::Strategy Strat = local_soc_solver::RevHybrid>
#else
          local_soc_solver::Strategy Strat = local_soc_solver::PureEnumerative>
#endif
class SOCLaw;

template <DenseIndexType Dimension, typename Scalar,
          local_soc_solver::Strategy Strat = local_soc_solver::PureEnumerative>
class HBLaw;

}  // namespace bogus

#endif
