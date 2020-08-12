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

#ifndef BOGUS_SOCLAW_HPP
#define BOGUS_SOCLAW_HPP

#include <vector>

#include "FischerBurmeister.hpp"
#include "LocalSOCSolver.hpp"

namespace bogus {

//! Non-smooth laws based on Second Order Cone complementarity. To be used
//! within as the first argument to GaussSeidel::solve().
/*!
        \tparam Dimension the dimension of the local problem. Specializations
   exist form dimension 2 and 3. \tparam Scalar the scalar type \tparam
   DeSaxceCOV Whether to perform the \cite DSF98 change of variable when solving
   the local problem. Should be \c true for modeling Coulomb friction, or \c
   false for standard SOC complementarity. \sa solveLocal() \tparam Strat
   local_soc_solver::Strategy for solving the local problems. Unavailable for
   dimensions other than 2 and 3.
  */
template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV,
          local_soc_solver::Strategy Strat>
class SOCLaw {
 public:
  typedef LocalProblemTraits<Dimension, Scalar> Traits;
  enum { dimension = Dimension };

  //! Constructor
  /*!
    \param n the size of the global problem ( number of contacts )
    \param mu array containing the apertures of each second order cone (
    friction coefficients )
    */
  SOCLaw(const unsigned n, const double *mu);

  //! \return \f$ \vert fb( mu, x, y ) \vert^2_2 \f$, where fb is the SOC
  //! Fischer-Burmeister function
  Scalar eval(const unsigned problemIndex, const typename Traits::Vector &x,
              const typename Traits::Vector &y, const double scaling) const {
    typedef FischerBurmeister<Traits::dimension, typename Traits::Scalar,
                              DeSaxceCOV>
        FBFunction;

    if (m_mu[problemIndex] < 0) return y.squaredNorm();

    typename Traits::Vector fb(x.rows());
    FBFunction::compute(m_mu[problemIndex], x * scaling, y, fb);

    return fb.squaredNorm();
  }

  //! Solves the local problem
  /*!
    \f[
          \left\{
            \begin{array}{rcl}
                  y &=& \mathrm{ <DS> } \left( A x + b \right ) \\
                  K_{ \frac 1 \mu } \ni y & \perp & x \in K_{ \mu }
            \end{array}
          \right.
    \f]
    where \f$ \mu \f$ is \c m_mu[\p problemIndex] and \c <DS> is the optional De
    Saxce change of variable.

    That is, if \p DeSaxceCOV is false then \c <DS> is the identity function,
    otherwise \f[ \mathrm{ <DS> }( x ) := x + \mu \vert x_T \vert \left( 1, 0,
    ... \right)^{\top} \f] \param scaling Used as a scaling factor for \p x when
    calculating the error function
  */
  bool solveLocal(const unsigned problemIndex, const typename Traits::Matrix &A,
                  const typename Traits::Vector &b, typename Traits::Vector &x,
                  const Scalar scaling) const;

  //! Projects x on \f$ K_{ \mu } \f$
  void projectOnConstraint(const unsigned problemIndex,
                           typename Traits::Vector &x) const;

  //! Computes the change of variable \p s(y) so that (x, y+s(y)) obeys an
  //! associated law.
  /*! ie \f$ y + s(y) \in - N_C(x) \f$.
   *  Here C = K_{mu}, and if \tparam DeSaxceCOV is true, \f$ s_i(y) =  mu_i
   * |y_{i,T}| n \f$
   */
  template <typename Segment>
  void dualityCOV(const unsigned problemIndex, const Segment &y,
                  typename Traits::Vector &s) const {
    if (DeSaxceCOV) {
      Traits::np(s) = m_mu[problemIndex] * Traits::tp(y).norm();
      Traits::tp(s).setZero();
    } else
      s.setZero();
  }

 private:
  const double *m_mu;
  const unsigned m_n;
  Scalar m_localTol;
};

//! Predefined non-smooth law for 2D Coulomb friction
typedef SOCLaw<2, double, true> Coulomb2D;
//! Predefined non-smooth law for 3D Coulomb friction
typedef SOCLaw<3, double, true> Coulomb3D;
//! Predefined non-smooth law for 2D SOC complementarity
typedef SOCLaw<2, double, false> SOC2D;
//! Predefined non-smooth law for 3D SOC complementarity
typedef SOCLaw<3, double, false> SOC3D;

}  // namespace bogus

#endif  // SOCLAW_HPP
