/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2016 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_PYRAMIDLAW_HPP
#define BOGUS_PYRAMIDLAW_HPP

#include <cmath>

namespace bogus {

//! Experimental and incomplete pyramidal local solver that can be used within
//! GaussSeidel and ProjectedGradient solvers
/*!
  \warning solveLocal() is currently implemented in a highly inefficient way.
  This NSLaw may be used to compare the results with SOCLaw, not to compare
  performance. \warning projectOnConstraint() only handles dimensions 2 and 3

  For demonstration purposes. Since blocks are 1x1, other libraries are probably
  more suited. \tparam Scalar the scalar type \tparam Dimension the dimension of
  the blocks of the global matrix
  */
template <DenseIndexType Dimension, typename Scalar, bool DeSaxceCOV>
class PyramidLaw {
 public:
  typedef LocalProblemTraits<Dimension, Scalar> Traits;
  enum { dimension = Dimension };

  //! Constructor
  /*!
    \param n the size of the global problem ( number of contacts )
    \param mu array containing the apertures of each second order cone (
    friction coefficients )
    */
  PyramidLaw(const unsigned n, const double *mu);

  //! \return \f$ \vert fb( mu, x, y ) \vert^2_2 \f$, where fb is the SOC
  //! Fischer-Burmeister function
  Scalar eval(const unsigned problemIndex, const typename Traits::Vector &x,
              const typename Traits::Vector &y) const;

  //! Solve for \f$ y(x) \in - N_C (x)\f$ with y(x) = Ax + b
  bool solveLocal(const unsigned problemIndex, const typename Traits::Matrix &A,
                  const typename Traits::Vector &b, typename Traits::Vector &x,
                  const Scalar scaling) const;

  //! Projects x on \f$ R^+ \f$
  void projectOnConstraint(const unsigned problemIndex,
                           typename Traits::Vector &x) const;

  //! Computes the change of variable \p s(y) so that (x, y+s(y)) obeys an
  //! associated law.
  /*! ie \f$ y + s(y) \in - N_C(x) \f$.
   *  Here C = K_{mu}, and if \tparam DeSaxceCOV is true, \f$ s_i(y) =  mu_i
   * \sum |y_{i,T,k}| n \f$
   */
  template <typename Segment>
  void dualityCOV(const unsigned problemIndex, const Segment &y,
                  typename Traits::Vector &s) const {
    if (DeSaxceCOV) {
      Traits::np(s) = m_mu[problemIndex] * Traits::tp(y).template lpNorm<1>();
      Traits::tp(s).setZero();
    } else
      s.setZero();
  }

 private:
  const double *m_mu;
  const unsigned m_n;
};

}  // namespace bogus

#endif
