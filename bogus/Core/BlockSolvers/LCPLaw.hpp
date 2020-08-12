/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_LCPLAW_HPP
#define BOGUS_LCPLAW_HPP

#include <cmath>

namespace bogus {

//! LCP local solver that can be used within GaussSeidel and ProjectedGradient
//! solvers
/*!
  For demonstration purposes. Since blocks are 1x1, other libraries are probably
  more suited. \tparam Scalar the scalar type \tparam Dimension the dimension of
  the blocks of the global matrix
  */
template <typename Scalar>
class LCPLaw {
 public:
  enum { dimension = 1 };

  typedef LocalProblemTraits<dimension, Scalar> Traits;

  //! Constructor
  LCPLaw() {}

  //! \return \f$ \vert fb( x, y ) \vert^2_2 \f$, where fb is the scalar
  //! Fischer-Burmeister function
  Scalar eval(const unsigned problemIndex, const typename Traits::Vector &x,
              const typename Traits::Vector &y) const;

  //! Solves the local problem
  /*!
          \f$ 0 \leq x \perp a x + b \geq 0 \f$
  */
  bool solveLocal(const unsigned problemIndex, const typename Traits::Matrix &A,
                  const typename Traits::Vector &b, typename Traits::Vector &x,
                  const Scalar scaling) const;

  //! Projects x on \f$ R^+ \f$
  void projectOnConstraint(const unsigned problemIndex,
                           typename Traits::Vector &x) const;

  //! This NSLaw is always associated, so dualityCOV is null.
  template <typename Segment>
  void dualityCOV(const unsigned, const Segment &,
                  typename Traits::Vector &s) const {
    s->setZero();
  }
};

}  // namespace bogus

#endif
