/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BANDMATRIXFWD_HH
#define BANDMATRIXFWD_HH

#include "Definitions.hh"

namespace strandsim {

template <typename ScalarT, IndexType kl, IndexType ku>
class BandMatrix;

typedef BandMatrix<Scalar, 10, 10> JacobianMatrixType;
typedef BandMatrix<Scalar, 15, 15> RestJacobianMatrixType;
typedef BandMatrix<Scalar, 1, 1> TriDiagonalMatrixType;

template <typename ScalarT, int kl>
class SymmetricBandMatrixSolver;
typedef SymmetricBandMatrixSolver<Scalar, 10> JacobianSolver;

}  // namespace strandsim

#endif  // BANDMATRIXFWD_HH
