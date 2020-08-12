/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_SOLVERS_FWD_HPP
#define BOGUS_BLOCK_SOLVERS_FWD_HPP

namespace bogus {

//! Configuration properties for Krylov solvers
namespace krylov {
enum Method {
  CG,        //!< Conjugate Gradient. \sa krylov::solvers::CG
  BiCG,      //!< BiConjugate Gradient \sa krylov::solvers::BiCG
  BiCGSTAB,  //!< BiConjugate Gradient Stabilized \sa krylov::solvers::BiCGSTAB
  CGS,       //!< Conjugate Gradient Squared \sa krylov::solvers::CGS
  GMRES,     //!< Generalized Minimal Residual \sa krylov::solvers::GMRES
  TFQMR  //!< Tranpose-free Quasi Minimal Residual \sa krylov::solvers::TFQMR
};

}  // namespace krylov

//! Options for ProjectedGradient solvers
namespace projected_gradient {
//! Variants of Projected Gradient algorithm
enum Variant {
  //! Standard projected gradient
  Standard,
  //! Projected gradient descent
  Descent,
  //! Projected gradient with conjugation of search direction
  Conjugated,
  //! Accelerated Projected Gradient Descent based on \cite Nesterov1983 and
  //! developed in \cite Heyn13
  APGD,
  //! Spectral Projected Gradient, loosely adapted from \cite Tasora13
  SPG
};
}  // namespace projected_gradient

template <typename MatrixType>
struct ProblemTraits;

template <typename Derived, typename BlockMatrixType>
class ConstrainedSolverBase;

template <typename BlockMatrixType>
class GaussSeidel;

template <typename BlockMatrixType, typename DiagonalType, bool PrecomputeDMt>
class ProductGaussSeidel;

template <typename BlockMatrixType>
class ProjectedGradient;

template <typename BlockMatrixType>
class ADMM;

template <typename BlockMatrixType>
class DualAMA;

template <typename MatrixType>
class TrivialPreconditioner;

template <typename BlockMatrixType, template <typename BlockMatrixT>
                                    class PreconditionerType =
                                        TrivialPreconditioner>
class Krylov;

}  // namespace bogus

#endif
