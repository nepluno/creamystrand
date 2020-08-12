/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_SOLVER_BASE_HPP
#define BOGUS_BLOCK_SOLVER_BASE_HPP

#include "../Block.fwd.hpp"
#include "../BlockSolvers.fwd.hpp"

#ifndef BOGUS_WITHOUT_EIGEN
#include "../Eigen/EigenBlockContainers.hpp"
#include "../Eigen/EigenProblemTraits.hpp"
#endif

#include "../Utils/Signal.hpp"

namespace bogus {

//! Base class for solvers that operate on BlockMatrixBase matrices
template <typename BlockMatrixType>
class BlockSolverBase {
 public:
  typedef BlockMatrixTraits<BlockMatrixType> BlockTraits;
  typedef typename BlockTraits::Scalar Scalar;
  typedef ProblemTraits<Scalar> GlobalProblemTraits;
  typedef Signal<unsigned, Scalar> CallBackType;

  virtual ~BlockSolverBase() {}

  //! For iterative solvers: sets the maximum number of iterations
  void setMaxIters(unsigned maxIters) { m_maxIters = maxIters; }
  unsigned maxIters() const { return m_maxIters; }

  //! For iterative solvers: sets the solver tolerance
  void setTol(Scalar tol) { m_tol = tol; }
  Scalar tol() const { return m_tol; }

  //! Callback hook; will be triggered every N iterations, depending on the
  //! solver
  /*! Useful to monitor the convergence of the solver. Can be connected to a
     function that takes an \c unsigned and a \c Scalar as parameters. The first
     argument will be the current iteration number, and the second the current
     residual. \sa Signal< unsigned, Scalar > */
  CallBackType &callback() { return m_callback; }
  const CallBackType &callback() const { return m_callback; }

  const BlockObjectBase<BlockMatrixType> &matrix() const { return *m_matrix; }

 protected:
  BlockSolverBase(const BlockObjectBase<BlockMatrixType> *matrix =
                      BOGUS_NULL_PTR(const BlockObjectBase<BlockMatrixType>),
                  unsigned maxIters = 0, Scalar tol = 0)
      : m_matrix(matrix), m_maxIters(maxIters), m_tol(tol) {}

  //! Pointer to the matrix of the system
  const BlockObjectBase<BlockMatrixType> *m_matrix;

  //! See setMaxIters()
  unsigned m_maxIters;
  //! See setTol()
  Scalar m_tol;

  CallBackType m_callback;
};

}  // namespace bogus

#endif
