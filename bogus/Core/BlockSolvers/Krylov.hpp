/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_KRYLOV_HPP
#define BOGUS_KRYLOV_HPP

#include <vector>

#include "BlockSolverBase.hpp"
#include "KrylovMethods.hpp"
#include "Preconditioners.hpp"

namespace bogus {

//! Preconditionned Krylov Solvers
/*!
  \tparam BlockMatrixT The type of system matrix, which should be a subclass of
  BlockObjectBase \tparam PreconditionerType The preconditioner type. It should
  accept BlockMatrixT as a template parameter. The default value,
  TrivialPreconditioner, means that no preconditioning will be done. \sa
  TrivialPreconditioner, DiagonalPreconditioner, DiagonalLUPreconditioner,
  DiagonalLDLTPreconditioner \sa krylov
  */

template <typename BlockMatrixType,
          template <typename BlockMatrixT> class PreconditionerType>
class Krylov : public BlockSolverBase<BlockMatrixType> {
 public:
  typedef BlockSolverBase<BlockMatrixType> Base;
  typedef PreconditionerType<BlockObjectBase<BlockMatrixType> >
      PreconditionerImplType;

  typedef typename Base::GlobalProblemTraits GlobalProblemTraits;
  typedef typename GlobalProblemTraits::Scalar Scalar;

  //! Default constructor -- you will have to call setMatrix() before using any
  //! of the solve() functions
  Krylov();
  //! Constructor with the system matrix -- initializes preconditioner
  explicit Krylov(const BlockObjectBase<BlockMatrixType> &matrix);

  //! Sets the system matrix and initializes the preconditioner
  Krylov &setMatrix(const BlockObjectBase<BlockMatrixType> &matrix);

  // For each value of the krylov::Method enum, create a MethodNameType typedef
  // and the asMethodName() -> MethodNameType and solve_MethodName -> Scalar
  // methods
#define BOGUS_PROCESS_KRYLOV_METHOD(MethodName)                               \
  typedef krylov::solvers::MethodName<                                        \
      BlockMatrixType, PreconditionerType<BlockObjectBase<BlockMatrixType> >, \
      GlobalProblemTraits>                                                    \
      MethodName##Type;                                                       \
                                                                              \
  MethodName##Type as##MethodName() const {                                   \
    return MethodName##Type(Base::m_matrix->derived(), Base::m_maxIters,      \
                            Base::m_tol, &m_preconditioner,                   \
                            &this->m_callback);                               \
  }                                                                           \
                                                                              \
  template <typename RhsT, typename ResT>                                     \
  Scalar solve_##MethodName(const RhsT &b, ResT &x) const {                   \
    return as##MethodName().solve(b, x);                                      \
  }

  BOGUS_KRYLOV_METHODS
#undef BOGUS_PROCESS_KRYLOV_METHOD

  //! Solve function that takes the method to use as an argument
  template <typename RhsT, typename ResT>
  Scalar solve(const RhsT &b, ResT &x,
               krylov::Method method = krylov::CG) const;

  const PreconditionerImplType &preconditioner() const {
    return m_preconditioner;
  }
  PreconditionerImplType &preconditioner() { return m_preconditioner; }

 protected:
  PreconditionerImplType m_preconditioner;
};

}  // namespace bogus

#endif
