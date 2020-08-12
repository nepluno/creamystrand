/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_KRYLOV_METHODS_HPP
#define BOGUS_KRYLOV_METHODS_HPP

#include "../BlockSolvers.fwd.hpp"
#include "../Utils/CppTools.hpp"
#include "../Utils/LinearSolverBase.hpp"
#include "../Utils/NumTraits.hpp"
#include "../Utils/Signal.hpp"

#define BOGUS_KRYLOV_METHODS            \
  BOGUS_PROCESS_KRYLOV_METHOD(CG)       \
  BOGUS_PROCESS_KRYLOV_METHOD(BiCG)     \
  BOGUS_PROCESS_KRYLOV_METHOD(BiCGSTAB) \
  BOGUS_PROCESS_KRYLOV_METHOD(CGS)      \
  BOGUS_PROCESS_KRYLOV_METHOD(GMRES)    \
  BOGUS_PROCESS_KRYLOV_METHOD(TFQMR)

namespace bogus {

//! Krylov methods implementation
namespace krylov {

/*!
  Base class for Krylov solvers implementations
        \note These implementations are not able to process multiple rhs
  simultaneously ; instead they will be solved for sequentially ( or in parallel
  if parallelizeRhs( true ) has been called )
*/
template <template <typename, typename, typename> class Method, typename Matrix,
          typename Preconditioner, typename Traits>
struct KrylovSolverBase
    : public LinearSolverBase<Method<Matrix, Preconditioner, Traits> > {
  typedef Method<Matrix, Preconditioner, Traits> Derived;
  typedef LinearSolverBase<Derived> Base;
  typedef typename Traits::Scalar Scalar;
  typedef Signal<unsigned, Scalar> SignalType;

  const Matrix *m_A;
  const Preconditioner *m_P;
  const SignalType *m_callback;

  Scalar m_tol;
  unsigned m_maxIters;

  KrylovSolverBase(const Matrix &A, unsigned maxIters, Scalar tol,
                   const Preconditioner *P, const SignalType *callback)
      : m_A(&A),
        m_P(P),
        m_callback(callback),
        m_tol(tol),
        m_maxIters(maxIters),
        m_parallelizeRhs(false),
        m_enableResCaching(false) {}

  KrylovSolverBase()
      : m_A(BOGUS_NULL_PTR(const Matrix)),
        m_P(BOGUS_NULL_PTR(const Preconditioner)),
        m_callback(BOGUS_NULL_PTR(const SignalType)),
        m_tol(0),
        m_maxIters(0),
        m_parallelizeRhs(false),
        m_enableResCaching(false) {}

  //! Returns the solution \b x of the linear system \b M \c * \b x \c = \c rhs
  template <typename RhsT>
  typename LinearSolverTraits<Derived>::template Result<RhsT>::Type solve(
      const RhsT &rhs) const {
    typedef typename LinearSolverTraits<Derived>::template Result<RhsT>::Type
        ReturnType;
    static ReturnType s_cachedRes;

    if (!m_enableResCaching) {
      ReturnType x(m_A->rows(), rhs.cols());
      x.setZero();
      Base::solve(rhs, x);
      return x;
    }

    if (s_cachedRes.rows() != m_A->rows() || s_cachedRes.cols() != rhs.cols()) {
      s_cachedRes.resize(m_A->rows(), rhs.cols());
      s_cachedRes.setZero();
    }

    Base::solve(rhs, s_cachedRes);
    return s_cachedRes;
  }

  //! Returns the solution \b x of the linear system \b M \c * \b x \c = \c rhs
  template <typename ResT, typename RhsT>
  Scalar solve(const RhsT &rhs, ResT &x) const {
    Scalar res = 0;
#ifndef BOGUS_DONT_PARALLELIZE
#pragma omp parallel for reduction(+ : res) if (m_parallelizeRhs)
#endif
    for (std::ptrdiff_t c = 0; c < (std::ptrdiff_t)rhs.cols(); ++c) {
      res += Base::derived().vectorSolve(rhs.col(c), x.col(c));
    }

    return res;
  }

  //! Whether to process multiple rhs in parallel
  Derived &parallelizeRhs(bool parallelize = true) {
    m_parallelizeRhs = parallelize;
    return Base::derived();
  }

  //! Whether to enable caching of solve(rhs) result for warmstarting purposes
  /*! \warning Not thread safe for parallel calls to solve(rhs) */
  Derived &enableResCaching(bool doCache = true) {
    m_enableResCaching = doCache;
    return Base::derived();
  }

 protected:
  bool m_parallelizeRhs;
  bool m_enableResCaching;
};

//! Namespace containing the implementations of the methods listed in Method, as
//! specializations of KrylovSolverBase
namespace solvers {

// Useful macros for common declarations

#define BOGUS_MAKE_KRYLOV_SOLVER_TYPEDEFS(MethodName)                   \
  typedef KrylovSolverBase<solvers::MethodName, Matrix, Preconditioner, \
                           Traits>                                      \
      Base;                                                             \
  typedef typename Traits::Scalar Scalar;                               \
                                                                        \
  using Base::m_A;                                                      \
  using Base::m_P;                                                      \
  using Base::m_maxIters;                                               \
  using Base::m_tol;                                                    \
  using Base::m_callback;

#define BOGUS_MAKE_KRYLOV_SOLVER_HEADER(MethodName)                          \
  BOGUS_MAKE_KRYLOV_SOLVER_TYPEDEFS(MethodName)                              \
  MethodName(const Matrix &A, unsigned maxIters,                             \
             Scalar tol = NumTraits<Scalar>::epsilon(),                      \
             const Preconditioner *P = BOGUS_NULL_PTR(const Preconditioner), \
             const typename Base::SignalType *callback =                     \
                 BOGUS_NULL_PTR(const typename Base::SignalType))            \
      : Base(A, maxIters, tol, P, callback) {}                               \
                                                                             \
  MethodName() : Base() {}                                                   \
                                                                             \
  template <typename RhsT, typename ResT>                                    \
  Scalar vectorSolve(const RhsT &b, ResT x) const;

// Conjugate Gradient

//! Solves ( m_A * \p x = \p b ) using the Conjugate Gradient algorithm
/*! For symmetric matrices only. Converges for positive definite linear systems.

                <b>Matrix-vector mults/iter: </b> 1
                <b>Preconditionner calls/iter: </b> 1
                <b>Storage requirements: </b> 4n
        */
template <
    typename Matrix, typename Preconditioner = TrivialPreconditioner<Matrix>,
    typename Traits = ProblemTraits<typename MatrixTraits<Matrix>::Scalar> >
struct CG : public KrylovSolverBase<CG, Matrix, Preconditioner, Traits> {
  BOGUS_MAKE_KRYLOV_SOLVER_HEADER(CG)
};

//! Solves ( m_A * \p x = \p b ) using the BiConjugate Gradient algorithm
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
                Requires ability to perform transpose multiplication and
   preconditioning

                <b>Matrix-vector mults/iter: </b> 2 ( inc. 1 transpose )
                <b>Preconditionner calls/iter: </b> 2 ( inc. 1 transpose )
                <b>Storage requirements: </b> 8n
        */
template <
    typename Matrix, typename Preconditioner = TrivialPreconditioner<Matrix>,
    typename Traits = ProblemTraits<typename MatrixTraits<Matrix>::Scalar> >
struct BiCG : public KrylovSolverBase<BiCG, Matrix, Preconditioner, Traits> {
  BOGUS_MAKE_KRYLOV_SOLVER_HEADER(BiCG)
};

//! Solves ( m_A * \p x = \p b ) using the BiConjugate Gradient stabilized
//! algorithm
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
                Supposedly less erratic convergence than BiCG method.

                <b>Matrix-vector mults/iter: </b> 2
                <b>Preconditionner calls/iter: </b> 2
                <b>Storage requirements: </b> 8n

        \warning This implmentation seems to suffer from suspiciously bad
   convergence.

        */
template <
    typename Matrix, typename Preconditioner = TrivialPreconditioner<Matrix>,
    typename Traits = ProblemTraits<typename MatrixTraits<Matrix>::Scalar> >
struct BiCGSTAB
    : public KrylovSolverBase<BiCGSTAB, Matrix, Preconditioner, Traits> {
  BOGUS_MAKE_KRYLOV_SOLVER_HEADER(BiCGSTAB)
};

//! Solves ( m_A * \p x = \p b ) using the Conjugate Gradient Squared algorithm
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
                Supposedly faster convergence than BiCG \a when \a converging.

                <b>Matrix-vector mults/iter: </b> 2
                <b>Preconditionner calls/iter: </b> 2
                <b>Storage requirements: </b> 7n
        */
template <
    typename Matrix, typename Preconditioner = TrivialPreconditioner<Matrix>,
    typename Traits = ProblemTraits<typename MatrixTraits<Matrix>::Scalar> >
struct CGS : public KrylovSolverBase<CGS, Matrix, Preconditioner, Traits> {
  BOGUS_MAKE_KRYLOV_SOLVER_HEADER(CGS)
};

//! Solves ( m_A * \p x = \p b ) using the (restarted) Generalized Minimum
//! Residual
/*!
                \param restart If non-zero, use the GMRES(m) restarted
   algorithm. Lower the storage cost, but slows-down or even forbid the
   convergence of the algorithm.

                Works for non-symmetric linear systems.
                Probably the more robust method for non symmetric systems, but
   with the highest storage cost.

                <b>Matrix-vector mults/iter: </b> 1
                <b>Preconditionner calls/iter: </b> 1
                <b>Other ops/iter: </b> 1 k*k triangular solve, 2 k*n m-v mult,
   1 k*k m-v mult <b>Storage requirements: </b> 2*m*( n + m )

        */
template <
    typename Matrix, typename Preconditioner = TrivialPreconditioner<Matrix>,
    typename Traits = ProblemTraits<typename MatrixTraits<Matrix>::Scalar> >
struct GMRES : public KrylovSolverBase<GMRES, Matrix, Preconditioner, Traits> {
  BOGUS_MAKE_KRYLOV_SOLVER_TYPEDEFS(GMRES)

  GMRES() : Base(), m_restart(0) {}

  GMRES(const Matrix &A, unsigned maxIters,
        Scalar tol = NumTraits<Scalar>::epsilon(),
        const Preconditioner *P = BOGUS_NULL_PTR(const Preconditioner),
        const typename Base::SignalType *callback =
            BOGUS_NULL_PTR(const typename Base::SignalType),
        unsigned restart = 0)
      : Base(A, maxIters, tol, P, callback), m_restart(restart) {}

  GMRES &setRestart(unsigned restart) {
    m_restart = restart;
    return *this;
  }

  template <typename RhsT, typename ResT>
  Scalar vectorSolve(const RhsT &b, ResT x) const;

 protected:
  unsigned m_restart;
};

//! Solves ( m_A * \p x = \p b ) using the transpose-free Quasi Minimal Reisual
//! method
/*! Works for non-symmetric linear systems. Convergence not guaranteed.
                Supposedly less erratic convergence than BiCG method, but faster
   than BiCGSTAB.

                <b>Matrix-vector mults/iter: </b> 2
                <b>Preconditionner calls/iter: </b> 2
                <b>Storage requirements: </b> 7n

                \warning This function returns an approximation of the residual
   instead of the real one
        */
template <
    typename Matrix, typename Preconditioner = TrivialPreconditioner<Matrix>,
    typename Traits = ProblemTraits<typename MatrixTraits<Matrix>::Scalar> >
struct TFQMR : public KrylovSolverBase<TFQMR, Matrix, Preconditioner, Traits> {
  BOGUS_MAKE_KRYLOV_SOLVER_HEADER(TFQMR)
};

}  // namespace solvers

}  // namespace krylov

// Sovlers traits ( in bogus namespace )

#define BOGUS_PROCESS_KRYLOV_METHOD(MethodName)                                \
  template <typename Matrix, typename Preconditioner, class Traits>            \
  struct LinearSolverTraits<                                                   \
      krylov::solvers::MethodName<Matrix, Preconditioner, Traits> > {          \
    typedef Matrix MatrixType;                                                 \
    template <typename RhsT>                                                   \
    struct Result {                                                            \
      typedef typename Traits::template MutableClone<RhsT>::Type Type;         \
    };                                                                         \
  };                                                                           \
  template <typename Matrix, typename Preconditioner, class Traits,            \
            typename RhsBlockT, bool TransposeLhs, bool TransposeRhs>          \
  struct BlockBlockProductTraits<                                              \
      krylov::solvers::MethodName<Matrix, Preconditioner, Traits>, RhsBlockT,  \
      TransposeLhs, TransposeRhs> {                                            \
    typedef                                                                    \
        typename BlockBlockProductTraits<Matrix, RhsBlockT, TransposeLhs,      \
                                         TransposeRhs>::ReturnType ReturnType; \
  };

BOGUS_KRYLOV_METHODS
#undef BOGUS_PROCESS_KRYLOV_METHOD

}  // namespace bogus

#endif
