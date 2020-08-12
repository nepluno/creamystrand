/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_EIGEN_SPARSE_LINEAR_SOLVERS_HPP
#define BOGUS_EIGEN_SPARSE_LINEAR_SOLVERS_HPP

/*  Depending on your version of Eigen, this file may make use of LGPL licensed
   code. See Eigen/src/Sparse/SimplicialCholesky.h for more information
*/

#include "SparseHeader.hpp"

#define BOGUS_WITH_EIGEN_STABLE_SPARSE_API
#ifdef BOGUS_WITH_EIGEN_STABLE_SPARSE_API

#include <Eigen/OrderingMethods>

#include "../Utils/LinearSolverBase.hpp"
#include "../Utils/NaiveSharedPtr.hpp"
#include "EigenLinearSolvers.hpp"

#ifndef EIGEN_MPL2_ONLY
#include <Eigen/SparseCholesky>
#define BOGUS_WITH_EIGEN_SPARSE_LDLT
#endif

#if EIGEN_VERSION_AT_LEAST(3, 1, 92)
#include <Eigen/SparseLU>
#define BOGUS_WITH_EIGEN_SPARSE_LU
#endif

#ifndef BOGUS_EIGEN_SPARSE_LDLT_IMPL
#define BOGUS_EIGEN_SPARSE_LDLT_IMPL(Derived) Eigen::SimplicialLDLT<Derived>
#endif

#ifndef BOGUS_EIGEN_SPARSE_LU_IMPL
#define BOGUS_EIGEN_SPARSE_LU_IMPL(Derived) \
  Eigen::SparseLU<Derived, Eigen::COLAMDOrdering<int> >
#endif

namespace bogus {

template <typename Derived, typename ImplType>
struct LinearSolverTraits<
    Factorization<Eigen::SparseMatrixBase<Derived>, ImplType> > {
  typedef typename Derived::PlainObject MatrixType;
  typedef ImplType FactType;

  template <typename RhsT>
  struct Result {
    typedef typename EigenSolveResult<FactType, RhsT>::Type Type;
  };
  template <typename RhsT>
  struct Result<Eigen::MatrixBase<RhsT> > {
    typedef typename Result<RhsT>::Type Type;
  };
};

template <typename Derived, typename FactType>
struct Factorization<Eigen::SparseMatrixBase<Derived>, FactType>
    : public LinearSolverBase<
          Factorization<Eigen::SparseMatrixBase<Derived>, FactType> > {
  typedef Eigen::SparseMatrixBase<Derived> MatrixType;
  typedef LinearSolverTraits<Factorization<MatrixType, FactType> > Traits;

  Factorization() {}
  template <typename OtherDerived>
  explicit Factorization(const Eigen::SparseMatrixBase<OtherDerived>& mat)
      : m_fact(new typename Traits::FactType(mat)) {}

  template <typename OtherDerived>
  Factorization& compute(const Eigen::SparseMatrixBase<OtherDerived>& mat) {
    m_fact.reset(new typename Traits::FactType(mat));
    return *this;
  }

  template <typename RhsT, typename ResT>
  void solve(const Eigen::MatrixBase<RhsT>& rhs, ResT& res) const {
    assert(m_fact);
    res = m_fact->solve(rhs);
  }

  template <typename RhsT>
  typename Traits::template Result<Eigen::MatrixBase<RhsT> >::Type solve(
      const Eigen::MatrixBase<RhsT>& rhs) const {
    assert(m_fact);
    return m_fact->solve(rhs);
  }

  const typename Traits::FactType& factorization() const { return *m_fact; }

 protected:
  BOGUS_SHARED_PTR(typename Traits::FactType, m_fact);
};

template <typename Derived, typename FactType>
struct EigenSparseFactorization
    : public Factorization<Eigen::SparseMatrixBase<Derived>, FactType> {
  EigenSparseFactorization() {}
  explicit EigenSparseFactorization(const Eigen::SparseMatrixBase<Derived>& mat)
      : Factorization<Eigen::SparseMatrixBase<Derived>, FactType>(mat) {}
};

// LDLT

#ifdef BOGUS_WITH_EIGEN_SPARSE_LDLT
#if !EIGEN_VERSION_AT_LEAST(3, 2, 90)
template <typename MatrixType, typename RhsType>
struct EigenSolveResult<Eigen::SimplicialLDLT<MatrixType>, RhsType> {
  typedef Eigen::SimplicialLDLT<MatrixType> FactType;
  typedef Eigen::internal::solve_retval<Eigen::SimplicialCholeskyBase<FactType>,
                                        RhsType>
      Type;
};
#endif
#endif

template <typename Derived>
struct LinearSolverTraits<LDLT<Eigen::SparseMatrixBase<Derived> > >
    : public LinearSolverTraits<
          Factorization<Eigen::SparseMatrixBase<Derived>,
                        BOGUS_EIGEN_SPARSE_LDLT_IMPL(Derived)> > {};
template <typename Derived>
struct LinearSolverTraits<LU<Eigen::SparseMatrixBase<Derived> > >
    : public LinearSolverTraits<
          Factorization<Eigen::SparseMatrixBase<Derived>,
                        BOGUS_EIGEN_SPARSE_LU_IMPL(Derived)> > {};

template <typename Derived>
struct LDLT<Eigen::SparseMatrixBase<Derived> >
    : public Factorization<Eigen::SparseMatrixBase<Derived>,
                           typename LinearSolverTraits<LDLT<
                               Eigen::SparseMatrixBase<Derived> > >::FactType> {
  typedef Eigen::SparseMatrixBase<Derived> MatrixType;
  typedef LinearSolverTraits<LDLT<MatrixType> > Traits;

  typedef Factorization<MatrixType, typename Traits::FactType> Base;

  LDLT() : Base() {}
  template <typename OtherDerived>
  explicit LDLT(const Eigen::SparseMatrixBase<OtherDerived>& mat) : Base(mat) {}
};

template <typename Scalar, int _Options = 0, typename _Index = int>
struct SparseLDLT : public LDLT<Eigen::SparseMatrixBase<
                        Eigen::SparseMatrix<Scalar, _Options, _Index> > > {
  SparseLDLT() {}
  template <typename OtherDerived>
  explicit SparseLDLT(const Eigen::SparseMatrixBase<OtherDerived>& mat)
      : LDLT<Eigen::SparseMatrixBase<
            Eigen::SparseMatrix<Scalar, _Options, _Index> > >(mat) {}
};

template <typename Derived>
struct LU<Eigen::SparseMatrixBase<Derived> >
    : public Factorization<Eigen::SparseMatrixBase<Derived>,
                           typename LinearSolverTraits<LU<
                               Eigen::SparseMatrixBase<Derived> > >::FactType> {
  typedef Eigen::SparseMatrixBase<Derived> MatrixType;
  typedef LinearSolverTraits<LU<MatrixType> > Traits;

  typedef Factorization<MatrixType, typename Traits::FactType> Base;

  LU() : Base() {}
  template <typename OtherDerived>
  explicit LU(const Eigen::SparseMatrixBase<OtherDerived>& mat) : Base(mat) {}
};

template <typename Scalar, int _Options = 0, typename _Index = int>
struct SparseLU : public LU<Eigen::SparseMatrixBase<
                      Eigen::SparseMatrix<Scalar, _Options, _Index> > > {
  SparseLU() {}
  template <typename OtherDerived>
  explicit SparseLU(const Eigen::SparseMatrixBase<OtherDerived>& mat)
      : LU<Eigen::SparseMatrixBase<
            Eigen::SparseMatrix<Scalar, _Options, _Index> > >(mat) {}
};

}  // namespace bogus

#endif  // EIGEN_STABLE_API

#endif  // HPP
