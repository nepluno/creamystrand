/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef LINEARSOLVER_HH_
#define LINEARSOLVER_HH_

#include "BandMatrix.hh"
#include "Definitions.hh"

namespace strandsim {

template <typename BandMatrixCompressedStorageT>
class DummyPreconditionner;
template <typename BandMatrixCompressedStorageT>
class JacobiPreconditionner;
template <typename BandMatrixCompressedStorageT>
class FactorizationPreconditionner;

// This template class is defined in two partial template specialisations to
// cover different MKL-LAPACK routine calls

template <typename ScalarT, bool symmetric, bool positiveDefinite, int kl,
          int ku>
class BandMatrixCompressedStorage {
 public:
  typedef BandMatrix<ScalarT, kl, ku> BandMatrixT;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> VectorT;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> MatrixT;

  enum PreconditionnerType { NoPreconditionner, Jacobi, Factorization };

  BandMatrixCompressedStorage();
  BandMatrixCompressedStorage(const BandMatrixT& A) { store(A); }

  ~BandMatrixCompressedStorage();

  void store(const BandMatrixT& A, bool keepOldFactorization = false);
  void desallocate();

  // x = A \ b . store( A ) should have been called.
  template <typename Derived>
  int solve(MatrixT& x, const Eigen::MatrixBase<Derived>& b) {
    return solveInPlace(x = b);
  }
  template <typename Derived>
  int solve(VectorT& x, const Eigen::MatrixBase<Derived>& b) {
    return solveInPlace(x = b);
  }

  // x = A \ x . store( A ) should have been called.
  int solveInPlace(MatrixT& x);
  int solveInPlace(VectorT& x);

  // x = A \ x . store( A ) and factorize() should have been called.
  int solveInPlaceFactorized(MatrixT& x) const;
  int solveInPlaceFactorized(VectorT& x) const;

  // y = alpha * A * x + beta * y.  store( A, false) should have been called.
  void multiply(VectorT& y, const VectorT& x, ScalarT alpha = 1.,
                ScalarT beta = 0.);
  void multiply(MatrixT& y, const MatrixT& x, ScalarT alpha = 1.,
                ScalarT beta = 0.);

  int factorize();

  ScalarT invConditionNumber();

  int computeEigenvalues(VectorT& lambda, MatrixT* vectors = NULL);

  // x = A \ b using input x as guess. store( A, false) should have been called.
  int conjugateGradient(VectorT& x, const VectorT& b, unsigned maxSteps,
                        ScalarT tol, PreconditionnerType = NoPreconditionner);

  template <typename PrecondT>
  int conjugateGradient(VectorT& x, const VectorT& b, unsigned maxSteps,
                        ScalarT tol, PrecondT& precond);

  const BandMatrixT& A() const {
    assert(m_A);
    return *m_A;
  }

  const int& cols() const { return m_cols; }
  const int& rows() const { return m_rows; }
  const int& frows() const { return m_frows; }
  int& cols() { return m_cols; }
  int& rows() { return m_rows; }
  int& frows() { return m_frows; }

  bool AisSet() const { return m_A != NULL; }

  bool factorized() const { return m_factorized; }

  const ScalarT* ab() const { return m_ab.data(); }
  ScalarT* ab() { return m_ab.data(); }
  const ScalarT* fab() const { return m_factorization.data(); }
  ScalarT* fab() { return m_factorization.data(); }

  const int* ipiv() const { return &m_ipiv[0]; }
  int* ipiv() { return &m_ipiv[0]; }

  static int s_kl;  // Just because dgbsv wants a pointer
  static int s_ku;

  enum StorageType { Upper, Full, FullWithOffset };
  StorageType storage() const { return m_storageType; }
  StorageType fstorage() const { return m_fstorageType; }

  Scalar norm() const { return m_ab.norm(); }
  Scalar fnorm() const { return m_factorization.norm(); }

  const MatrixT& bands() const { return m_ab; }
  const MatrixT& factorizationBands() const { return m_factorization; }

  void needsAb();
  void needsFab();

 private:
  const BandMatrixT* m_A;

  MatrixT m_ab;
  MatrixT m_factorization;

  std::vector<int> m_ipiv;

  int m_cols;
  int m_rows;
  int m_frows;
  StorageType m_storageType;
  StorageType m_fstorageType;

  bool m_factorized;
};

// Preconditionners

template <typename BMCST>
class DummyPreconditionner {
 public:
  void operator()(typename BMCST::VectorT& lhs,
                  const typename BMCST::VectorT& rhs) {
    lhs = rhs;
  }
};

template <typename BMCST>
class JacobiPreconditionner {
 public:
  JacobiPreconditionner(BMCST& A) : m_invDiag(A.A().diagonal().array()) {}

  void operator()(typename BMCST::VectorT& lhs,
                  const typename BMCST::VectorT& rhs) {
    lhs = m_invDiag.cwiseProduct(rhs);
  }

 private:
  typename BMCST::VectorT m_invDiag;
};

template <typename BMCST>
class FactorizationPreconditionner {
 public:
  FactorizationPreconditionner(BMCST& A) : m_A(A) {}

  void operator()(typename BMCST::VectorT& lhs,
                  const typename BMCST::VectorT& rhs) {
    m_A.solve(lhs, rhs);
  }

 private:
  BMCST& m_A;
};

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
struct BandMatrixMKLCaller;
template <typename ScalarT, bool sym, bool pd, int kl, int ku>
struct BandMatrixCG;

typedef BandMatrixCompressedStorage<Scalar, true, false,
                                    JacobianMatrixType::LowerBands,
                                    JacobianMatrixType::UpperBands>
    JacobianStorage;
typedef BandMatrixCompressedStorage<Scalar, true, true,
                                    JacobianMatrixType::LowerBands,
                                    JacobianMatrixType::UpperBands>
    SPDJacobianStorage;
typedef BandMatrixCompressedStorage<Scalar, true, true, 15, 15>
    RestJacobianStorage;
typedef BandMatrixCompressedStorage<Scalar, true, false, 1, 1>
    TriDiagonalStorage;

}  // namespace strandsim

#endif /* LINEARSOLVER_HH_ */
