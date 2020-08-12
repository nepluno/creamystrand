/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "LinearSolver.hh"

#include <mkl.h>

#include "../Utils/TextLog.hh"
#include "BandMatrix.hh"

namespace strandsim {

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
int BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::s_kl = kl;
template <typename ScalarT, bool sym, bool pd, int kl, int ku>
int BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::s_ku = ku;
static const int s_one = 1;

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
BandMatrixCompressedStorage<ScalarT, sym, pd, kl,
                            ku>::BandMatrixCompressedStorage()
    : m_A(NULL), m_factorized(false) {}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
BandMatrixCompressedStorage<ScalarT, sym, pd, kl,
                            ku>::~BandMatrixCompressedStorage() {}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
void BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::store(
    const BandMatrixT &A, bool keepOldFactorization) {
  assert(A.rows() == A.cols());

  m_cols = A.cols();
  m_A = &A;

  m_rows = 0;

  if (!keepOldFactorization) {
    m_frows = 0;
    m_factorized = false;
  }
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
void BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::desallocate() {
  m_ab = MatrixT();
  m_factorization = MatrixT();
  m_factorized = false;
  m_rows = 0;
  m_cols = 0;
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
void BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::needsAb() {
  if (!m_rows) {
    assert(m_A);

    if (sym) {
      m_storageType = Upper;
      m_rows = ku + 1;
    } else {
      m_storageType = Full;
      m_rows = ku + kl + 1;
    }

    const std::vector<ScalarT> &data = m_A->getData();

    m_ab = Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic,
                         Eigen::RowMajor>::Map(&data[0], m_rows, m_cols);
  }
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
void BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::needsFab() {
  if (!m_frows) {
    assert(m_A);
    const std::vector<ScalarT> &data = m_A->getData();

    if (sym && pd) {
      m_fstorageType = Upper;
      m_frows = ku + 1;

      // Copy the data from A to m_factorization
      m_factorization =
          Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic,
                        Eigen::RowMajor>::Map(data.data(), m_frows, m_cols);

    } else {
      m_fstorageType = FullWithOffset;
      m_frows = 2 * kl + ku + 1;

      m_factorization.resize(m_frows, m_cols);
      m_ipiv.assign(m_cols, 0);

      // Copy the data from A to m_factorization
      m_factorization.block(0, 0, kl, m_cols).setZero();
      m_factorization.block(kl, 0, kl + ku + 1, m_cols) =
          Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic,
                        Eigen::RowMajor>::Map(data.data(), kl + ku + 1, m_cols);
    }
  }
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
int BandMatrixCompressedStorage<ScalarT, sym, pd, kl,
                                ku>::solveInPlaceFactorized(VectorT &x) const {
  assert(m_factorized);
  return BandMatrixMKLCaller<ScalarT, sym, pd, kl, ku>::solveFactorized(*this,
                                                                        x);
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
int BandMatrixCompressedStorage<ScalarT, sym, pd, kl,
                                ku>::solveInPlaceFactorized(MatrixT &x) const {
  assert(m_factorized);
  return BandMatrixMKLCaller<ScalarT, sym, pd, kl, ku>::solveFactorized(*this,
                                                                        x);
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
int BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::solveInPlace(
    MatrixT &x) {
  assert(x.rows() == cols());

  if (m_factorized) {
    return BandMatrixMKLCaller<ScalarT, sym, pd, kl, ku>::solveFactorized(*this,
                                                                          x);
  } else {
    needsFab();
    int info = BandMatrixMKLCaller<ScalarT, sym, pd, kl, ku>::solve(*this, x);
    m_factorized = !info;

    return info;
  }
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
int BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::solveInPlace(
    VectorT &x) {
  assert(x.rows() == cols());

  if (m_factorized) {
    return BandMatrixMKLCaller<ScalarT, sym, pd, kl, ku>::solveFactorized(*this,
                                                                          x);
  } else {
    needsFab();
    int info = BandMatrixMKLCaller<ScalarT, sym, pd, kl, ku>::solve(*this, x);
    m_factorized = !info;

    return info;
  }
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
void BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::multiply(
    VectorT &y, const VectorT &x, ScalarT alpha, ScalarT beta) {
  assert(y.rows() == cols());
  assert(x.rows() == cols());

  needsAb();
  BandMatrixMKLCaller<ScalarT, sym, pd, kl, ku>::multiply(*this, x, y, alpha,
                                                          beta);
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
void BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::multiply(
    MatrixT &y, const MatrixT &x, ScalarT alpha, ScalarT beta) {
  assert(y.cols() == x.cols());

  needsAb();

  VectorT tmpx(x.rows()), tmpy(y.rows());

  for (unsigned i = 0; i < y.cols(); ++i) {
    tmpy = y.col(i);
    tmpx = x.col(i);
    BandMatrixMKLCaller<ScalarT, sym, pd, kl, ku>::multiply(*this, tmpx, tmpy,
                                                            alpha, beta);
    y.col(i) = tmpy;
  }
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
int BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::factorize() {
  m_frows = 0;  // reload matrix data
  needsFab();

  int info = BandMatrixMKLCaller<ScalarT, sym, pd, kl, ku>::factorize(*this);
  m_factorized = !info;

  return info;
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
ScalarT
BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::invConditionNumber() {
  if (!m_factorized) factorize();

  if (m_factorized) {
    needsAb();

    // compute |A|_inf = max_i sum_j | Aij |
    const ScalarT aNorm = m_ab.cwiseAbs().colwise().sum().maxCoeff();

    return BandMatrixMKLCaller<ScalarT, sym, pd, kl, ku>::invConditionNumber(
        *this, aNorm);
  }

  return 0.;
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
int BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::computeEigenvalues(
    VectorT &lambda, MatrixT *vectors) {
  lambda.resize(cols());

  needsAb();
  return BandMatrixMKLCaller<ScalarT, sym, pd, kl, ku>::computeEigenvalues(
      *this, lambda, vectors);
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
int BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::conjugateGradient(
    Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> &b, unsigned maxSteps,
    ScalarT tol, PreconditionnerType pt) {
  // Allows to instantiate stuff
  switch (pt) {
    case Jacobi: {
      JacobiPreconditionner<
          BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku> >
          jp(*this);
      return conjugateGradient(x, b, maxSteps, tol, jp);
    }
    case Factorization: {
      FactorizationPreconditionner<
          BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku> >
          fp(*this);
      return conjugateGradient(x, b, maxSteps, tol, fp);
    }
    default: {
      DummyPreconditionner<
          BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku> >
          p;
      return conjugateGradient(x, b, maxSteps, tol, p);
    }
  };
}

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
template <typename PrecondT>
int BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku>::conjugateGradient(
    Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> &b, unsigned maxSteps,
    ScalarT tol, PrecondT &precond) {
  assert(b.rows() == cols());
  assert(x.rows() == cols());

  needsAb();
  return BandMatrixCG<ScalarT, sym, pd, kl, ku>::solve(*this, x, b, maxSteps,
                                                       tol * tol, precond);
}

// (Bi) COnjugate Gradient

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
struct BandMatrixCG {
  typedef BandMatrixCompressedStorage<ScalarT, sym, pd, kl, ku> Matrix;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> Vector;

  template <typename PrecondT>
  static int solve(Matrix &A, Vector &x, const Vector &b, int maxSteps,
                   ScalarT sqrTol, PrecondT &precond) {
    std::cerr << " Not implemented " << std::endl;
    return -1;
  }
};

// Symmetric positive definite

template <typename ScalarT, int kl, int ku>
struct BandMatrixCG<ScalarT, true, true, kl, ku> {
  typedef BandMatrixCompressedStorage<ScalarT, true, true, kl, ku> Matrix;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> Vector;

  template <typename PrecondT>
  static int solve(Matrix &A, Vector &x, const Vector &b, int maxSteps,
                   ScalarT sqrTol, PrecondT &precond) {
    return CG(A, x, b, maxSteps, sqrTol, precond);
  }

  template <typename PrecondT>
  static int CG(Matrix &A, Vector &x, const Vector &b, int maxSteps,
                ScalarT sqrTol, PrecondT &precond) {
    // Damped Jacobi precond
    // const Vector invDiag ( .5 * ( Vector( A.A().diagonal().array().inverse()
    // ) +  Vector::Ones( b.rows() ) ) ); const Vector invDiag ( ) );

    ScalarT alpha, beta;

    Vector r(b);
    A.multiply(r, x, -1., 1.);  // r = b - Ax

    Vector rD(b.rows()), pA(b.rows());

    precond(rD, r);
    Vector p(rD);

    ScalarT rDr, old_rDr = r.dot(rD);

    for (auto i = 0; i < maxSteps; ++i) {
      // std::cout << r.squaredNorm() << std::endl;

      A.multiply(pA, p);
      alpha = old_rDr / p.dot(pA);

      x += alpha * p;
      r -= alpha * pA;

      if (r.squaredNorm() < sqrTol) break;

      precond(rD, r);
      rDr = r.dot(rD);
      beta = rDr / old_rDr;
      old_rDr = rDr;

      p = rD + beta * p;
    }

    return 0;
  }
};

//  Symmetric undefinite

template <typename ScalarT, int kl, int ku>
struct BandMatrixCG<ScalarT, true, false, kl, ku> {
  typedef BandMatrixCompressedStorage<ScalarT, true, false, kl, ku> Matrix;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> Vector;

  template <typename PrecondT>
  static int solve(Matrix &A, Vector &x, const Vector &b, int maxSteps,
                   ScalarT sqrTol, PrecondT &precond) {
    return BiCGSTAB(A, x, b, maxSteps, sqrTol, precond);
  }

  template <typename PrecondT>
  static int CGNE(Matrix &A, Vector &x, const Vector &b, int maxSteps,
                  ScalarT sqrTol, PrecondT &precond) {
    // Damped Jacobi precond
    // const Vector invDiag ( .5 * ( Vector( A.A().diagonal().array().inverse()
    // ) +  Vector::Ones( b.rows() ) ) ); const Vector invDiag ( ) );

    ScalarT alpha, beta;

    Vector s(b);
    A.multiply(s, x, -1., 1.);  // s = b - Ax

    Vector r(b.rows());
    A.multiply(r, s);
    Vector p(r);
    Vector q(b.rows());
    A.multiply(q, p);

    ScalarT rDr, old_rDr = r.squaredNorm();

    for (unsigned i = 0; i < maxSteps; ++i) {
      // std::cout << i << " " << r.squaredNorm() << std::endl ;

      alpha = old_rDr / q.squaredNorm();

      s -= alpha * q;
      x += alpha * p;

      A.multiply(r, s);

      rDr = r.squaredNorm();
      beta = rDr / old_rDr;

      if (rDr < sqrTol) break;
      old_rDr = rDr;

      p = r + beta * p;
      A.multiply(q, p);
    }

    return 0;
  }

  template <typename PrecondT>
  static int BiCGSTAB(Matrix &A, Vector &x, const Vector &b, int maxSteps,
                      ScalarT sqrTol, PrecondT &precond) {
    ScalarT rho = 1., alpha = 1., omega = 1.;

    Vector r(b);
    A.multiply(r, x, -1., 1.);  // r = b - Ax

    const Vector r_hat(r);  // \hat{r} = r_0

    Vector nu(Vector::Zero(b.rows()));
    Vector p(nu);

    ScalarT rhom1, beta;
    Vector s(b.rows()), t(b.rows()), w, y, z;

    for (auto i = 0; i < maxSteps; ++i) {
      // std::cout << r.squaredNorm() << std::endl ;

      rhom1 = rho;
      rho = r_hat.dot(r);
      beta = (rho / rhom1) * (alpha / omega);

      p = r + beta * (p - omega * nu);
      precond(y, p);
      A.multiply(nu, y);

      alpha = rho / (r_hat.dot(nu));

      s = r - alpha * nu;
      precond(z, s);
      A.multiply(t, z);
      precond(w, t);

      omega = z.dot(w) / w.squaredNorm();

      x += alpha * y + omega * z;
      r = s - omega * t;

      if (r.squaredNorm() < sqrTol) break;
    }

    return 0;
  }
};

// Non-symmetric

template <typename ScalarT, int kl, int ku>
struct BandMatrixCG<ScalarT, false, false, kl, ku> {
  typedef BandMatrixCompressedStorage<ScalarT, false, false, kl, ku> Matrix;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> Vector;

  template <typename PrecondT>
  static int solve(Matrix &A, Vector &x, const Vector &b, int maxSteps,
                   ScalarT sqrTol, PrecondT &precond) {
    return BiCGSTAB(A, x, b, maxSteps, sqrTol);
  }

  static int BiCGSTAB(Matrix &A, Vector &x, const Vector &b, int maxSteps,
                      ScalarT sqrTol) {
    ScalarT rho = 1., alpha = 1., omega = 1.;

    Vector r(b);
    A.multiply(r, x, -1., 1.);  // r = b - Ax

    Vector r_hat(r);  // \hat{r} = r_0

    Vector nu(Vector::Zero(b.rows()));
    Vector p(nu);

    ScalarT rhom1, beta;
    Vector s(b.rows()), t(b.rows());

    for (unsigned i = 0; i < maxSteps; ++i) {
      rhom1 = rho;
      rho = r_hat.dot(r);
      beta = (rho / rhom1) * (alpha / omega);
      p = r + beta * (p - omega * nu);
      A.multiply(nu, p);
      alpha = rho / (r_hat.dot(nu));

      s = r - alpha * nu;
      A.multiply(t, s);
      omega = s.dot(t) / t.squaredNorm();
      x += alpha * p + omega * s;

      r = s - omega * t;

      if (r.squaredNorm() < sqrTol) break;
      // std::cout << r.squaredNorm() << std::endl ;
    }

    return 0;
  }
};

// Direct solve / multiply

template <typename ScalarT, bool sym, bool pd, int kl, int ku>
struct BandMatrixMKLCaller {};

////
// Specialized MKL BLAS/LAPACK calls

// General banded

template <int kl, int ku>
struct BandMatrixMKLCaller<double, false, false, kl, ku> {
  typedef double ScalarT;
  typedef BandMatrixCompressedStorage<ScalarT, false, false, kl, ku> Matrix;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> Mat;

  template <typename MatrixT>
  static int solve(Matrix &A, MatrixT &x) {
    assert(A.fstorage() == Matrix::FullWithOffset);

    /*       std::cerr << " Not implemented " << std::endl;
           return -1;*/

    // That probably  won't work, FIXME
    int info;
    const int nrhs = x.cols();
    /*      dgbsv_( &A.cols(), &Matrix::s_kl, &Matrix::s_ku, &nrhs, A.fab(),
       &A.frows(), A.ipiv(), x.data(), &A.cols(), &info );*/
    info =
        LAPACKE_dgbsv(LAPACK_COL_MAJOR, A.cols(), Matrix::s_kl, Matrix::s_ku,
                      nrhs, A.fab(), A.frows(), A.ipiv(), x.data(), A.cols());

    return info;
  }

  static int factorize(Matrix &A) {
    std::cerr << " Not implemented " << std::endl;
    return -1;
  }

  template <typename MatrixT>
  static int solveFactorized(const Matrix &A, MatrixT &x) {
    std::cerr << " Not implemented " << std::endl;
    return -1;
  }

  static Scalar invConditionNumber(Matrix &A, ScalarT aNorm) {
    std::vector<ScalarT> work(3 * A.cols());
    std::vector<int> ewok(A.cols());

    int info;
    ScalarT res;
    /*     dgbcon_( "I", &A.cols(), &Matrix::s_kl, &Matrix::s_ku, A.fab(),
       &A.frows(), A.ipiv(), &aNorm, &res, &work[0], &ewok[0], &info );*/
    info =
        LAPACKE_dgbcon(LAPACK_COL_MAJOR, 'I', A.cols(), Matrix::s_kl,
                       Matrix::s_ku, A.fab(), A.frows(), A.ipiv(), aNorm, &res);

    return res;
  }

  static int multiply(Matrix &A, const Vector &x, Vector &y, ScalarT alpha,
                      ScalarT beta) {
    assert(A.storage() == Matrix::Full);

    /*      dgbmv( "n", &A.cols(), &A.cols(), &Matrix::s_kl, &Matrix::s_ku,
       &alpha, A.ab(), &A.rows(), x.data(), &s_one, &beta, y.data(), &s_one );*/
    cblas_dgbmv(CblasColMajor, CblasNoTrans, A.cols(), A.cols(), Matrix::s_kl,
                Matrix::s_ku, alpha, A.fab(), A.frows(), x.data(), s_one, beta,
                y.data(), s_one);
    return -1;
  }

  static int computeEigenvalues(Matrix &A, Vector &lambda, Mat *vectors) {
    std::cerr << " Not implemented " << std::endl;
    return -1;
  }
};

template <int kl, int ku>
struct BandMatrixMKLCaller<float, false, false, kl, ku> {
  typedef float ScalarT;
  typedef BandMatrixCompressedStorage<ScalarT, false, false, kl, ku> Matrix;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> Mat;

  template <typename MatrixT>
  static int solve(Matrix &A, MatrixT &x) {
    assert(A.fstorage() == Matrix::FullWithOffset);

    /*       std::cerr << " Not implemented " << std::endl;
           return -1;*/

    // That probably  won't work, FIXME
    int info;
    const int nrhs = x.cols();
    /*     sgbsv_( &A.cols(), &Matrix::s_kl, &Matrix::s_ku, &nrhs, A.fab(),
       &A.frows(), A.ipiv(), x.data(), &A.cols(), &info );*/
    info =
        LAPACKE_sgbsv(LAPACK_COL_MAJOR, A.cols(), Matrix::s_kl, Matrix::s_ku,
                      nrhs, A.fab(), A.frows(), A.ipiv(), x.data(), A.cols());

    return info;
  }

  static int factorize(Matrix &A) {
    std::cerr << " Not implemented " << std::endl;
    return -1;
  }

  template <typename MatrixT>
  static int solveFactorized(const Matrix &A, MatrixT &x) {
    std::cerr << " Not implemented " << std::endl;
    return -1;
  }

  static Scalar invConditionNumber(Matrix &A, ScalarT aNorm) {
    std::vector<ScalarT> work(3 * A.cols());
    std::vector<int> ewok(A.cols());

    int info;
    ScalarT res;
    /*      sgbcon_( "I", &A.cols(), &Matrix::s_kl, &Matrix::s_ku, A.fab(),
       &A.frows(), A.ipiv(), &aNorm, &res, &work[0], &ewok[0], &info );*/
    info = LAPACKE_sgbcon(LAPACK_COL_MAJOR, 'I', A.cols(), Matrix::s_kl,
                          Matrix::s_ku, A.fab(), A.frows(), A.ipiv(), &aNorm,
                          &res);

    return res;
  }

  static int multiply(Matrix &A, const Vector &x, Vector &y, ScalarT alpha,
                      ScalarT beta) {
    assert(A.storage() == Matrix::Full);

    /*       sgbmv( "n", &A.cols(), &A.cols(), &Matrix::s_kl, &Matrix::s_ku,
       &alpha, A.ab(), &A.rows(), x.data(), &s_one, &beta, y.data(), &s_one );*/
    cblas_sgbmv(CblasColMajor, CblasNoTrans, A.cols(), A.cols(), Matrix::s_kl,
                Matrix::s_ku, alpha, A.ab(), A.rows(), x.data(), s_one, beta,
                y.data(), s_one);
    return 0;
  }

  static int computeEigenvalues(Matrix &A, Vector &lambda, Mat *vectors) {
    std::cerr << " Not implemented " << std::endl;
    return -1;
  }
};

// Symmetric non-positive definite

template <int kl, int ku>
struct BandMatrixMKLCaller<double, true, false, kl, ku> {
  typedef double ScalarT;
  typedef BandMatrixCompressedStorage<ScalarT, true, false, kl, ku> Matrix;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> Mat;

  template <typename MatrixT>
  static int solve(Matrix &A, MatrixT &x) {
    assert(A.fstorage() == Matrix::FullWithOffset);

    int info;
    const int nrhs = x.cols();
    /*        dgbsv_( &A.cols(), &Matrix::s_kl, &Matrix::s_ku, &nrhs, A.fab(),
       &A.frows(), A.ipiv(), x.data(), &A.cols(), &info );*/
    info =
        LAPACKE_dgbsv(LAPACK_COL_MAJOR, A.cols(), Matrix::s_kl, Matrix::s_ku,
                      nrhs, A.fab(), A.frows(), A.ipiv(), x.data(), A.cols());

    return info;
  }

  static int factorize(Matrix &A) {
    assert(A.fstorage() == Matrix::FullWithOffset);

    int info;
    /*     dgbtrf_( &A.cols(), &A.cols(), &Matrix::s_kl, &Matrix::s_ku, A.fab(),
       &A.frows(), A.ipiv(), &info );*/
    info = LAPACKE_dgbtrf(LAPACK_COL_MAJOR, A.cols(), A.cols(), Matrix::s_kl,
                          Matrix::s_ku, A.fab(), A.frows(), A.ipiv());

    return info;
  }

  template <typename MatrixT>
  static int solveFactorized(const Matrix &A, MatrixT &x) {
    assert(A.fstorage() == Matrix::FullWithOffset && A.factorized());

    int info;
    const int nrhs = x.cols();
    /*     dgbtrs_( "N", &A.cols(), &Matrix::s_kl, &Matrix::s_ku, &nrhs,
       A.fab(), &A.frows(), A.ipiv(), x.data(), &A.cols(), &info );*/
    info = LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', A.cols(), Matrix::s_kl,
                          Matrix::s_ku, nrhs, A.fab(), A.frows(), A.ipiv(),
                          x.data(), A.cols());

    return info;
  }

  static Scalar invConditionNumber(Matrix &A, ScalarT aNorm) {
    std::vector<ScalarT> work(3 * A.cols());
    std::vector<int> ewok(A.cols());

    int info;
    ScalarT res;
    /*      dgbcon_( "I", &A.cols(), &Matrix::s_kl, &Matrix::s_ku, A.fab(),
       &A.frows(), A.ipiv(), &aNorm, &res, &work[0], &ewok[0], &info );*/
    info =
        LAPACKE_dgbcon(LAPACK_COL_MAJOR, 'I', A.cols(), Matrix::s_kl,
                       Matrix::s_ku, A.fab(), A.frows(), A.ipiv(), aNorm, &res);

    return res;
  }

  static int multiply(Matrix &A, const Vector &x, Vector &y, ScalarT alpha,
                      ScalarT beta) {
    assert(A.storage() == Matrix::Upper);

    // dsbmv( "u", &A.cols(), &Matrix::s_ku, &alpha, A.ab(), &A.rows(),
    // x.data(), &s_one, &beta, y.data(), &s_one );
    cblas_dsbmv(CblasColMajor, CblasUpper, A.cols(), Matrix::s_ku, alpha,
                A.ab(), A.rows(), x.data(), s_one, beta, y.data(), s_one);

    return 0;
  }

  static int computeEigenvalues(Matrix &A, Vector &lambda, Mat *vectors) {
    assert(A.storage() == Matrix::Upper);

    char mode = vectors ? 'v' : 'n';
    ScalarT *z = vectors ? vectors->data() : NULL;
    int ldz = vectors ? A.cols() : 1;

    int info;
    /*
            int lwork( -1 );
            std::vector<ScalarT> work( 1 );
            int liwork( -1 );
            std::vector<int> ewok( 1 );

            // The first call to dsbevd is the workspace query, just computes
       the optimal lwork and liwork. dsbevd_( &mode, "u", &A.cols(),
       &Matrix::s_ku, A.ab(), &A.rows(), lambda.data(), z, &ldz, &work[0],
       &lwork, &ewok[0], &liwork, &info );*/
    info = LAPACKE_dsbevd(LAPACK_COL_MAJOR, mode, 'u', A.cols(), Matrix::s_ku,
                          A.ab(), A.rows(), lambda.data(), z, ldz);
    /*
            lwork = work[0];
            work.resize( lwork );
            liwork = ewok[0];
            ewok.resize( liwork );

            // The second call is the eigenvalue computation
            dsbevd_( &mode, "u", &A.cols(), &Matrix::s_ku, A.ab(), &A.rows(),
       lambda.data(), z, &ldz, &work[0], &lwork, &ewok[0], &liwork, &info );*/

    return info;
  }
};

template <int kl, int ku>
struct BandMatrixMKLCaller<float, true, false, kl, ku> {
  typedef float ScalarT;
  typedef BandMatrixCompressedStorage<ScalarT, true, false, kl, ku> Matrix;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> Mat;

  template <typename MatrixT>
  static int solve(Matrix &A, MatrixT &x) {
    assert(A.fstorage() == Matrix::FullWithOffset);

    int info;
    const int nrhs = x.cols();
    /*     sgbsv_( &A.cols(), &Matrix::s_kl, &Matrix::s_ku, &nrhs, A.fab(),
       &A.frows(), A.ipiv(), x.data(), &A.cols(), &info );*/
    info =
        LAPACKE_sgbsv(LAPACK_COL_MAJOR, A.cols(), Matrix::s_kl, Matrix::s_ku,
                      nrhs, A.fab(), A.frows(), A.ipiv(), x.data(), A.cols());

    return info;
  }

  static int factorize(Matrix &A) {
    assert(A.fstorage() == Matrix::FullWithOffset);

    int info;
    /*     sgbtrf_( &A.cols(), &A.cols(), &Matrix::s_kl, &Matrix::s_ku, A.fab(),
       &A.frows(), A.ipiv(), &info );*/
    info = LAPACKE_sgbtrf(LAPACK_COL_MAJOR, A.cols(), A.cols(), Matrix::s_kl,
                          Matrix::s_ku, A.fab(), A.frows(), A.ipiv());
    return info;
  }

  template <typename MatrixT>
  static int solveFactorized(const Matrix &A, MatrixT &x) {
    assert(A.fstorage() == Matrix::FullWithOffset && A.factorized());

    int info;
    const int nrhs = x.cols();
    /*       sgbtrs_( "N", &A.cols(), &Matrix::s_kl, &Matrix::s_ku, &nrhs,
       A.fab(), &A.frows(), A.ipiv(), x.data(), &A.cols(), &info );*/
    info = LAPACKE_sgbtrs(LAPACK_COL_MAJOR, 'N', A.cols(), Matrix::s_kl,
                          Matrix::s_ku, nrhs, A.fab(), A.frows(), A.ipiv(),
                          x.data(), A.cols());
    return info;
  }

  static Scalar invConditionNumber(Matrix &A, ScalarT aNorm) {
    std::vector<ScalarT> work(3 * A.cols());
    std::vector<int> ewok(A.cols());

    int info;
    ScalarT res;
    /*      sgbcon_( "I", &A.cols(), &Matrix::s_kl, &Matrix::s_ku, A.fab(),
       &A.frows(), A.ipiv(), &aNorm, &res, &work[0], &ewok[0], &info );*/
    info =
        LAPACKE_sgbcon(LAPACK_COL_MAJOR, 'I', A.cols(), Matrix::s_kl,
                       Matrix::s_ku, A.fab(), A.frows(), A.ipiv(), aNorm, &res);
    return res;
  }

  static int multiply(Matrix &A, const Vector &x, Vector &y, ScalarT alpha,
                      ScalarT beta) {
    assert(A.storage() == Matrix::Upper);

    /*       ssbmv( "u", &A.cols(), &Matrix::s_ku, &alpha, A.ab(), &A.rows(),
       x.data(), &s_one, &beta, y.data(), &s_one );*/
    cblas_ssbmv(CblasColMajor, CblasUpper, A.cols(), Matrix::s_ku, alpha,
                A.ab(), A.rows(), x.data(), s_one, beta, y.data(), s_one);
    return 0;
  }

  static int computeEigenvalues(Matrix &A, Vector &lambda, Mat *vectors) {
    assert(A.storage() == Matrix::Upper);

    char mode = vectors ? 'v' : 'n';
    ScalarT *z = vectors ? vectors->data() : NULL;
    int ldz = vectors ? A.cols() : 1;

    int info;
    /*
            int lwork( -1 );
            std::vector<ScalarT> work( 1 );
            int liwork( -1 );
            std::vector<int> ewok( 1 );

            ssbevd_( &mode, "u", &A.cols(), &Matrix::s_ku, A.ab(), &A.rows(),
       lambda.data(), z, &ldz, &work[0], &lwork, &ewok[0], &liwork, &info );*/
    info = LAPACKE_ssbevd(LAPACK_COL_MAJOR, mode, 'u', A.cols(), Matrix::s_ku,
                          A.ab(), A.rows(), lambda.data(), z, ldz);
    /*
            lwork = work[0];
            work.resize( lwork );
            liwork = ewok[0];
            ewok.resize( liwork );

            ssbevd_( &mode, "u", &A.cols(), &Matrix::s_ku, A.ab(), &A.rows(),
       lambda.data(), z, &ldz, &work[0], &lwork, &ewok[0], &liwork, &info );*/

    return info;
  }
};

// Symmetric Positive definite

template <int kl, int ku>
struct BandMatrixMKLCaller<double, true, true, kl, ku> {
  typedef double ScalarT;
  typedef BandMatrixCompressedStorage<ScalarT, true, true, kl, ku> Matrix;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> Mat;

  template <typename MatrixT>
  static int solve(Matrix &A, MatrixT &x) {
    assert(A.fstorage() == Matrix::Upper);

    int info;
    const int nrhs = x.cols();
    /*       dpbsv_( "u", &A.cols(), &Matrix::s_ku, &nrhs, A.fab(), &A.frows(),
       x.data(), &A.cols(), &info );*/
    info = LAPACKE_dpbsv(LAPACK_COL_MAJOR, 'u', A.cols(), Matrix::s_ku, nrhs,
                         A.fab(), A.frows(), x.data(), A.cols());

    return info;
  }

  static int factorize(Matrix &A) {
    assert(A.fstorage() == Matrix::Upper);

    int info;
    /*      dpbtrf_( "u", &A.cols(), &Matrix::s_ku, A.fab(), &A.frows(), &info
     * );*/
    info = LAPACKE_dpbtrf(LAPACK_COL_MAJOR, 'u', A.cols(), Matrix::s_ku,
                          A.fab(), A.frows());

    return info;
  }

  template <typename MatrixT>
  static int solveFactorized(const Matrix &A, MatrixT &x) {
    assert(A.fstorage() == Matrix::Upper && A.factorized());

    int info;
    const int nrhs = x.cols();
    /*       dpbtrs_( "u", &A.cols(), &Matrix::s_ku, &nrhs, A.fab(), &A.frows(),
       x.data(), &A.cols(), &info );*/
    info = LAPACKE_dpbtrs(LAPACK_COL_MAJOR, 'u', A.cols(), Matrix::s_ku, nrhs,
                          A.fab(), A.frows(), x.data(), A.cols());
    return info;
  }

  static Scalar invConditionNumber(Matrix &A, ScalarT aNorm) {
    std::vector<ScalarT> work(3 * A.cols());
    std::vector<int> ewok(A.cols());

    int info;
    ScalarT res;
    /*      dpbcon_( "U", &A.cols(), &Matrix::s_ku, A.fab(), &A.frows(), &aNorm,
       &res, &work[0], &ewok[0], &info );*/
    LAPACKE_dpbcon(LAPACK_COL_MAJOR, 'u', A.cols(), Matrix::s_ku, A.fab(),
                   A.frows(), aNorm, &res);
    return res;
  }

  static int multiply(Matrix &A, const Vector &x, Vector &y, ScalarT alpha,
                      ScalarT beta) {
    assert(A.storage() == Matrix::Upper);

    // dsbmv( "u", &A.cols(), &Matrix::s_ku, &alpha, A.ab(), &A.rows(),
    // x.data(), &s_one, &beta,
    // y.data(), &s_one );
    cblas_dsbmv(CblasColMajor, CblasUpper, A.cols(), Matrix::s_ku, alpha,
                A.ab(), A.rows(), x.data(), s_one, beta, y.data(), s_one);
    return 0;
  }

  static int computeEigenvalues(Matrix &A, Vector &lambda, Mat *vectors) {
    assert(A.storage() == Matrix::Upper);

    char mode = vectors ? 'v' : 'n';
    ScalarT *z = vectors ? vectors->data() : NULL;
    int ldz = vectors ? A.cols() : 1;

    int info;

    /*       int lwork( -1 );
           std::vector<ScalarT> work( 1 );
           int liwork( -1 );
           std::vector<int> ewok( 1 );

           dsbevd_( &mode, "u", &A.cols(), &Matrix::s_ku, A.ab(), &A.rows(),
       lambda.data(), z, &ldz, &work[0], &lwork, &ewok[0], &liwork, &info );*/

    info = LAPACKE_dsbevd(LAPACK_COL_MAJOR, mode, 'u', A.cols(), Matrix::s_ku,
                          A.ab(), A.rows(), lambda.data(), z, ldz);

    /*       lwork = work[0];
           work.resize( lwork );
           liwork = ewok[0];
           ewok.resize( liwork );

           dsbevd_( &mode, "u", &A.cols(), &Matrix::s_ku, A.ab(), &A.rows(),
       lambda.data(), z, &ldz, &work[0], &lwork, &ewok[0], &liwork, &info );*/

    return info;
  }
};

template <int kl, int ku>
struct BandMatrixMKLCaller<float, true, true, kl, ku> {
  typedef float ScalarT;
  typedef BandMatrixCompressedStorage<ScalarT, true, true, kl, ku> Matrix;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> Mat;

  template <typename MatrixT>
  static int solve(Matrix &A, MatrixT &x) {
    assert(A.fstorage() == Matrix::Upper);

    int info;
    const int nrhs = x.cols();
    /*     spbsv_( "u", &A.cols(), &Matrix::s_ku, &nrhs, A.fab(), &A.frows(),
       x.data(), &A.cols(), &info );*/
    info = LAPACKE_spbsv(LAPACK_COL_MAJOR, 'u', A.cols(), Matrix::s_ku, nrhs,
                         A.fab(), A.frows(), x.data(), A.cols());

    return info;
  }

  static int factorize(Matrix &A) {
    assert(A.fstorage() == Matrix::Upper);

    int info;
    /*       spbtrf_( "u", &A.cols(), &Matrix::s_ku, A.fab(), &A.frows(), &info
     * );*/
    info = LAPACKE_spbtrf(LAPACK_COL_MAJOR, 'u', A.cols(), Matrix::s_ku,
                          A.fab(), A.frows());

    return info;
  }

  template <typename MatrixT>
  static int solveFactorized(const Matrix &A, MatrixT &x) {
    assert(A.fstorage() == Matrix::Upper && A.factorized());

    int info;
    const int nrhs = x.cols();
    /*        spbtrs_( "u", &A.cols(), &Matrix::s_ku, &nrhs, A.fab(),
       &A.frows(), x.data(), &A.cols(), &info );*/
    info = LAPACKE_spbtrs(LAPACK_COL_MAJOR, 'u', A.cols(), Matrix::s_ku, nrhs,
                          A.fab(), A.frows(), x.data(), A.cols());

    return info;
  }

  static Scalar invConditionNumber(Matrix &A, ScalarT aNorm) {
    std::vector<ScalarT> work(3 * A.cols());
    std::vector<int> ewok(A.cols());

    int info;
    ScalarT res;
    /*       spbcon_( "U", &A.cols(), &Matrix::s_ku, A.fab(), &A.frows(),
       &aNorm, &res, &work[0], &ewok[0], &info );*/
    LAPACKE_spbcon(LAPACK_COL_MAJOR, 'u', A.cols(), Matrix::s_ku, A.fab(),
                   A.frows(), aNorm, &res);

    return res;
  }

  static int multiply(Matrix &A, const Vector &x, Vector &y, ScalarT alpha,
                      ScalarT beta) {
    assert(A.storage() == Matrix::Upper);

    /*      ssbmv( "u", &A.cols(), &Matrix::s_ku, &alpha, A.ab(), &A.rows(),
       x.data(), &s_one, &beta, y.data(), &s_one );*/
    cblas_ssbmv(CblasColMajor, CblasUpper, A.cols(), Matrix::s_ku, alpha,
                A.ab(), A.rows(), x.data(), s_one, beta, y.data(), s_one);

    return 0;
  }

  static int computeEigenvalues(Matrix &A, Vector &lambda, Mat *vectors) {
    assert(A.storage() == Matrix::Upper);

    char mode = vectors ? 'v' : 'n';
    ScalarT *z = vectors ? vectors->data() : NULL;
    int ldz = vectors ? A.cols() : 1;

    int info;
    /*
int lwork( -1 );
std::vector<ScalarT> work( 1 );
int liwork( -1 );
std::vector<int> ewok( 1 );

ssbevd_( &mode, "u", &A.cols(), &Matrix::s_ku, A.ab(), &A.rows(), lambda.data(),
z, &ldz, &work[0], &lwork, &ewok[0], &liwork, &info );*/

    LAPACKE_ssbevd(LAPACK_COL_MAJOR, mode, 'u', A.cols(), Matrix::s_ku, A.ab(),
                   A.rows(), lambda.data(), z, ldz);

    /*       lwork = work[0];
           work.resize( lwork );
           liwork = ewok[0];
           ewok.resize( liwork );

           ssbevd_( &mode, "u", &A.cols(), &Matrix::s_ku, A.ab(), &A.rows(),
       lambda.data(), z, &ldz, &work[0], &lwork, &ewok[0], &liwork, &info );*/

    return info;
  }
};

// template class BandMatrixCompressedStorage< Scalar, false, false, 10, 10 > ;
template class BandMatrixCompressedStorage<Scalar, true, false, 10, 10>;
template class BandMatrixCompressedStorage<Scalar, true, true, 10, 10>;
template class BandMatrixCompressedStorage<Scalar, true, true, 15, 15>;
template class BandMatrixCompressedStorage<Scalar, true, false, 1, 1>;

}  // namespace strandsim
