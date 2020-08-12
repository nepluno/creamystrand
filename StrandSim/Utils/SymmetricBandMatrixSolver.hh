/**
 * \copyright 2014 Danny Kaufman
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SYMMETRICBANDMATRIXSOLVER_HH
#define SYMMETRICBANDMATRIXSOLVER_HH

#include <boost/serialization/split_member.hpp>

#include "../Core/BandMatrix.hh"
#include "../Core/Definitions.hh"
#include "../Core/LinearSolver.hh"

namespace strandsim {

template <typename ScalarT, int kl>
class SymmetricBandMatrixSolver {
 public:
  typedef BandMatrixCompressedStorage<ScalarT, true, true, kl, kl> SPDStorageT;
  typedef BandMatrixCompressedStorage<ScalarT, true, false, kl, kl>
      NonSPDStorageT;
  typedef BandMatrix<ScalarT, kl, kl> BandMatrixT;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> VectorT;
  typedef Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> MatrixT;

  SymmetricBandMatrixSolver()
      : m_notSPD(false), m_invScaling(1.), m_allocatedBandMatrix(NULL) {}

  SymmetricBandMatrixSolver(const BandMatrixT& A)
      : m_notSPD(false), m_invScaling(1.), m_allocatedBandMatrix(NULL) {
    store(A);
  }

  ~SymmetricBandMatrixSolver() { delete m_allocatedBandMatrix; }

  void store(const BandMatrixT& A) {
    setScaling(1.);

    if (m_notSPD) {
      m_nonSPDStorage.desallocate();
    }

    m_SPDStorage.store(A);
    m_notSPD = m_SPDStorage.factorize();
    if (m_notSPD) {
      m_SPDStorage.desallocate();
      m_nonSPDStorage.store(A);
      m_nonSPDStorage.factorize();
    }
  }

  void setScaling(Scalar scaling) { m_invScaling = 1. / scaling; }

  bool notSPD() const { return m_notSPD; }

  // x = A \ b . store( A ) should have been called.
  template <typename Derived, typename OtherDerived>
  int solve(Eigen::MatrixBase<Derived>& x,
            const Eigen::MatrixBase<OtherDerived>& b) const {
    x = m_invScaling * b;

    if (m_notSPD) return m_nonSPDStorage.solveInPlaceFactorized(x.derived());
    return m_SPDStorage.solveInPlaceFactorized(x.derived());
  }

  SPDStorageT& SPDStorage() { return m_SPDStorage; }
  NonSPDStorageT& nonSPDStorage() { return m_nonSPDStorage; }

  const BandMatrixT& matrix() const { return m_SPDStorage.A(); }

  template <class Archive>
  void load(Archive& ar, const unsigned int version) {
    delete m_allocatedBandMatrix;
    m_allocatedBandMatrix = new BandMatrixT();

    // invScaling will be overwritten in store(), so we retrieve it in a
    // temporary
    Scalar invScaling;
    ar >> invScaling;
    ar >> *m_allocatedBandMatrix;

    store(*m_allocatedBandMatrix);
    m_invScaling = invScaling;
  }
  template <class Archive>
  void save(Archive& ar, const unsigned int version) const {
    ar << m_invScaling;
    ar << matrix();
  }
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    boost::serialization::split_member(ar, *this, version);
  }

 private:
  SymmetricBandMatrixSolver(const SymmetricBandMatrixSolver&);
  SymmetricBandMatrixSolver& operator=(const SymmetricBandMatrixSolver&);

  SPDStorageT m_SPDStorage;
  NonSPDStorageT m_nonSPDStorage;

  bool m_notSPD;
  Scalar m_invScaling;

  BandMatrixT* m_allocatedBandMatrix;  // If loaded from archive
};

typedef SymmetricBandMatrixSolver<Scalar, JacobianMatrixType::LowerBands>
    JacobianSolver;

}  // namespace strandsim

#endif  // SYMMETRICBANDMATRIXSOLVER_HH
