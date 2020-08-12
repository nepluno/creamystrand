/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BANDMATRIX_HH_
#define BANDMATRIX_HH_

#include "BandMatrixFwd.hh"

namespace strandsim {

template <typename ScalarT, IndexType kl, IndexType ku>
class BandMatrix {
 public:
  static const int LowerBands = kl;
  static const int UpperBands = ku;

  explicit BandMatrix(const IndexType rows = 0, const IndexType cols = 0)
      : m_rows(rows), m_cols(cols), m_size((kl + ku + 1) * cols) {
    m_data.resize(m_size);
  }

  ~BandMatrix() {}

  void resize(const IndexType rows, const IndexType cols) {
    m_rows = rows;
    m_cols = cols;
    m_size = (kl + ku + 1) * cols;
    m_data.resize(m_size);
  }

  const ScalarT operator()(const IndexType i, const IndexType j) const {
    if (!indicesValid(i, j)) return 0.0;

    return m_data[(ku + i - j) * m_cols + j];
  }

  ScalarT& operator()(IndexType i, IndexType j) {
    assert(indicesValid(i, j));

    return m_data[(ku + i - j) * m_cols + j];
  }

  VecXx::ConstMapType diagonal() const {
    assert(m_cols == m_rows);
    return VecXx::Map(&m_data[ku * m_cols], m_cols);
  }

  VecXx::MapType diagonal() {
    assert(m_cols == m_rows);
    return VecXx::Map(&m_data[ku * m_cols], m_cols);
  }

  // Add lambda * Identity
  void addConstantDiagonal(const Scalar lambda) {
    diagonal().array() += lambda;
  }

  // Returns column i
  VecXx col(int i) const {
    VecXx col(m_rows);
    col.setZero();

    const unsigned start = std::max(0, (int)i - (int)ku);
    const unsigned len =
        std::min((unsigned)(kl + ku + 1), (unsigned)(m_rows - 1 - start));

    Eigen::Stride<1, Eigen::Dynamic> stride(1, m_cols);
    col.segment(start, len) = VecXx::Map(&m_data[i], len, stride);

    return col;
  }

  template <int nfixed>
  void fixFirstDOFs(int startDof = 0) {
    assert(startDof + nfixed <= m_rows && startDof + nfixed <= m_cols);

    for (int i = startDof; i < startDof + nfixed; i++) {
      for (int j = std::max(0, (int)i - (int)kl);
           j <= std::min((int)m_cols - 1, (int)i + (int)ku); ++j)
        (*this)(i, j) = 0.0;
      for (int k = std::max(0, (int)i - (int)ku);
           k <= std::min((int)m_rows - 1, (int)i + (int)kl); k++)
        (*this)(k, i) = 0.0;
      (*this)(i, i) = 1.0;
    }
  }

  // Add the local Jacobian to the diagonal part of matrix, top left corner at
  // "start" position on the diagonal.
  void localDiagonalAdd(int start, const ScalarT& localJ) {
    m_data[start + m_cols * ku] += localJ;
  }
  // Add the local Jacobian to the banded matrix, top left corner at "start"
  // position on the diagonal.
  template <IndexType localSize>
  void localStencilAdd(
      int start, const Eigen::Matrix<ScalarT, localSize, localSize>& localJ) {
    start += m_cols * ku;

    for (int i = 0; i < localSize; ++i) {
      for (int j = 0; j < localSize; ++j) {
        m_data[start] += localJ(i, j);
        start += 1 - (int)m_cols;
      }
      start += localSize * ((int)m_cols - 1) + m_cols;
    }
  }

  // Add the local Jacobian, with a middle row and a middle column of zeros
  // inserted, to the banded matrix, top left corner at "start" position on the
  // diagonal. Parameter localSize must be an even number.
  template <IndexType localSize>
  void edgeStencilAdd(
      IndexType start,
      const Eigen::Matrix<ScalarT, localSize, localSize>& localJ) {
    start += m_cols * ku;

    for (int i = 0; i < localSize / 2; ++i) {
      for (int j = 0; j < localSize / 2; ++j) {
        m_data[start] += localJ(i, j);
        start += 1 - (int)m_cols;
      }
      start += 1 - m_cols;

      for (int j = localSize / 2; j < localSize; ++j) {
        m_data[start] += localJ(i, j);
        start += 1 - (int)m_cols;
      }
      start += (localSize + 1) * ((int)m_cols - 1) + (int)m_cols;
    }
    start += m_cols;

    for (int i = localSize / 2; i < localSize; ++i) {
      for (int j = 0; j < localSize / 2; ++j) {
        m_data[start] += localJ(i, j);
        start += 1 - (int)m_cols;
      }
      start += 1 - m_cols;

      for (int j = localSize / 2; j < localSize; ++j) {
        m_data[start] += localJ(i, j);
        start += 1 - (int)m_cols;
      }
      start += (localSize + 1) * ((int)m_cols - 1) + (int)m_cols;
    }
  }

  BandMatrix<ScalarT, kl, ku>& operator*=(const ScalarT multiplier) {
    for (typename std::vector<ScalarT>::iterator i = m_data.begin();
         i != m_data.end(); ++i)
      *i *= multiplier;

    return *this;
  }

  BandMatrix<ScalarT, kl, ku> operator+(
      const BandMatrix<ScalarT, kl, ku>& rhs) const {
    BandMatrix<ScalarT, kl, ku> copy(*this);
    return (copy += rhs);
  }

  BandMatrix<ScalarT, kl, ku> operator*(
      const BandMatrix<ScalarT, kl, ku>& rhs) const {
    BandMatrix<ScalarT, kl, ku> copy(*this);
    return (copy *= rhs);
  }

  BandMatrix<ScalarT, kl, ku>& operator+=(
      const BandMatrix<ScalarT, kl, ku>& rhs) {
    assert(rhs.rows() == rows() && rhs.cols() == cols());
    assert(m_size == rhs.m_size);

    VecXx::Map(&m_data[0], m_size) += VecXx::Map(&rhs.m_data[0], m_size);

    return *this;
  }

  template <typename MatrixT>
  BandMatrix<ScalarT, kl, ku>& operator+=(const MatrixT& rhs) {
    assert(rhs.rows() == rows() && rhs.cols() == cols());

    for (IndexType i = 0; i < rows(); ++i) {
      IndexType s = std::max(0, ((int)i) - (int)kl);
      IndexType e = std::min((int)cols(), (int)(i + ku + 1));

      for (IndexType j = s; j < e; ++j) {
        m_data[(ku + i - j) * m_cols + j] += rhs.coeff(i, j);
      }
    }

    return *this;
  }

  /* Very unoptimized Matrix/matrix product  - only for testing purposes*/
  BandMatrix<ScalarT, kl, ku>& operator*=(
      const BandMatrix<ScalarT, kl, ku>& rhs) {
    assert(rhs.rows() == rows() && rhs.cols() == cols());

    // FIXME
    // For now, do the stupidest possible thing

    MatXx A = MatXx::Zero(rows(), cols()), B = MatXx::Zero(rows(), cols());

    for (IndexType i = 0; i < rows(); ++i) {
      for (IndexType j = 0; j < cols(); ++j) {
        if (indicesValid(i, j))  // Should theorically not be used, but gcc
                                 // seems to use the wrong operator()
        {
          A(i, j) = (*this)(i, j);
          B(i, j) = rhs(i, j);
        }
      }
    }

    A *= B;

    for (IndexType i = 0; i < rows(); ++i) {
      for (IndexType j = 0; j < cols(); ++j) {
        if (indicesValid(i, j)) {
          (*this)(i, j) = A(i, j);
        }
      }
    }

    return *this;
  }

  void setZero() { m_data.assign(m_size, 0); }

  void multiply(VecXx& y, ScalarT s, const VecXx& x) const {
    assert(y.size() == m_rows);
    assert(x.size() == m_cols);

    for (int i = 0; i < m_rows; ++i) {
      int lower = std::max((int)i - (int)kl, 0);
      int upper = std::min((int)i + (int)ku + 1, (int)m_cols);
      const ScalarT* val = &m_data[(ku + i - lower) * m_cols + lower];
      ScalarT sum = 0.0;
      for (int j = lower; j < upper; ++j, val += (1 - (int)m_cols)) {
        sum += (*val) * x[j];
      }
      y[i] += s * sum;
    }
  }

  IndexType rows() const { return m_rows; }

  IndexType cols() const { return m_cols; }

  Scalar norm() {
    return MatXx::Map(m_data.data(), kl + ku + 1, m_cols).norm();
  }

  Scalar diagNorm() { return diagonal().array().abs().sum() / m_cols; }

  // NOTE: Assumes not symmetric if lower band is not same size as upper (could
  // be a bunch of zeros in bigger band)
  bool isApproxSymmetric(ScalarT eps) const {
    if (m_rows != m_cols || kl != ku) return false;

    for (int i = 0; i < m_rows; ++i)
      for (int j = i + 1; j <= i + ku; ++j)
        if (fabs((*this)(i, j) - (*this)(j, i)) > eps) return false;

    return true;
  }

  const std::vector<ScalarT>& getData() const { return m_data; }

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& m_rows;
    ar& m_cols;
    ar& m_size;
    ar& m_data;
  }

 private:
  bool indicesValid(IndexType r, IndexType c) const {
    return r < m_rows && c < m_cols && r - c <= kl && c - r <= ku;
  }

  IndexType m_rows;
  IndexType m_cols;
  IndexType m_size;             // Number of non-zero entries
  std::vector<ScalarT> m_data;  // Storage of non-zero entries
};

template <typename ScalarT, IndexType kl, IndexType ku>
std::ostream& operator<<(std::ostream& os,
                         const BandMatrix<ScalarT, kl, ku>& M) {
  os << '{';
  for (int i = 0; i < M.rows() - 1; i++) {
    os << '{';
    for (int j = 0; j < M.cols() - 1; j++) os << M(i, j) << ", ";
    os << M(i, M.cols() - 1) << "}, \n";
  }
  os << '{';
  for (int j = 0; j < M.cols() - 1; j++) os << M(M.rows() - 1, j) << ", ";
  os << M(M.rows() - 1, M.cols() - 1) << "}} \n";

  return os;
}

}  // namespace strandsim

#endif /* BANDMATRIX_HH_ */
