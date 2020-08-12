/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ELASTICSTRANDUTILS_HH_
#define ELASTICSTRANDUTILS_HH_
#ifdef WIN32
#include <boost/math/special_functions/fpclassify.hpp>
#endif
#include <limits>
#include <stdexcept>

#include "Definitions.hh"

namespace strandsim {

/**
 * \brief Tests if a matrix is symmetric
 */
template <typename MatrixT>
inline bool isSymmetric(const MatrixT& A) {
  for (int i = 0; i < A.rows(); ++i)
    for (int j = i + 1; j < A.cols(); ++j)
      if (!isSmall(A(i, j) - A(j, i))) {
        std::cerr << "isSymmetric failing by " << fabs(A(i, j) - A(j, i))
                  << '\n';
        return false;
      }
  return true;
}

/**
 * \brief Projects v on the plane normal to n and normalize it
 */
inline void orthoNormalize(Vec3x& v, const Vec3x& n) {
  assert(isApproxUnit(n));

  v -= v.dot(n) * n;
  v.normalize();
}

/**
 * \brief Finds an arbitrary unit vector orthogonal to u
 *
 * \param [in] u vector in dimension n
 * \return An arbitrary vector orthogonal to u
 */
template <int n>
inline Eigen::Matrix<Scalar, n, 1> findNormal(
    const Eigen::Matrix<Scalar, n, 1>& u) {
  assert(u.norm() != 0);
  Eigen::Matrix<Scalar, n, 1> v = Eigen::Matrix<Scalar, n, 1>::Zero();

  int maxCoordinate = 0;
  for (int i = 0; i < n; ++i) {
    if (isSmall(u[i])) {
      v[i] = 1;
      goto finished;
    }
    if (fabs(u[i]) > fabs(u[maxCoordinate])) maxCoordinate = i;
  }
  {
    const int otherCoordinate = (maxCoordinate + 1) % n;
    v[otherCoordinate] = u[maxCoordinate];
    v[maxCoordinate] = -u[otherCoordinate];
  }
  v.normalize();

finished:
  assert(isSmall(u.dot(v)));

  return v;
}

/**
 * \brief Computes the signed angle from one vector to another.
 *
 * \param [in] u first input vector
 * \param [in] v second input vector
 * \param [in] n orientation vector
 *
 * The sign of the angle is positive if u.cross(v) is in the same half-space as
 * n.
 */
inline Scalar signedAngle(const Vec3x& u, const Vec3x& v, const Vec3x& n) {
  Vec3x w = u.cross(v);
  Scalar angle = atan2(w.norm(), u.dot(v));
  if (n.dot(w) < 0) return -angle;
  return angle;
}

/**
 * \brief Rotates a vector
 *
 * \param [in] v vector to rotate
 * \param [in] z normalized vector on the rotation axis
 * \param [in] theta rotation angle
 *
 */
template <typename ScalarT>
inline void rotateAxisAngle(typename Eigen::Matrix<ScalarT, 3, 1>& v,
                            const typename Eigen::Matrix<ScalarT, 3, 1>& z,
                            const ScalarT theta) {
  assert(isApproxUnit(z));

  if (theta == 0) return;

  const ScalarT c = cos(theta);
  const ScalarT s = sin(theta);

  v = c * v + s * z.cross(v) + z.dot(v) * (1.0 - c) * z;
}

/**
 * \brief Outer product of two vectors.
 */
template <int n>
inline Eigen::Matrix<Scalar, n, n> outerProd(
    const Eigen::Matrix<Scalar, n, 1>& a,
    const Eigen::Matrix<Scalar, n, 1>& b) {
  return a * b.transpose();
}

/**
 * \brief Outer product of two vectors.
 */
template <int n>
inline Eigen::Matrix<Scalar, n, n> symPart(
    const Eigen::Matrix<Scalar, n, n>& A) {
  return (A + A.transpose()) * 0.5;
}

/**
 * \brief Computes u^T B v, assuming B is symmetric 2x2 and u, v are 2x1
 * vectors.
 */
inline Scalar innerBProduct(const Mat2x& B, const Vec2x& u, const Vec2x& v) {
  assert(isSymmetric(B));

  //    return u[0] * ( B( 0, 0 ) * v[0] + B( 0, 1 ) * v[1] )
  //            + u[1] * ( B( 1, 0 ) * v[0] + B( 1, 1 ) * v[1] ); // Bad
  return B(0, 0) * u[0] * v[0] + B(0, 1) * (u[0] * v[1] + u[1] * v[0]) +
         B(1, 1) * u[1] * v[1];  // Good
}

/**
 * \brief Computes Q B Q^T, assuming B is symmetric 2x2 and Q is nx2. The result
 * is then (exactly) symmetric nxn.
 */
template <int n>
inline void symBProduct(Eigen::Matrix<Scalar, n, n>& result, const Mat2x& B,
                        const Eigen::Matrix<Scalar, n, 2>& Q) {
  assert(isSymmetric(B));

  for (int i = 0; i < n; ++i) {
    const Vec2x& Qrow_i = Q.row(i);
    result(i, i) = innerBProduct(B, Qrow_i, Qrow_i);
    for (int j = 0; j < i; ++j)
      result(i, j) = result(j, i) = innerBProduct(B, Qrow_i, Q.row(j));
  }
}

/**
 * \brief Computes Q B Q^T, assuming B is symmetric 2x2 and Q is nx2. The result
 * is then (exactly) symmetric nxn.
 */
template <int n>
inline void symBProductAdd(Eigen::Matrix<Scalar, n, n>& result, const Mat2x& B,
                           const Eigen::Matrix<Scalar, n, 2>& Q) {
  assert(isSymmetric(B));

  for (int i = 0; i < n; ++i) {
    const Vec2x& Qrow_i = Q.row(i);
    result(i, i) += innerBProduct(B, Qrow_i, Qrow_i);
    for (int j = 0; j < i; ++j) {
      const Scalar val = innerBProduct(B, Qrow_i, Q.row(j));
      result(i, j) += val;
      result(j, i) += val;
    }
  }
}
/**
 * \brief Angular interpolation by a factor t between two vectors v0 and v1
 * @param v0
 * @param v1
 * @param t
 * @return
 */
template <typename VectorT, typename ScalarT>
inline VectorT vectorSlerp(const VectorT& v0, const VectorT& v1, ScalarT t) {
  const ScalarT angle = std::acos(clamp(v0.dot(v1), -1., 1.));
  if (isSmall(angle)) return v0;
  const ScalarT invSin = 1. / std::sin(angle);
  return invSin * (std::sin((1. - t) * angle) * v0 + std::sin(t * angle) * v1);
}

/**
 * \brief Applies the (projective space) transform M to x identified to (x,1)
 * @param M
 * @param x
 * @return
 */
template <typename ScalarT>
Eigen::Matrix<ScalarT, 3, 1> transformPoint(
    const Eigen::Matrix<ScalarT, 4, 4>& M,
    const Eigen::Matrix<ScalarT, 3, 1>& x) {
  return M.template block<3, 3>(0, 0) * x + M.template block<3, 1>(0, 3);
}

/**
 \return an angle theta in [-pi, pi] such that theta == \p angle mod (2pi)
 */
inline Scalar clamp2Pi(Scalar angle) {
  Scalar theta = angle;
  while (theta > M_PI) theta -= 2. * M_PI;
  while (theta <= -M_PI) theta += 2. * M_PI;

  return theta;
}

}  // namespace strandsim

#endif /* ELASTICSTRANDUTILS_HH_ */
