/**
 * \copyright 2014 Danny Kaufman
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSIM_GEOMETRY_HH
#define STRANDSIM_GEOMETRY_HH

#include "../Core/Definitions.hh"

namespace strandsim {

// Closest point on a segment to a vertex
Vec3x ClosestPtPointSegment(const Vec3x& point, const Vec3x& first,
                            const Vec3x& last);

Vec3x ClosestPtPointSegment(bool& extremum, const Vec3x& point,
                            const Vec3x& first, const Vec3x& last);

// Closest point on a triangle to a vertex
Vec3x ClosestPtPointTriangle(const Vec3x& p, const Vec3x& a, const Vec3x& b,
                             const Vec3x& c);

// Computes the squared distance between and closest points of two edges.
double ClosestPtSegmentSegment(const Vec3x& p1, const Vec3x& q1,
                               const Vec3x& p2, const Vec3x& q2, double& s,
                               double& t, Vec3x& c1, Vec3x& c2);

// Computes the barycentric coordiantes of a point wrt a triangle
void computeBarycentricCoordinates(const Vec3x& a, const Vec3x& b,
                                   const Vec3x& c, const Vec3x& p, double& u,
                                   double& v, double& w);

// Ci: center of rectangle i ; xi and yi : extents on boths axis
// CPi: closest point on rectangle i
Scalar SquareDistRectangleToRectangle(const Vec3x& C1, const Vec3x& x1,
                                      const Vec3x& y1, const Vec3x& C2,
                                      const Vec3x& x2, const Vec3x& y2,
                                      Vec3x& CP1, Vec3x& CP2);

// /!\ Will return false if segment and rectangle are coplanar
bool intersectionSegmentRectangle(const Vec3x& s_edge_0, const Vec3x& s_edge_1,
                                  const Vec3x& r_center,
                                  const Vec3x& r_extent_1,
                                  const Vec3x& r_extent_2, Scalar& t);

// Templated segment/point and segment/segment distances

template <typename VecT, typename ScalarT>
ScalarT SquareDistPointToSegment(const VecT& point, const VecT& first,
                                 const VecT& last) {
  const VecT dfirst = (point - first);
  const ScalarT sqDistFirst =
      dfirst[0] * dfirst[0] + dfirst[1] * dfirst[1] + dfirst[2] * dfirst[2];
  if (isSmall(sqDistFirst)) return sqDistFirst;

  const ScalarT vx = last[0] - first[0];
  const ScalarT vy = last[1] - first[1];
  const ScalarT vz = last[2] - first[2];
  // Squared norm of segment
  const ScalarT len2 = (vy * vy + vx * vx + vz * vz);
  // Abscissa of proj on line
  const ScalarT dtp = (vx * dfirst[0] + vy * dfirst[1] + vz * dfirst[2]) / len2;

  if (dtp <= 0) {
    return sqDistFirst;
  }
  const ScalarT dtp2 = dtp * dtp;

  if (dtp2 >= len2) {
    const VecT dlast = (point - last);
    const ScalarT sqDistLast =
        dlast[0] * dlast[0] + dlast[1] * dlast[1] + dlast[2] * dlast[2];

    return sqDistLast;
  }
  return sqDistFirst - dtp2;
}

template <typename VecT, typename ScalarT, typename InVecT>
ScalarT SquareDistSegmentToSegment(const Eigen::MatrixBase<InVecT>& p0,
                                   const Eigen::MatrixBase<InVecT>& p1,
                                   const Eigen::MatrixBase<InVecT>& q0,
                                   const Eigen::MatrixBase<InVecT>& q1,
                                   ScalarT& s, ScalarT& t) {
  const VecT dp = p1 - p0;  // Direction vector of segment S1
  const VecT dq = q1 - q0;  // Direction vector of segment S2
  const VecT r = p0 - q0;

  const ScalarT a =
      dp.dot(dp);  // Squared length of segment S1, always nonnegative
  const ScalarT e =
      dq.dot(dq);  // Squared length of segment S2, always nonnegative
  const ScalarT f = dq.dot(r);

  const ScalarT c = dp.dot(r);
  const ScalarT b = dp.dot(dq);

  const ScalarT denom = a * e - b * b;

  if (strandsim::isSmall(denom))  // parallel
  {
    const ScalarT s0 = -c / a;
    const ScalarT s1 = (b - c) / a;

    s = 0.5;  // FIXME
    t = 0.5;  // FIXME

    if (s0 < 0) {
      s = 0;
      if (s1 < s0) {
        t = 0;
        return r.squaredNorm();
      } else if (s1 < 0) {
        t = 1;
        return (p0 - q1).squaredNorm();
      }
    } else if (s0 > 1) {
      s = 1;
      if (s1 > s0) {
        t = 0;
        return (p1 - q0).squaredNorm();
      } else if (s1 > 1) {
        t = 1;
        return (p1 - q1).squaredNorm();
      }
    }

    return (r - c * dp / a).squaredNorm();
  }

  const ScalarT s_ = (b * f - c * e) / denom;
  const ScalarT t_ = (b * s_ + f) / e;

  s = strandsim::clamp<ScalarT>(s_, 0, 1);
  t = strandsim::clamp<ScalarT>(t_, 0, 1);

  return (p0 + s * dp - q0 - t * dq).squaredNorm();
}

template <typename VecT, typename ScalarT, typename InVecT>
ScalarT DistSegmentToSegment(const Eigen::MatrixBase<InVecT>& p0,
                             const Eigen::MatrixBase<InVecT>& p1,
                             const Eigen::MatrixBase<InVecT>& q0,
                             const Eigen::MatrixBase<InVecT>& q1) {
  ScalarT s, t;
  std::sqrt(
      SquareDistSegmentToSegment<VecT, ScalarT, InVecT>(p0, p1, q0, q1, s, t));
}

template <typename VecT, typename ScalarT, typename InVecT>
ScalarT HausdorffDistSegmentToSegment(const Eigen::MatrixBase<InVecT>& p0,
                                      const Eigen::MatrixBase<InVecT>& p1,
                                      const Eigen::MatrixBase<InVecT>& q0,
                                      const Eigen::MatrixBase<InVecT>& q1) {
  const VecT dp = p1 - p0;  // Direction vector of segment S1
  const VecT dq = q1 - q0;  // Direction vector of segment S2
  const VecT r = p0 - q0;

  const ScalarT a =
      dp.dot(dp);  // Squared length of segment S1, always nonnegative
  const ScalarT e =
      dq.dot(dq);  // Squared length of segment S2, always nonnegative

  const ScalarT f = dq.dot(r);
  const ScalarT c = dp.dot(r);
  const ScalarT b = dp.dot(dq);

  const ScalarT s0 = strandsim::clamp(-c / a, 0.f, 1.f);
  const ScalarT s1 = strandsim::clamp((b - c) / a, 0.f, 1.f);

  const ScalarT t0 = strandsim::clamp(f / e, 0.f, 1.f);
  const ScalarT t1 = strandsim::clamp((b + f) / e, 0.f, 1.f);

  return std::sqrt(std::max(std::max((p0 + s0 * dp - q0).squaredNorm(),
                                     (p0 + s1 * dp - q1).squaredNorm()),
                            std::max((q0 + t0 * dq - p0).squaredNorm(),
                                     (q0 + t1 * dq - p1).squaredNorm())));
}

}  // namespace strandsim

#endif  // STRANDSIM_GEOMETRY_HH
