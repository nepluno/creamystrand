/**
 * \copyright 2014 Danny Kaufman
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Distances.hh"

#define ABSURDLY_LARGE_DISTANCE (1.e99)

namespace strandsim {

Vec3x ClosestPtPointSegment(const Vec3x& point, const Vec3x& first,
                            const Vec3x& last) {
  if (isClose(point, first)) return first;
  if (isClose(point, last)) return last;

  const Scalar vx = last[0] - first[0];
  const Scalar vy = last[1] - first[1];
  const Scalar vz = last[2] - first[2];
  const Scalar inv = 1. / (vy * vy + vx * vx + vz * vz);
  const Scalar dtp = vx * point[0] + vy * point[1] + vz * point[2];
  const Scalar flx = first[2] * last[1] - first[1] * last[2];
  const Scalar fly = first[0] * last[2] - first[2] * last[0];
  const Scalar flz = first[1] * last[0] - first[0] * last[1];

  // First project on the line
  Vec3x projonline((vz * fly - vy * flz + vx * dtp) * inv,
                   (vx * flz - vz * flx + vy * dtp) * inv,
                   (vy * flx - vx * fly + vz * dtp) * inv);

  // Check if we are outside the interval
  if ((projonline - last).dot(first - last) <= 0)
    return last;
  else if ((projonline - first).dot(last - first) <= 0)
    return first;
  else
    return projonline;
}

Vec3x ClosestPtPointSegment(bool& extremum, const Vec3x& point,
                            const Vec3x& first, const Vec3x& last) {
  if (isClose(point, first)) {
    extremum = true;
    return first;
  }
  if (isClose(point, last)) {
    extremum = true;
    return last;
  }

  const Scalar vx = last[0] - first[0];
  const Scalar vy = last[1] - first[1];
  const Scalar vz = last[2] - first[2];
  const Scalar inv = 1. / (vy * vy + vx * vx + vz * vz);
  const Scalar dtp = vx * point[0] + vy * point[1] + vz * point[2];
  const Scalar flx = first[2] * last[1] - first[1] * last[2];
  const Scalar fly = first[0] * last[2] - first[2] * last[0];
  const Scalar flz = first[1] * last[0] - first[0] * last[1];

  // First project on the line
  Vec3x projonline((vz * fly - vy * flz + vx * dtp) * inv,
                   (vx * flz - vz * flx + vy * dtp) * inv,
                   (vy * flx - vx * fly + vz * dtp) * inv);

  // Check if we are outside the interval
  if ((projonline - last).dot(first - last) <= 0) {
    extremum = true;
    return last;
  } else if ((projonline - first).dot(last - first) <= 0) {
    extremum = true;
    return first;
  } else {
    extremum = false;
    return projonline;
  }
}

// Adapted from Christer Ericson, "Real Time Collision Detection"
Vec3x ClosestPtPointTriangle(const Vec3x& p, const Vec3x& a, const Vec3x& b,
                             const Vec3x& c) {
  Vec3x result;

  // Check if P in vertex region outside A
  const Vec3x ab = b - a;
  const Vec3x ac = c - a;
  const Vec3x ap = p - a;
  double d1 = ab.dot(ap);
  double d2 = ac.dot(ap);
  if (d1 <= 0.0 && d2 <= 0.0) {
    result = a;  // barycentric coordinates (1,0,0)
    // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 1 p = " << p << "
    // a = " << a << " b = " << b << " c = " << c << " result = " << result <<
    // '\n';
    return result;
  }

  // Check if P in vertex region outside B
  const Vec3x bp = p - b;
  double d3 = ab.dot(bp);
  double d4 = ac.dot(bp);
  if (d3 >= 0.0 && d4 <= d3) {
    result = b;  // barycentric coordinates (0,1,0)
    // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 2 p = " << p << "
    // a = " << a << " b = " << b << " c = " << c << " result = " << result <<
    // '\n';
    return result;
  }

  // Check if P in edge region of AB, if so return projection of P onto AB
  double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    double v = d1 / (d1 - d3);
    result = a + v * ab;  // barycentric coordinates (1-v,v,0)
    // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 3 p = " << p << "
    // a = " << a << " b = " << b << " c = " << c << " result = " << result <<
    // '\n';
    return result;
  }

  // Check if P in vertex region outside C
  const Vec3x cp = p - c;
  double d5 = ab.dot(cp);
  double d6 = ac.dot(cp);
  if (d6 >= 0.0 && d5 <= d6) {
    result = c;  // barycentric coordinates (0,0,1)
    // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 4 p = " << p << "
    // a = " << a << " b = " << b << " c = " << c << " result = " << result <<
    // '\n';
    return result;
  }

  // Check if P in edge region of AC, if so return projection of P onto AC
  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    double w = d2 / (d2 - d6);
    result = a + w * ac;  // barycentric coordinates (1-w,0,w)
    // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 5 p = " << p << "
    // a = " << a << " b = " << b << " c = " << c << " result = " << result <<
    // '\n';
    return result;
  }

  // Check if P in edge region of BC, if so return projection of P onto BC
  double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    result = b + w * (c - b);  // barycentric coordinates (0,1-w,w)
    // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 6 p = " << p << "
    // a = " << a << " b = " << b << " c = " << c << " result = " << result <<
    // '\n';
    return result;
  }

  // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
  double denom = 1.0 / (va + vb + vc);
  double v = vb * denom;
  double w = vc * denom;
  result =
      a + ab * v + ac * w;  // = u*a + v*b + w*c, u = va * denom = 1.0f - v - w

  // TraceStream(g_log, "") << "ClosestPtPointTriangle: CASE 7 p = " << p << " a
  // = " << a << " b = " << b << " c = " << c << " result = " << result << '\n';

  return result;
}

// Adapted from Christer Ericson, "Real Time Collision Detection"
// Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
// S2(t)=P2+t*(Q2-P2), returning s and t. Function result is squared
// distance between between S1(s) and S2(t).
// TODO: Explore behavior in degenerate case more closely.
double ClosestPtSegmentSegment(const Vec3x& p1, const Vec3x& q1,
                               const Vec3x& p2, const Vec3x& q2, double& s,
                               double& t, Vec3x& c1, Vec3x& c2) {
  double EPSILON = 1.0e-12;

  Vec3x d1 = q1 - p1;  // Direction vector of segment S1
  Vec3x d2 = q2 - p2;  // Direction vector of segment S2
  Vec3x r = p1 - p2;
  double a = d1.dot(d1);  // Squared length of segment S1, always nonnegative
  double e = d2.dot(d2);  // Squared length of segment S2, always nonnegative
  double f = d2.dot(r);

  // Check if either or both segments degenerate into points
  if (a <= EPSILON && e <= EPSILON) {
    // Both segments degenerate into points
    s = t = 0.0;
    c1 = p1;
    c2 = p2;
    return (c1 - c2).dot(c1 - c2);
  }
  if (a <= EPSILON) {
    // First segment degenerates into a point
    s = 0.0;
    t = f / e;  // s = 0 => t = (b*s + f) / e = f / e
    t = clamp(t, 0.0, 1.0);
  } else {
    double c = d1.dot(r);
    if (e <= EPSILON) {
      // Second segment degenerates into a point
      t = 0.0;
      s = clamp(-c / a, 0.0, 1.0);  // t = 0 => s = (b*t - c) / a = -c / a
    } else {
      // The general nondegenerate case starts here
      double b = d1.dot(d2);
      double denom = a * e - b * b;  // Always nonnegative

      // If segments not parallel, compute closest point on L1 to L2, and
      // clamp to segment S1. Else pick arbitrary s (here 0)
      if (denom != 0.0) {
        s = clamp((b * f - c * e) / denom, 0.0, 1.0);
      } else
        s = 0.0;

      // Compute point on L2 closest to S1(s) using
      // t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f) / e
      t = (b * s + f) / e;

      // If t in [0,1] done. Else clamp t, recompute s for the new value
      // of t using s = Dot((P2+D2*t)-P1,D1) / Dot(D1,D1)= (t*b - c) / a
      // and clamp s to [0, 1]
      if (t < 0.0) {
        t = 0.0;
        s = clamp(-c / a, 0.0, 1.0);
      } else if (t > 1.0) {
        t = 1.0;
        s = clamp((b - c) / a, 0.0, 1.0);
      }
    }
  }

  c1 = p1 + d1 * s;
  c2 = p2 + d2 * t;
  return (c1 - c2).dot(c1 - c2);
}

// Adapted from Christer Ericson, "Real Time Collision Detection"
// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
void computeBarycentricCoordinates(const Vec3x& a, const Vec3x& b,
                                   const Vec3x& c, const Vec3x& p, double& u,
                                   double& v, double& w) {
  const Vec3x ab = b - a;
  const Vec3x ac = c - a;
  const Vec3x ap = p - a;
  const double abab = ab.dot(ab);
  const double abac = ab.dot(ac);
  const double acac = ac.dot(ac);
  const double apab = ap.dot(ab);
  const double apac = ap.dot(ac);
  const double denom = abab * acac - abac * abac;

  if (fabs(denom) > std::numeric_limits<double>::epsilon()) {
    v = (acac * apab - abac * apac) / denom;
    w = (abab * apac - abac * apab) / denom;
    u = 1. - v - w;
  } else  // Points are aligned: next best thing,
          // compute the barycentric coordinates among the two closest points to
          // the projection of p on the line
  {
    const Vec3x cb = b - c;
    const Vec3x cp = p - c;
    const double cbcb = cb.dot(cb);
    const double cpcb = cp.dot(cb);
    const double sab = apab / abab;
    const double sac = apac / acac;
    const double scb = cpcb / cbcb;
    const double xab = sab * (1 - sab);
    const double xac = sac * (1 - sac);
    const double xcb = scb * (1 - scb);

    if (acac > abab && acac > cbcb) {
      if (xab > xcb) {
        u = 1. - sab;
        v = sab;
        w = 0.;
      } else {
        u = 0;
        v = scb;
        w = 1 - scb;
      }
    } else if (abab > cbcb) {
      if (xac > xcb) {
        u = 1. - sac;
        v = 0.;
        w = sac;
      } else {
        u = 0.;
        v = scb;
        w = 1 - scb;
      }
    } else {
      if (xac > xab) {
        u = 1 - sac;
        v = 0.;
        w = sac;
      } else {
        u = 1 - sab;
        v = sab;
        w = 0.;
      }
    }
  }
}

// Utility finctions for rect/rect distance -- may be exposed later

void assignIfCloser(Scalar tentativeDist, const Vec3x& tentativeP,
                    Scalar& minDist, Vec3x& CP) {
  if (tentativeDist < minDist) {
    minDist = tentativeDist;
    CP = tentativeP;
  }
}

void assignIfCloser(Scalar tentativeDist, const Vec3x& tentativeP1,
                    const Vec3x& tentativeP2, Scalar& minDist, Vec3x& CP1,
                    Vec3x& CP2) {
  if (tentativeDist < minDist) {
    minDist = tentativeDist;
    CP1 = tentativeP1;
    CP2 = tentativeP2;
  }
}

Scalar sqDistPointSeg(const Vec3x& P, const Vec3x& S, const Vec3x& D,
                      Vec3x& CP) {
  CP = ClosestPtPointSegment(P, S, S + D);
  return (P - CP).squaredNorm();
}

Scalar sqDistSegSeg(const Vec3x& S1, const Vec3x& D1, const Vec3x& S2,
                    const Vec3x& D2, Vec3x& CP1, Vec3x& CP2) {
  Scalar s, t;
  return ClosestPtSegmentSegment(S1, S1 + D1, S2, S2 + D2, s, t, CP1, CP2);

  //    Scalar s,t ;
  //    Scalar sqDist = SquareDistSegmentToSegment< Vec3x, Scalar, Vec3x >
  //            ( S1, S1+D1, S2, S2+D2, s, t ) ;

  //    if( s < 0 || t < 0 ) return ABSURDLY_LARGE_DISTANCE ;

  //    CP1 = S1 + s*D1 ;
  //    CP2 = S2 + t*D2 ;

  //    return sqDist ;
}

Scalar sqDistPointRect(const Vec3x& P, const Vec3x& C, const Vec3x& x,
                       const Vec3x& y, Vec3x& CP) {
  const Vec3x& D = P - C;

  const Scalar xn2 = x.squaredNorm();
  const Scalar yn2 = y.squaredNorm();

  const Scalar Dx = D.dot(x);
  const Scalar Dy = D.dot(y);

  if (std::fabs(Dx) < xn2 && std::fabs(Dy) < yn2) {
    // Point projects inside rectangle
    CP = C + Dx * x / xn2 + Dy * y / xn2;

    return (P - CP).squaredNorm();
  }

  Scalar minDist = ABSURDLY_LARGE_DISTANCE;
  Vec3x tCP;

  assignIfCloser(sqDistPointSeg(P, C - x - y, 2 * x, tCP), tCP, minDist, CP);
  assignIfCloser(sqDistPointSeg(P, C - x + y, 2 * x, tCP), tCP, minDist, CP);
  assignIfCloser(sqDistPointSeg(P, C - y - x, 2 * y, tCP), tCP, minDist, CP);
  assignIfCloser(sqDistPointSeg(P, C - y + x, 2 * y, tCP), tCP, minDist, CP);

  return minDist;
}

bool intersectionSegmentRectangle(const Vec3x& s_edge_0, const Vec3x& s_edge_1,
                                  const Vec3x& r_center,
                                  const Vec3x& r_extent_1,
                                  const Vec3x& r_extent_2, Scalar& t) {
  const Vec3x& n = r_extent_1.cross(r_extent_2);
  const Vec3x& D0 = s_edge_0 - r_center;
  const Vec3x& D1 = s_edge_1 - r_center;

  const Scalar D0n = D0.dot(n);
  const Scalar D1n = D1.dot(n);

  // both ends are on the same side ; or segment coplanar to rectangle
  if (D0n * D1n > -SMALL_NUMBER<Scalar>()) return false;

  // time of intersection with plane
  const Scalar ti = D0n / (D0n - D1n);
  const Vec3x& Pi = D0 + ti * (D1 - D0);

  const Scalar xn2 = r_extent_1.squaredNorm();
  const Scalar yn2 = r_extent_2.squaredNorm();

  const Scalar Dx = Pi.dot(r_extent_1);
  const Scalar Dy = Pi.dot(r_extent_2);

  if (std::fabs(Dx) < xn2 && std::fabs(Dy) < yn2) {
    t = ti;
    return true;
  }

  return false;
}

Scalar sqDistSegRect(const Vec3x& S, const Vec3x& D, const Vec3x& C,
                     const Vec3x& x, const Vec3x& y, Vec3x& CP1, Vec3x& CP2) {
  Scalar t;

  if (intersectionSegmentRectangle(S, S + D, C, x, y, t)) {
    // Segment and rectangle are intersecting
    CP1 = CP2 = S + t * D;
    return 0;
  }

  Scalar minDist = ABSURDLY_LARGE_DISTANCE;
  Vec3x P1, P2;

  // Point Rect
  assignIfCloser(sqDistPointRect(S, C, x, y, P2), S, P2, minDist, CP1, CP2);
  assignIfCloser(sqDistPointRect(S + D, C, x, y, P2), S + D, P2, minDist, CP1,
                 CP2);

  // Seg Seg
  assignIfCloser(sqDistSegSeg(S, D, C - x - y, 2 * x, P1, P2), P1, P2, minDist,
                 CP1, CP2);
  assignIfCloser(sqDistSegSeg(S, D, C - x + y, 2 * x, P1, P2), P1, P2, minDist,
                 CP1, CP2);
  assignIfCloser(sqDistSegSeg(S, D, C - y - x, 2 * y, P1, P2), P1, P2, minDist,
                 CP1, CP2);
  assignIfCloser(sqDistSegSeg(S, D, C - y + x, 2 * y, P1, P2), P1, P2, minDist,
                 CP1, CP2);

  return minDist;
}

// Requires xi . yi = 0
Scalar SquareDistRectangleToRectangle(const Vec3x& C1, const Vec3x& x1,
                                      const Vec3x& y1, const Vec3x& C2,
                                      const Vec3x& x2, const Vec3x& y2,
                                      Vec3x& CP1, Vec3x& CP2) {
  Scalar minDist = ABSURDLY_LARGE_DISTANCE;
  Vec3x P1, P2;

  // Small dist from each edge of R1 to R2 ( and vice versa )
  assignIfCloser(sqDistSegRect(C1 + x1 - y1, 2 * y1, C2, x2, y2, P1, P2), P1,
                 P2, minDist, CP1, CP2);
  if (minDist == 0.) return minDist;
  assignIfCloser(sqDistSegRect(C1 - x1 - y1, 2 * y1, C2, x2, y2, P1, P2), P1,
                 P2, minDist, CP1, CP2);
  if (minDist == 0.) return minDist;
  assignIfCloser(sqDistSegRect(C1 + y1 - x1, 2 * x1, C2, x2, y2, P1, P2), P1,
                 P2, minDist, CP1, CP2);
  if (minDist == 0.) return minDist;
  assignIfCloser(sqDistSegRect(C1 - y1 - x1, 2 * x1, C2, x2, y2, P1, P2), P1,
                 P2, minDist, CP1, CP2);
  if (minDist == 0.) return minDist;
  assignIfCloser(sqDistSegRect(C2 + x2 - y2, 2 * y2, C1, x1, y1, P1, P2), P1,
                 P2, minDist, CP2, CP1);
  if (minDist == 0.) return minDist;
  assignIfCloser(sqDistSegRect(C2 - x2 - y2, 2 * y2, C1, x1, y1, P1, P2), P1,
                 P2, minDist, CP2, CP1);
  if (minDist == 0.) return minDist;
  assignIfCloser(sqDistSegRect(C2 + y2 - x2, 2 * x2, C1, x1, y1, P1, P2), P1,
                 P2, minDist, CP2, CP1);
  if (minDist == 0.) return minDist;
  assignIfCloser(sqDistSegRect(C2 - y2 - x2, 2 * x2, C1, x1, y1, P1, P2), P1,
                 P2, minDist, CP2, CP1);

  return minDist;
}

}  // namespace strandsim
