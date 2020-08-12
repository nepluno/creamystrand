/**
 * \copyright 2010 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "CollisionUtils.hh"

#include <math.h>

#include "../Core/ElasticStrand.hh"
#include "../Dynamic/ImplicitStepper.hh"
#include "../Dynamic/StrandDynamicTraits.hh"
#include "../Utils/Distances.hh"
#include "../Utils/MathUtilities.hh"
#include "OrientedBoundingBox.hh"
#include "ProximityCollision.hh"
#ifndef M_PI_2
#define M_PI_2 M_PI / 2
#endif
#define COS_PARALLEL_ENOUGH 0.86602540378  // cos( Pi/6 )
namespace strandsim {

// Adapted from some code on Robert Bridson's website, I believe

void addUnique(std::vector<double>& a, double e) {
  for (int i = 0; i < a.size(); ++i)
    if (a[i] == e) return;
  a.push_back(e);
}

void addUnique(double* a, unsigned& a_size, double e) {
  for (int i = 0; i < a_size; ++i)
    if (a[i] == e) return;
  a[a_size++] = e;
}

inline void compare_and_swap(double& a, double& b) {
  if (a > b) std::swap(a, b);
}

//! Sorting networks for a_size <= 4
void sort(double* a, unsigned a_size) {
  switch (a_size) {
    case 4:
      compare_and_swap(a[0], a[2]);
      compare_and_swap(a[1], a[3]);
      compare_and_swap(a[0], a[1]);
      compare_and_swap(a[2], a[3]);
      compare_and_swap(a[1], a[2]);
      break;
    case 3:
      compare_and_swap(a[0], a[1]);
      compare_and_swap(a[0], a[2]);
      compare_and_swap(a[1], a[2]);
      break;
    case 2:
      compare_and_swap(a[0], a[1]);
      break;
    default:
      break;
  }
}

double triple(const Vec3x& a, const Vec3x& b, const Vec3x& c) {
  return a[0] * (b[1] * c[2] - b[2] * c[1]) +
         a[1] * (b[2] * c[0] - b[0] * c[2]) +
         a[2] * (b[0] * c[1] - b[1] * c[0]);
}

double signed_volume(const Vec3x& x0, const Vec3x& x1, const Vec3x& x2,
                     const Vec3x& x3) {
  // Equivalent to triple(x1-x0, x2-x0, x3-x0), six times the signed volume of
  // the tetrahedron. But, for robustness, we want the result (up to sign) to be
  // independent of the ordering. And want it as accurate as possible... But all
  // that stuff is hard, so let's just use the common assumption that all
  // coordinates are >0, and do something reasonably accurate in fp.

  // This formula does almost four times too much multiplication, but if the
  // coordinates are non-negative it suffers in a minimal way from cancellation
  // error.
  return (x0[0] * (x1[1] * x3[2] + x3[1] * x2[2] + x2[1] * x1[2]) +
          x1[0] * (x2[1] * x3[2] + x3[1] * x0[2] + x0[1] * x2[2]) +
          x2[0] * (x3[1] * x1[2] + x1[1] * x0[2] + x0[1] * x3[2]) +
          x3[0] * (x1[1] * x2[2] + x2[1] * x0[2] + x0[1] * x1[2]))

         - (x0[0] * (x2[1] * x3[2] + x3[1] * x1[2] + x1[1] * x2[2]) +
            x1[0] * (x3[1] * x2[2] + x2[1] * x0[2] + x0[1] * x3[2]) +
            x2[0] * (x1[1] * x3[2] + x3[1] * x0[2] + x0[1] * x1[2]) +
            x3[0] * (x2[1] * x1[2] + x1[1] * x0[2] + x0[1] * x2[2]));
}

// All roots returned in interval [0,1]. Assumed geometry followed a linear
// trajectory between x and xnew.
void getCoplanarityTimes(const Vec3x& x0, const Vec3x& x1, const Vec3x& x2,
                         const Vec3x& x3, const Vec3x& xnew0,
                         const Vec3x& xnew1, const Vec3x& xnew2,
                         const Vec3x& xnew3, double* times, double* errors,
                         unsigned& num_times) {
  const double tol = 1e-8;
  num_times = 0;

  // cubic coefficients, A*t^3+B*t^2+C*t+D (for t in [0,1])
  const Vec3x x03 = x0 - x3;
  const Vec3x x13 = x1 - x3;
  const Vec3x x23 = x2 - x3;
  const Vec3x v03 = (xnew0 - xnew3) - x03;
  const Vec3x v13 = (xnew1 - xnew3) - x13;
  const Vec3x v23 = (xnew2 - xnew3) - x23;

  double A = triple(v03, v13, v23);
  double B =
      triple(x03, v13, v23) + triple(v03, x13, v23) + triple(v03, v13, x23);
  double C =
      triple(x03, x13, v23) + triple(x03, v13, x23) + triple(v03, x13, x23);
  double D = triple(x03, x13, x23);

  const double convergence_tol =
      tol * (std::fabs(A) + std::fabs(B) + std::fabs(C) + std::fabs(D));

  // find intervals to check, or just solve it if it reduces to a quadratic
  // =============================
  double interval_times[4];
  int interval_times_size = 0;

  double discriminant =
      B * B - 3 * A * C;  // of derivative of cubic, 3*A*t^2+2*B*t+C, divided by
                          // 4 for convenience
  if (discriminant <= 0) {  // monotone cubic: only one root in [0,1] possible
    // so we just
    interval_times[0] = 0;
    interval_times[1] = 1;
    interval_times_size = 2;
  } else {         // positive discriminant, B!=0
    if (A == 0) {  // the cubic is just a quadratic, B*t^2+C*t+D
                   // ========================================
      discriminant = C * C - 4 * B * D;  // of the quadratic
      if (discriminant <= 0) {
        double t = -C / (2 * B);
        if (t >= -tol && t <= 1 + tol) {
          t = clamp(t, 0., 1.);
          double val = std::fabs(signed_volume(
              (1 - t) * x0 + t * xnew0, (1 - t) * x1 + t * xnew1,
              (1 - t) * x2 + t * xnew2, (1 - t) * x3 + t * xnew3));
          if (val < convergence_tol) {
            times[num_times++] = t;
          }
        }
      } else {  // two separate real roots
        double t0, t1;
        if (C > 0)
          t0 = (-C - std::sqrt(discriminant)) / (2 * B);
        else
          t0 = (-C + std::sqrt(discriminant)) / (2 * B);
        t1 = D / (B * t0);
        if (t1 < t0) std::swap(t0, t1);
        if (t0 >= -tol && t0 <= 1 + tol) {
          times[num_times++] = clamp(t0, 0., 1.);
        }
        if (t1 >= -tol && t1 <= 1 + tol) {
          addUnique(times, num_times, clamp(t1, 0., 1.));
        }
      }

      if (errors) {
        for (int i = 0; i < num_times; ++i) {
          double ti = times[i];
          double val = std::fabs(signed_volume(
              (1 - ti) * x0 + ti * xnew0, (1 - ti) * x1 + ti * xnew1,
              (1 - ti) * x2 + ti * xnew2, (1 - ti) * x3 + ti * xnew3));
          errors[i] = val;
        }
      }

      return;
    } else {  // cubic is not monotone: divide up [0,1] accordingly
              // =====================================
      double t0, t1;
      if (B > 0)
        t0 = (-B - std::sqrt(discriminant)) / (3 * A);
      else
        t0 = (-B + std::sqrt(discriminant)) / (3 * A);
      t1 = C / (3 * A * t0);
      if (t1 < t0) std::swap(t0, t1);

      interval_times[interval_times_size++] = 0;
      if (t0 > 0 && t0 < 1) interval_times[interval_times_size++] = t0;
      if (t1 > 0 && t1 < 1) interval_times[interval_times_size++] = t1;

      interval_times[interval_times_size++] = 1;
    }
  }

  // look for roots in indicated intervals
  // ============================================================== evaluate
  // coplanarity more accurately at each endpoint of the intervals
  std::vector<double> interval_values(interval_times_size);
  for (int i = 0; i < interval_times_size; ++i) {
    double t = interval_times[i];
    interval_values[i] =
        signed_volume((1 - t) * x0 + t * xnew0, (1 - t) * x1 + t * xnew1,
                      (1 - t) * x2 + t * xnew2, (1 - t) * x3 + t * xnew3);
  }
  // first look for interval endpoints that are close enough to zero, without a
  // sign change
  for (int i = 0; i < interval_times_size; ++i) {
    if (interval_values[i] == 0) {
      times[num_times++] = interval_times[i];
    } else if (std::fabs(interval_values[i]) < convergence_tol) {
      if ((i == 0 || (interval_values[i - 1] >= 0 && interval_values[i] >= 0) ||
           (interval_values[i - 1] <= 0 && interval_values[i] <= 0)) &&
          (i == interval_times_size - 1 ||
           (interval_values[i + 1] >= 0 && interval_values[i] >= 0) ||
           (interval_values[i + 1] <= 0 && interval_values[i] <= 0))) {
        times[num_times++] = interval_times[i];
      }
    }
  }
  // and then search in intervals with a sign change
  for (int i = 1; i < interval_times_size; ++i) {
    double tlo = interval_times[i - 1], thi = interval_times[i], tmid;
    double vlo = interval_values[i - 1], vhi = interval_values[i], vmid;
    if ((vlo < 0 && vhi > 0) || (vlo > 0 && vhi < 0)) {
      // start off with secant approximation (in case the cubic is actually
      // linear)
      double alpha = vhi / (vhi - vlo);
      tmid = alpha * tlo + (1 - alpha) * thi;
      for (int iteration = 0; iteration < 50; ++iteration) {
        vmid = signed_volume(
            (1 - tmid) * x0 + tmid * xnew0, (1 - tmid) * x1 + tmid * xnew1,
            (1 - tmid) * x2 + tmid * xnew2, (1 - tmid) * x3 + tmid * xnew3);
        if (std::fabs(vmid) < 1e-2 * convergence_tol) break;
        if ((vlo < 0 && vmid > 0) ||
            (vlo > 0 && vmid < 0)) {  // if sign change between lo and mid
          thi = tmid;
          vhi = vmid;
        } else {  // otherwise sign change between hi and mid
          tlo = tmid;
          vlo = vmid;
        }
        if (iteration % 2)
          alpha =
              0.5;  // sometimes go with bisection to guarantee we make progress
        else
          alpha =
              vhi /
              (vhi -
               vlo);  // other times go with secant to hopefully get there fast
        tmid = alpha * tlo + (1 - alpha) * thi;
      }
      times[num_times++] = tmid;
    }
  }
  sort(times, num_times);

  if (errors) {
    for (int i = 0; i < num_times; ++i) {
      double ti = times[i];
      double val = std::fabs(signed_volume(
          (1 - ti) * x0 + ti * xnew0, (1 - ti) * x1 + ti * xnew1,
          (1 - ti) * x2 + ti * xnew2, (1 - ti) * x3 + ti * xnew3));
      errors[i] = val;
    }
  }
}

void getIntersectionPoint(const Vec3x& x_edge_0, const Vec3x& x_edge_1,
                          const Vec3x& x_face_0, const Vec3x& x_face_1,
                          const Vec3x& x_face_2, double* times, double* errors,
                          unsigned& num_times) {
  const double tol = 1e-12;
  num_times = 0;

  Vec3x x03 = x_edge_0 - x_face_2;
  Vec3x x13 = x_face_0 - x_face_2;
  Vec3x x23 = x_face_1 - x_face_2;
  Vec3x v03 = x_edge_1 - x_face_2 - x03;

  double C = triple(v03, x13, x23);
  double D = triple(x03, x13, x23);

  const double convergence_tol =
      tol * (std::fabs(0.0) + std::fabs(0.0) + std::fabs(C) + std::fabs(D));

  // find intervals to check, or just solve it if it reduces to a quadratic
  // =============================
  const int interval_times_size = 2;
  double interval_times[interval_times_size] = {0., 1.};

  // look for roots in indicated intervals
  // ============================================================== evaluate
  // coplanarity more accurately at each endpoint of the intervals
  double interval_values[interval_times_size];

  for (int i = 0; i < interval_times_size; ++i) {
    double t = interval_times[i];
    interval_values[i] = signed_volume((1 - t) * x_edge_0 + t * x_edge_1,
                                       x_face_0, x_face_1, x_face_2);
  }
  // first look for interval endpoints that are close enough to zero, without a
  // sign change
  for (int i = 0; i < interval_times_size; ++i) {
    if (interval_values[i] == 0) {
      times[num_times++] = interval_times[i];
    } else if (std::fabs(interval_values[i]) < convergence_tol) {
      if ((i == 0 || (interval_values[i - 1] >= 0 && interval_values[i] >= 0) ||
           (interval_values[i - 1] <= 0 && interval_values[i] <= 0)) &&
          (i == interval_times_size - 1 ||
           (interval_values[i + 1] >= 0 && interval_values[i] >= 0) ||
           (interval_values[i + 1] <= 0 && interval_values[i] <= 0))) {
        times[num_times++] = interval_times[i];
      }
    }
  }
  // and then search in intervals with a sign change
  for (int i = 1; i < interval_times_size; ++i) {
    double tlo = interval_times[i - 1], thi = interval_times[i], tmid;
    double vlo = interval_values[i - 1], vhi = interval_values[i], vmid;
    if ((vlo < 0 && vhi > 0) || (vlo > 0 && vhi < 0)) {
      // start off with secant approximation (in case the cubic is actually
      // linear)
      double alpha = vhi / (vhi - vlo);
      tmid = alpha * tlo + (1 - alpha) * thi;
      for (int iteration = 0; iteration < 50; ++iteration) {
        vmid = signed_volume((1 - tmid) * x_edge_0 + tmid * x_edge_1,
                             (1 - tmid) * x_face_0 + tmid * x_face_0,
                             (1 - tmid) * x_face_1 + tmid * x_face_1,
                             (1 - tmid) * x_face_2 + tmid * x_face_2);
        if (std::fabs(vmid) < 1e-2 * convergence_tol) break;
        if ((vlo < 0 && vmid > 0) ||
            (vlo > 0 && vmid < 0)) {  // if sign change between lo and mid
          thi = tmid;
          vhi = vmid;
        } else {  // otherwise sign change between hi and mid
          tlo = tmid;
          vlo = vmid;
        }
        if (iteration % 2)
          alpha =
              0.5;  // sometimes go with bisection to guarantee we make progress
        else
          alpha =
              vhi /
              (vhi -
               vlo);  // other times go with secant to hopefully get there fast
        tmid = alpha * tlo + (1 - alpha) * thi;
      }
      times[num_times++] = tmid;
    }
  }
  sort(times, num_times);

  if (errors) {
    for (int i = 0; i < num_times; ++i) {
      double ti = times[i];
      double val = std::fabs(signed_volume((1 - ti) * x_edge_0 + ti * x_edge_1,
                                           x_face_0, x_face_1, x_face_2));
      errors[i] = val;
    }
  }
}

void buildFrame(const Vec3x& n_hat, Vec3x& t1, Vec3x& t2) {
  // build collision frame basis (n_hat, t1, t2)
  Vec3x tmp(n_hat(0), -n_hat(2), n_hat(1));
  if (n_hat(2) * n_hat(2) < 1e-12 && n_hat(1) * n_hat(1) < 1e-12)
    tmp << -n_hat(1), n_hat(0), n_hat(2);
  t1 = tmp.cross(n_hat);

  if (t1.norm() < 1e-8) {
    ErrorStream(g_log, "") << "tangent degenerate: tmp =" << tmp
                           << ", n = " << n_hat << ", t1 =" << t1;
    return;
  }

  t1.normalize();
  t2 = t1.cross(n_hat);
  //    t2 = t2 / t2.norm();
  assert(isApproxUnit(t2));
}

Scalar getEllipticExternalCollisionOffset(const ElasticStrand& strand,
                                          const unsigned edge,
                                          const Vec3x& normal) {
  const CollisionParameters& params = strand.collisionParameters();
  return params.externalCollisionsRadius(edge);
}

bool analyseRoughRodRodCollision(const ProximityCollisionDatabase& database,
                                 const ElasticStrand* sP,
                                 const ElasticStrand* sQ, const int iP,
                                 const int iQ, const Scalar& contactAngle,
                                 Vec3x& depl, Scalar& s, Scalar& t, Scalar& d,
                                 Scalar& adhesion, Scalar& tilde_yield,
                                 Scalar& tilde_eta, Scalar& tilde_power,
                                 Scalar& relative_vel, bool& do_soc_solve) {
  const CollisionParameters& cpP = sP->collisionParameters();
  const CollisionParameters& cpQ = sQ->collisionParameters();

  const Vec3x& P0 = sP->getVertex(iP);
  const Vec3x& P1 = sP->getVertex(iP + 1);
  const Vec3x& Q0 = sQ->getVertex(iQ);
  const Vec3x& Q1 = sQ->getVertex(iQ + 1);

  //                std::cout << "Possible Coll: " << sP->m_globalIndex << " / "
  //                << iP
  //                                     << " vs " << sQ->m_globalIndex << " / "
  //                                     << iQ << std::endl ;

  // check flow info

  // use # collisions from last step to determine how to share the flow volume
  const int iSP = sP->getGlobalIndex();
  const int iSQ = sQ->getGlobalIndex();

  int nP = std::max(1, database.numCollisions(iSP, iP));
  int nQ = std::max(1, database.numCollisions(iSQ, iQ));

  Scalar aP0 = sP->getCurrentFlowDOFArea(iP);
  Scalar aP1 = sP->getCurrentFlowDOFArea(iP + 1);

  Scalar aQ0 = sQ->getCurrentFlowDOFArea(iQ);
  Scalar aQ1 = sQ->getCurrentFlowDOFArea(iQ + 1);

  const Scalar rP = cpP.selfCollisionsRadius(iP);
  const Scalar rQ = cpQ.selfCollisionsRadius(iQ);
  const Scalar BCRad = rP + rQ;

  Scalar sqDist =
      SquareDistSegmentToSegment<Vec3x, Scalar, Vec3x>(P0, P1, Q0, Q1, s, t);
  // Required to determnisticity -- are x87 registers sometimes used ?
  s = (float)s;
  t = (float)t;

  if (s < 0 || t < 0) return false;  // see FIXME in DistSegmentToSegment

  Scalar max_dist;

  const Vec3x& u0 = sP->getStepper()->velocities().segment<3>(iP * 4);
  const Vec3x& u1 = sP->getStepper()->velocities().segment<3>((iP + 1) * 4);
  const Vec3x& v0 = sQ->getStepper()->velocities().segment<3>(iQ * 4);
  const Vec3x& v1 = sQ->getStepper()->velocities().segment<3>((iQ + 1) * 4);

  const Vec3x us = (1. - s) * u0 + s * u1;
  const Vec3x vt = (1. - t) * v0 + t * v1;

  Vec3x PC = ((1. - s) * P0 + s * P1);
  Vec3x QC = ((1. - t) * Q0 + t * Q1);
  depl = (PC - QC);

  const Scalar n2depl = depl.squaredNorm();
  if (isSmall(n2depl)) return false;

  depl /= std::sqrt(n2depl);

  if (depl.dot((P1 - P0).normalized()) > COS_PARALLEL_ENOUGH ||
      depl.dot((Q1 - Q0).normalized()) < -COS_PARALLEL_ENOUGH)
    return false;

  relative_vel = (us - vt).dot(depl);
  const Scalar tangential_vel =
      sqrt(std::max(0., (us - vt).squaredNorm() - relative_vel * relative_vel));

  Scalar aP = (1. - s) * aP0 + s * aP1;
  Scalar aQ = (1. - t) * aQ0 + t * aQ1;

  Scalar aP_each = aP / (Scalar)nP;
  Scalar aQ_each = aQ / (Scalar)nQ;

  Scalar aPQ = aP_each + aQ_each;

  // separating, use wet dist to determine whether we'd enforce the constraint
  Scalar bc0 = cyl_h_from_area(rP, rP, aP_each) +
               cyl_h_from_area(rQ, rQ, aQ_each) + BCRad;
  Scalar bc1 = (1. + 0.5 * contactAngle) * sqrt(aPQ) + BCRad;

  if (database.connectedLastStep(iSP, iP, iSQ, iQ)) {
    // use [Lian et al.]'s criterion
    max_dist = std::max(bc0, bc1);
  } else {
    // check height only
    max_dist = bc0;
  }

  // volume on bridge
  if (sqDist > max_dist * max_dist) return false;

  if (aPQ < 1e-12 || relative_vel < -1e-4) {
    // dry or approaching, use solid to determine whether we should enforce the
    // constraint
    do_soc_solve = sqDist <= BCRad * BCRad;
  } else {
    // separating, do SOC solve anyways
    do_soc_solve = true;
  }

  // compute adhesion
  Scalar sqDist0 = std::min((P0 - Q0).squaredNorm(), (P0 - Q1).squaredNorm());
  Scalar sqDist1 = std::min((P1 - Q0).squaredNorm(), (P1 - Q1).squaredNorm());

  Scalar l0, l1;

  if (sqDist0 < max_dist * max_dist) {
    l0 = 1.0;
  } else {
    l0 = clamp((max_dist * max_dist - sqDist) / (sqDist0 - sqDist), 0.0, 1.0);
  }

  if (sqDist1 < max_dist * max_dist) {
    l1 = 1.0;
  } else {
    l1 = clamp((max_dist * max_dist - sqDist) / (sqDist1 - sqDist), 0.0, 1.0);
  }

  d = std::sqrt(sqDist);

  if (aPQ < 1e-12) {
    adhesion = 0.;
    tilde_yield = 0.;
    tilde_eta = 0.;
    tilde_power = 1.0;
  } else {
    if (&cpP == &cpQ) {
      const Scalar r =
          sqrt(((sP->getRadiusA(iP) * sP->getRadiusB(iP) * (1. - s) +
                 sP->getRadiusA(iP + 1) * sP->getRadiusB(iP + 1) * s) +
                (sP->getRadiusA(iQ) * sP->getRadiusB(iQ) * (1. - t) +
                 sP->getRadiusA(iQ + 1) * sP->getRadiusB(iQ + 1) * t)) *
               0.5);

      VecXx color = sP->getFlowNewComponents(iP) + sQ->getFlowNewComponents(iQ);
      make_gibbs_simplex(color);

      const Scalar eta = sP->getFlowConsistencyIndex(color);
      const Scalar n = sP->getFlowBehaviorIndex(color);
      const Scalar tilde_sigma_Y = sP->getFlowYieldStress(color) * 0.8164965809;

      const Scalar diff_vel = pow(fabs(relative_vel) / std::max(r * 2.0, d), n);
      const Scalar diff_vel_horizontal =
          pow(tangential_vel / std::max(r * 2.0, d), n);

      const Scalar Ac = 0.75 * M_PI * r;
      const Scalar elastic_force = Ac * (eta * diff_vel + tilde_sigma_Y);

      Scalar fP = std::max(0., cpP.adhesionForce(aPQ, d, color) +
                                   elastic_force);  // unit: dyn / cm
      Scalar coeff = ((s * l0 + (1. - s) * l1) * (P0 - P1).norm() +
                      (t * l0 + (1. - t) * l1) * (Q0 - Q1).norm()) *
                     0.5;

      adhesion = coeff * fP;

      tilde_yield = coeff * Ac * tilde_sigma_Y;
      tilde_eta = coeff * Ac * (eta / pow(std::max(r * 2.0, d), n));
      tilde_power = n;
    } else {
      const Scalar r =
          sqrt(((sP->getRadiusA(iP) * sP->getRadiusB(iP) * (1. - s) +
                 sP->getRadiusA(iP + 1) * sP->getRadiusB(iP + 1) * s) +
                (sP->getRadiusA(iQ) * sP->getRadiusB(iQ) * (1. - t) +
                 sP->getRadiusA(iQ + 1) * sP->getRadiusB(iQ + 1) * t)) *
               0.5);

      VecXx color = sP->getFlowNewComponents(iP) + sQ->getFlowNewComponents(iQ);
      make_gibbs_simplex(color);

      const Scalar eta = sP->getFlowConsistencyIndex(color);
      const Scalar n = sP->getFlowBehaviorIndex(color);
      const Scalar tilde_sigma_Y = sP->getFlowYieldStress(color) * 0.8164965809;

      const Scalar diff_vel = pow(fabs(relative_vel) / std::max(r * 2.0, d), n);
      const Scalar diff_vel_horizontal =
          pow(tangential_vel / std::max(r * 2.0, d), n);

      const Scalar Ac = 0.75 * M_PI * r;
      const Scalar elastic_force = Ac * (eta * diff_vel + tilde_sigma_Y);

      Scalar fP = std::max(0., cpP.adhesionForce(aPQ, d, color) +
                                   elastic_force);  // unit: dyn / cm
      Scalar fQ =
          std::max(0., cpQ.adhesionForce(aPQ, d, color) + elastic_force);

      adhesion = ((s * l0 + (1. - s) * l1) * (P0 - P1).norm() * fP +
                  (t * l0 + (1. - t) * l1) * (Q0 - Q1).norm() * fQ) *
                 0.5;  // unit: dyn

      Scalar coeff = ((s * l0 + (1. - s) * l1) * (P0 - P1).norm() +
                      (t * l0 + (1. - t) * l1) * (Q0 - Q1).norm()) *
                     0.5;

      tilde_yield = coeff * Ac * tilde_sigma_Y;
      tilde_eta = coeff * Ac * (eta / pow(std::max(r * 2.0, d), n));
      tilde_power = n;
    }
  }

  return true;
}

}  // namespace strandsim
