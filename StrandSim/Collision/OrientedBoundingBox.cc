/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "OrientedBoundingBox.hh"

#include <cmath>

#include "../Core/ElasticStrand.hh"

namespace strandsim {

OrientedBoundingBox::OrientedBoundingBox(
    const ElasticStrand& strand, const unsigned edge,
    CollisionParameters::CollisionType type) {
  init(strand.getCurrentState(), edge, strand.collisionParameters(), type);
}

OrientedBoundingBox::OrientedBoundingBox(
    const StrandState& geometry, const unsigned edge,
    const CollisionParameters& cp, CollisionParameters::CollisionType type) {
  init(geometry, edge, cp, type);
}

void OrientedBoundingBox::init(const StrandState& geometry, const unsigned edge,
                               const CollisionParameters& cp,
                               CollisionParameters::CollisionType type) {
  const Vec3x& P0 = geometry.getVertex(edge);
  const Vec3x& P1 = geometry.getVertex(edge + 1);

  m_center = .5 * (P0 + P1);

  m_axis[0] = (P1 - m_center);
  m_extents[0] = m_axis[0].norm();
  m_axis[0] /= m_extents[0];

  m_axis[1] = geometry.getMaterialFrame1(edge);
  m_axis[2] = geometry.getMaterialFrame2(edge);

  m_extents[1] = cp.collisionsRadius(type, edge, M_PI / 2);
  m_extents[2] = cp.collisionsRadius(type, edge, 0);
}

void OrientedBoundingBox::getAABB(Vec3x& min, Vec3x& max) {
  Vec3x origin = m_center;
  for (unsigned k = 0; k < 3; ++k) {
    origin -= m_axis[k] * m_extents[k];
  }
  min = origin;
  max = origin;

  Vec3x edges[3];
  for (unsigned k = 0; k < 3; ++k) {
    edges[k] = 2. * m_axis[k] * m_extents[k];
  }

  Vec3x vertex;
  for (int i = 1; i < 8; ++i) {
    vertex = origin;
    for (unsigned k = 0; k < 3; ++k) {
      vertex += ((i & (1 << k)) >> k) * edges[k];
    }

    min = min.cwiseMin(vertex);
    max = max.cwiseMax(vertex);
  }
}

// Adapated from
// http://www.geometrictools.com/LibMathematics/Intersection/Wm5IntrBox3Box3.cpp
// Geometric Tools, LLC
// Copyright (c) 1998-2012
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

bool OrientedBoundingBox::areIntersecting(const OrientedBoundingBox& obb0,
                                          const OrientedBoundingBox& obb1) {
  typedef Scalar Real;

  // Cutoff for cosine of angles between box axes.  This is used to catch
  // the cases when at least one pair of axes are parallel.  If this
  // happens, there is no need to test for separation along the
  // Cross(A[i],B[j]) directions.
  const Real cutoff =
      (Real)1 - SMALL_NUMBER<float>();  // Math<Real>::ZERO_TOLERANCE;
  bool existsParallelPair = false;
  int i;

  // Convenience variables.
  const Vec3x* A = obb0.m_axis;
  const Vec3x* B = obb1.m_axis;
  const Real* EA = obb0.m_extents;
  const Real* EB = obb1.m_extents;

  // Compute difference of box centers, D = C1-C0.
  const Vec3x D = obb1.m_center - obb0.m_center;

  Real C[3][3];     // matrix C = A^T B, c_{ij} = Dot(A_i,B_j)
  Real AbsC[3][3];  // |c_{ij}|
  Real AD[3];       // Dot(A_i,D)
  Real r0, r1, r;   // interval radii and distance between centers
  Real r01;         // = R0 + R1

  // axis C0+t*A0
  for (i = 0; i < 3; ++i) {
    C[0][i] = A[0].dot(B[i]);
    AbsC[0][i] = std::fabs(C[0][i]);
    if (AbsC[0][i] > cutoff) {
      existsParallelPair = true;
    }
  }
  AD[0] = A[0].dot(D);
  r = std::fabs(AD[0]);
  r1 = EB[0] * AbsC[0][0] + EB[1] * AbsC[0][1] + EB[2] * AbsC[0][2];
  r01 = EA[0] + r1;
  if (r > r01) {
    return false;
  }

  // axis C0+t*A1
  for (i = 0; i < 3; ++i) {
    C[1][i] = A[1].dot(B[i]);
    AbsC[1][i] = std::fabs(C[1][i]);
    if (AbsC[1][i] > cutoff) {
      existsParallelPair = true;
    }
  }
  AD[1] = A[1].dot(D);
  r = std::fabs(AD[1]);
  r1 = EB[0] * AbsC[1][0] + EB[1] * AbsC[1][1] + EB[2] * AbsC[1][2];
  r01 = EA[1] + r1;
  if (r > r01) {
    return false;
  }

  // axis C0+t*A2
  for (i = 0; i < 3; ++i) {
    C[2][i] = A[2].dot(B[i]);
    AbsC[2][i] = std::fabs(C[2][i]);
    if (AbsC[2][i] > cutoff) {
      existsParallelPair = true;
    }
  }
  AD[2] = A[2].dot(D);
  r = std::fabs(AD[2]);
  r1 = EB[0] * AbsC[2][0] + EB[1] * AbsC[2][1] + EB[2] * AbsC[2][2];
  r01 = EA[2] + r1;
  if (r > r01) {
    return false;
  }

  // axis C0+t*B0
  r = std::fabs(B[0].dot(D));
  r0 = EA[0] * AbsC[0][0] + EA[1] * AbsC[1][0] + EA[2] * AbsC[2][0];
  r01 = r0 + EB[0];
  if (r > r01) {
    return false;
  }

  // axis C0+t*B1
  r = std::fabs(B[1].dot(D));
  r0 = EA[0] * AbsC[0][1] + EA[1] * AbsC[1][1] + EA[2] * AbsC[2][1];
  r01 = r0 + EB[1];
  if (r > r01) {
    return false;
  }

  // axis C0+t*B2
  r = std::fabs(B[2].dot(D));
  r0 = EA[0] * AbsC[0][2] + EA[1] * AbsC[1][2] + EA[2] * AbsC[2][2];
  r01 = r0 + EB[2];
  if (r > r01) {
    return false;
  }

  // At least one pair of box axes was parallel, so the separation is
  // effectively in 2D where checking the "edge" normals is sufficient for
  // the separation of the boxes.
  if (existsParallelPair) {
    return true;
  }

  // axis C0+t*A0xB0
  r = std::fabs(AD[2] * C[1][0] - AD[1] * C[2][0]);
  r0 = EA[1] * AbsC[2][0] + EA[2] * AbsC[1][0];
  r1 = EB[1] * AbsC[0][2] + EB[2] * AbsC[0][1];
  r01 = r0 + r1;
  if (r > r01) {
    return false;
  }

  // axis C0+t*A0xB1
  r = std::fabs(AD[2] * C[1][1] - AD[1] * C[2][1]);
  r0 = EA[1] * AbsC[2][1] + EA[2] * AbsC[1][1];
  r1 = EB[0] * AbsC[0][2] + EB[2] * AbsC[0][0];
  r01 = r0 + r1;
  if (r > r01) {
    return false;
  }

  // axis C0+t*A0xB2
  r = std::fabs(AD[2] * C[1][2] - AD[1] * C[2][2]);
  r0 = EA[1] * AbsC[2][2] + EA[2] * AbsC[1][2];
  r1 = EB[0] * AbsC[0][1] + EB[1] * AbsC[0][0];
  r01 = r0 + r1;
  if (r > r01) {
    return false;
  }

  // axis C0+t*A1xB0
  r = std::fabs(AD[0] * C[2][0] - AD[2] * C[0][0]);
  r0 = EA[0] * AbsC[2][0] + EA[2] * AbsC[0][0];
  r1 = EB[1] * AbsC[1][2] + EB[2] * AbsC[1][1];
  r01 = r0 + r1;
  if (r > r01) {
    return false;
  }

  // axis C0+t*A1xB1
  r = std::fabs(AD[0] * C[2][1] - AD[2] * C[0][1]);
  r0 = EA[0] * AbsC[2][1] + EA[2] * AbsC[0][1];
  r1 = EB[0] * AbsC[1][2] + EB[2] * AbsC[1][0];
  r01 = r0 + r1;
  if (r > r01) {
    return false;
  }

  // axis C0+t*A1xB2
  r = std::fabs(AD[0] * C[2][2] - AD[2] * C[0][2]);
  r0 = EA[0] * AbsC[2][2] + EA[2] * AbsC[0][2];
  r1 = EB[0] * AbsC[1][1] + EB[1] * AbsC[1][0];
  r01 = r0 + r1;
  if (r > r01) {
    return false;
  }

  // axis C0+t*A2xB0
  r = std::fabs(AD[1] * C[0][0] - AD[0] * C[1][0]);
  r0 = EA[0] * AbsC[1][0] + EA[1] * AbsC[0][0];
  r1 = EB[1] * AbsC[2][2] + EB[2] * AbsC[2][1];
  r01 = r0 + r1;
  if (r > r01) {
    return false;
  }

  // axis C0+t*A2xB1
  r = std::fabs(AD[1] * C[0][1] - AD[0] * C[1][1]);
  r0 = EA[0] * AbsC[1][1] + EA[1] * AbsC[0][1];
  r1 = EB[0] * AbsC[2][2] + EB[2] * AbsC[2][0];
  r01 = r0 + r1;
  if (r > r01) {
    return false;
  }

  // axis C0+t*A2xB2
  r = std::fabs(AD[1] * C[0][2] - AD[0] * C[1][2]);
  r0 = EA[0] * AbsC[1][2] + EA[1] * AbsC[0][2];
  r1 = EB[0] * AbsC[2][1] + EB[1] * AbsC[2][0];
  r01 = r0 + r1;
  if (r > r01) {
    return false;
  }

  return true;
}

} /* namespace strandsim */
