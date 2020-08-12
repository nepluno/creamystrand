/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "DegreesOfFreedom.hh"

#include <mkl.h>

#include "../Core/ElasticStrandParameters.hh"
#include "../Core/ElasticStrandUtils.hh"
#include "../Utils/MathUtilities.hh"
//#else
// extern "C" void vdSinCos(const int n, const double a[], double r1[], double
// r2[]); #endif

//#define CURVATURE_USE_SINTHETA

#ifndef CURVATURE_USE_SINTHETA
#define CURVATURE_USE_TANTHETA
#endif

namespace strandsim {
Height::Height(AreaDOFs& dofs, const ElasticStrandParameters& parameters)
    : DependencyNode<std::vector<Scalar> >(0, dofs.get().size()),
      m_area_dofs(dofs),
      m_parameters(parameters) {
  m_area_dofs.addDependent(this);
}

void Height::compute() {
  m_value.resize(m_size);
  const VecXx& dofs = m_area_dofs.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    Scalar ra = m_parameters.getRadiusA(vtx, m_size);
    Scalar rb = m_parameters.getRadiusB(vtx, m_size);
    m_value[vtx] = cyl_h_from_area(ra, rb, dofs[vtx]);
  }

  // check_isnan("height_compute", m_value);

  setDependentsDirty();
}

void Edges::compute() {
  m_value.resize(m_size);
  const VecXx& dofs = m_dofs.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    assert(vtx >= 0 && vtx < m_size);
    m_value[vtx] = dofs.segment<3>(4 * (vtx + 1)) - dofs.segment<3>(4 * vtx);
    assert(m_value[vtx][0] < 1e+2);
    assert(m_value[vtx][1] < 1e+2);
    assert(m_value[vtx][2] < 1e+2);
  }

  setDependentsDirty();
}

void JumpEdges::compute() {
  m_value.resize(m_size);
  const VecXx& dofs = m_dofs.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    if (vtx == m_firstValidIndex || vtx == size() - 1)
      m_value[vtx].setZero();
    else {
      m_value[vtx] =
          dofs.segment<3>(4 * (vtx + 1)) - dofs.segment<3>(4 * (vtx - 1));
    }
  }

  setDependentsDirty();
}

void Lengths::compute() {
  m_value.resize(m_size);
  const Vec3xArray& edges = m_edges.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    assert(vtx >= 0 && vtx < m_size);
    m_value[vtx] = edges[vtx].norm();
    // assert( !isSmall(m_value[vtx]) ); // Commented-out assert, as it be may
    // thrown while we're checking stuff
    assert(m_value[vtx] < 1e+2);
  }

  // check_isnan("length_compute", m_value);

  setDependentsDirty();
}

void NeighborDistances::compute() {
  m_value.resize(m_size);
  const Vec3xArray& edges = m_jump_edges.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    if (vtx == m_firstValidIndex || vtx == size() - 1)
      m_value[vtx] = 0.;
    else {
      m_value[vtx] = edges[vtx].norm();
    }
  }

  setDependentsDirty();
}

void DynamicVoronoiLengths::compute() {
  m_value.resize(m_size);

  m_value[0] = m_lengths[0] * 0.5;

  const int num_edges = m_lengths.size();

  for (IndexType vtx = 1; vtx < num_edges; ++vtx) {
    assert(vtx >= 0 && vtx < m_size);
    m_value[vtx] = (m_lengths[vtx - 1] + m_lengths[vtx]) * 0.5;
    // assert( !isSmall(m_value[vtx]) ); // Commented-out assert, as it be may
    // thrown while we're checking stuff
  }
  m_value[num_edges] = 0.5 * m_lengths[num_edges - 1];

  // check_isnan("dvl_compute", m_value);

  setDependentsDirty();
}

void LaplaceHeight::compute() {
  m_value.resize(m_size);

  const IndexType num_vtx = size();
  for (IndexType vtx = m_firstValidIndex; vtx < num_vtx; ++vtx) {
    Scalar h0, h1, h2;
    if (vtx == m_firstValidIndex)
      h0 = 0.0;
    else
      h0 = m_heights[vtx - 1];
    if (vtx == num_vtx - 1)
      h2 = 0.0;
    else
      h2 = m_heights[vtx + 1];
    h1 = m_heights[vtx];

    const Scalar l = m_voronoi_lengths[vtx];

    m_value[vtx] = (h0 + h2 - h1 * 2.0) / (l * l);
  }

  // check_isnan("lh_compute", m_value);

  setDependentsDirty();
}

void Tangents::compute() {
  m_value.resize(m_size);
  const Vec3xArray& edges = m_edges.get();
  const std::vector<Scalar>& lengths = m_lengths.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    assert(vtx >= 0 && vtx < m_size);
    m_value[vtx] = edges[vtx] / lengths[vtx];
    assert(m_value[vtx][0] < 1e+2);
    assert(m_value[vtx][1] < 1e+2);
    assert(m_value[vtx][2] < 1e+2);
  }

  setDependentsDirty();
}

void CurvatureBinormals::compute() {
  m_value.resize(m_size);
  const Vec3xArray& tangents = m_tangents.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    const Vec3x& t1 = tangents[vtx - 1];
    const Vec3x& t2 = tangents[vtx];

#ifdef CURVATURE_USE_SINTHETA
    Scalar denominator = (t1 + t2).norm();
#else
    Scalar denominator = 1.0 + t1.dot(t2);
#endif

    if (denominator < 1e-12) {
      std::ostringstream oss;
      oss << "CurvatureBinormals::compute() denominator == " << denominator
          << " at vertex " << vtx;
      oss << " t1 = " << t1 << " t2 = " << t2;
      WarningStream(g_log, "") << oss.str();

      denominator = 1e-12;
    }

    m_value[vtx] = 2.0 * t1.cross(t2) / denominator;
  }

  setDependentsDirty();
}

void TrigThetas::compute() {
  const VecXx& dofs = m_dofs.get();
  const IndexType numThetas = m_dofs.getNumEdges();
  m_value.first.resize(numThetas);
  m_value.second.resize(numThetas);

  // Extract thetas in their own vector for mkl_vlm
  const Eigen::Map<const VecXx, Eigen::Unaligned, Eigen::InnerStride<4> >
      thetasMap(dofs.data() + 3, numThetas);
  const VecXx thetaVec(thetasMap);
  // Compute their sine and cosine
  assert(typeid(double) == typeid(VecXx::Scalar));
  vdSinCos(
      numThetas, thetaVec.data(), m_value.first.data(),
      m_value.second.data());  // FIXME this won't compile if Scalar != double

  setDependentsDirty();
}

}  // namespace strandsim
