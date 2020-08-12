/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "StrandState.hh"

#include "../Utils/Distances.hh"
#include "../Utils/TextLog.hh"
#include "BandMatrix.hh"
#include "ElasticStrandParameters.hh"
#include "ElasticStrandUtils.hh"

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Sparse>

namespace strandsim {
StrandState::StrandState(const VecXx &dofs, const VecXx &area_dofs,
                         BendingMatrixBase &bendingMatrixBase,
                         const ElasticStrandParameters &parameters)
    : SimulatedObjectState(),
      m_numVertices((dofs.size() + 1) / 4),  //
      m_area_dofs(area_dofs),
      m_dofs(dofs),     //
      m_edges(m_dofs),  //
      m_jump_edges(m_dofs),
      m_neighbor_distances(m_jump_edges),
      m_lengths(m_edges),  //
      m_dynamic_voronoi_lengths(m_lengths),
      m_height(m_area_dofs, parameters),
      m_laplace_height(m_height, m_dynamic_voronoi_lengths),
      m_tangents(m_edges, m_lengths),                      //
      m_referenceFrames1(m_tangents),                      //
      m_referenceFrames2(m_tangents, m_referenceFrames1),  //
      m_referenceTwists(m_tangents, m_referenceFrames1),   //
      m_twists(m_referenceTwists, m_dofs),                 //
      m_curvatureBinormals(m_tangents),                    //
      m_trigThetas(m_dofs),                                //
      m_materialFrames1(m_trigThetas, m_referenceFrames1,
                        m_referenceFrames2),  //
      m_materialFrames2(m_trigThetas, m_referenceFrames1,
                        m_referenceFrames2),                                 //
      m_kappas(m_curvatureBinormals, m_materialFrames1, m_materialFrames2),  //
      m_gradKappas(m_lengths, m_tangents, m_curvatureBinormals,
                   m_materialFrames1, m_materialFrames2, m_kappas),  //
      m_gradTwists(m_lengths, m_curvatureBinormals),                 //
      m_gradTwistsSquared(m_gradTwists),                             //
      m_hessKappas(m_lengths, m_tangents, m_curvatureBinormals,
                   m_materialFrames1, m_materialFrames2, m_kappas),  //
      m_hessTwists(m_tangents, m_lengths, m_curvatureBinormals),     //
      m_bendingProducts(bendingMatrixBase, m_gradKappas) {}

StrandState::StrandState(const VecXx &dofs,
                         BendingMatrixBase &bendingMatrixBase,
                         const ElasticStrandParameters &parameters)
    : SimulatedObjectState(),
      m_numVertices((dofs.size() + 1) / 4),  //
      m_area_dofs(VecXx::Zero((dofs.size() + 1) / 4)),
      m_dofs(dofs),     //
      m_edges(m_dofs),  //
      m_jump_edges(m_dofs),
      m_neighbor_distances(m_jump_edges),
      m_lengths(m_edges),  //
      m_dynamic_voronoi_lengths(m_lengths),
      m_height(m_area_dofs, parameters),
      m_laplace_height(m_height, m_dynamic_voronoi_lengths),
      m_tangents(m_edges, m_lengths),                      //
      m_referenceFrames1(m_tangents),                      //
      m_referenceFrames2(m_tangents, m_referenceFrames1),  //
      m_referenceTwists(m_tangents, m_referenceFrames1),   //
      m_twists(m_referenceTwists, m_dofs),                 //
      m_curvatureBinormals(m_tangents),                    //
      m_trigThetas(m_dofs),                                //
      m_materialFrames1(m_trigThetas, m_referenceFrames1,
                        m_referenceFrames2),  //
      m_materialFrames2(m_trigThetas, m_referenceFrames1,
                        m_referenceFrames2),                                 //
      m_kappas(m_curvatureBinormals, m_materialFrames1, m_materialFrames2),  //
      m_gradKappas(m_lengths, m_tangents, m_curvatureBinormals,
                   m_materialFrames1, m_materialFrames2, m_kappas),  //
      m_gradTwists(m_lengths, m_curvatureBinormals),                 //
      m_gradTwistsSquared(m_gradTwists),                             //
      m_hessKappas(m_lengths, m_tangents, m_curvatureBinormals,
                   m_materialFrames1, m_materialFrames2, m_kappas),  //
      m_hessTwists(m_tangents, m_lengths, m_curvatureBinormals),     //
      m_bendingProducts(bendingMatrixBase, m_gradKappas) {}

StrandState::StrandState(const size_t dofSize,
                         BendingMatrixBase &bendingMatrixBase,
                         const ElasticStrandParameters &parameters)
    : m_numVertices((dofSize + 1) / 4),  //
      m_area_dofs(VecXx::Zero((dofSize + 1) / 4)),
      m_dofs(VecXx(dofSize)),  //
      m_edges(m_dofs),         //
      m_jump_edges(m_dofs),
      m_neighbor_distances(m_jump_edges),
      m_lengths(m_edges),  //
      m_dynamic_voronoi_lengths(m_lengths),
      m_height(m_area_dofs, parameters),
      m_laplace_height(m_height, m_dynamic_voronoi_lengths),
      m_tangents(m_edges, m_lengths),                      //
      m_referenceFrames1(m_tangents),                      //
      m_referenceFrames2(m_tangents, m_referenceFrames1),  //
      m_referenceTwists(m_tangents, m_referenceFrames1),   //
      m_twists(m_referenceTwists, m_dofs),                 //
      m_curvatureBinormals(m_tangents),                    //
      m_trigThetas(m_dofs),                                //
      m_materialFrames1(m_trigThetas, m_referenceFrames1,
                        m_referenceFrames2),  //
      m_materialFrames2(m_trigThetas, m_referenceFrames1,
                        m_referenceFrames2),                                 //
      m_kappas(m_curvatureBinormals, m_materialFrames1, m_materialFrames2),  //
      m_gradKappas(m_lengths, m_tangents, m_curvatureBinormals,
                   m_materialFrames1, m_materialFrames2, m_kappas),  //
      m_gradTwists(m_lengths, m_curvatureBinormals),                 //
      m_gradTwistsSquared(m_gradTwists),                             //
      m_hessKappas(m_lengths, m_tangents, m_curvatureBinormals,
                   m_materialFrames1, m_materialFrames2, m_kappas),  //
      m_hessTwists(m_tangents, m_lengths, m_curvatureBinormals),     //
      m_bendingProducts(bendingMatrixBase, m_gradKappas) {}

StrandState::~StrandState() {}

void StrandState::resizeSelf() {
  const size_t ndofs = getDegreesOfFreedom().size();
  assert(ndofs % 4 == 3);
  // dofs are 3 per vertex, one per edge
  assert(ndofs > 3);
  // minimum two vertices per rod

  m_totalForce.resize(ndofs);
  m_totalJacobian->resize(static_cast<IndexType>(ndofs),
                          static_cast<IndexType>(ndofs));
}

void StrandState::freeCachedQuantities() {
  m_curvatureBinormals.free();
  m_trigThetas.free();
  m_gradKappas.free();
  m_gradTwists.free();
  m_gradTwistsSquared.free();
  m_hessKappas.free();
  m_hessTwists.free();
  m_bendingProducts.free();
}

bool StrandState::hasSmallForces(const Scalar lTwoTol,
                                 const Scalar lInfTol) const {
  return ((m_totalForce.norm() / m_numVertices <= lTwoTol) ||
          (m_totalForce.lpNorm<Eigen::Infinity>() <= lInfTol));
}

void StrandState::computeDeformationGradient(const unsigned edge,
                                             const Scalar alpha,
                                             const VecXx &velocities,
                                             SparseRowMatx *&pH) const {
  if (pH) {
    pH->resize(3, getDegreesOfFreedom().rows());
  } else {
    pH = new SparseRowMatx(3, getDegreesOfFreedom().rows());
  }

  SparseRowMatx &H = *pH;
  H.reserve(alpha > 0. ? 6 : 3);

  for (unsigned k = 0; k < 3; ++k) {
    H.startVec(k);
    const unsigned col = 4 * edge + k;
    H.insertBackByOuterInner(k, col) = (1. - alpha);

    if (alpha > 0.) {
      H.insertBackByOuterInner(k, col + 4) = alpha;
    }
  }
  H.finalize();
}

Vec3x StrandState::closestPoint(const Vec3x &x) const {
  Scalar mindist = std::numeric_limits<Scalar>::max();
  Vec3x winner;

  for (int vtx = 0; vtx < m_numVertices - 1; ++vtx) {
    Vec3x y = ClosestPtPointSegment(x, getVertex(vtx), getVertex(vtx + 1));
    Scalar dist = (y - x).squaredNorm();
    if (dist < mindist) {
      mindist = dist;
      winner = y;
    }
  }

  return winner;
}

std::pair<Vec3x, Vec3x> StrandState::getAABB() const {
  std::pair<Vec3x, Vec3x> aaBB;
  Vec3x &min = aaBB.first;
  Vec3x &max = aaBB.second;

  min = getVertex(0);
  max = getVertex(0);

  for (unsigned i = 1; i < m_numVertices; ++i) {
    const Vec3x &v = getVertex(i);
    min = min.cwiseMin(v);
    max = max.cwiseMax(v);
  }

  return aaBB;
}

std::pair<Vec3x, Vec3x> StrandState::getAABB(unsigned elementID) const {
  std::pair<Vec3x, Vec3x> aaBB;
  Vec3x &min = aaBB.first;
  Vec3x &max = aaBB.second;

  const Vec3x start = getVertex(elementID);
  if (1 + elementID == m_numVertices) {
    min = max = start;
  } else {
    const Vec3x end = getVertex(elementID + 1);

    min = start.cwiseMin(end);
    max = start.cwiseMax(end);
  }

  return aaBB;
}

std::pair<Vec3x, Vec3x> StrandState::getAABB(const Mat4x &transform) const {
  std::pair<Vec3x, Vec3x> aaBB;
  Vec3x &min = aaBB.first;
  Vec3x &max = aaBB.second;

  min = transformPoint(transform, getVertex(0));
  max = min;

  for (unsigned i = 1; i < m_numVertices; ++i) {
    const Vec3x &vi = transformPoint(transform, getVertex(i));

    min = min.cwiseMin(vi);
    max = max.cwiseMax(vi);
  }

  return aaBB;
}

std::pair<Vec3x, Vec3x> StrandState::getAABB(unsigned elementID,
                                             const Mat4x &transform) const {
  std::pair<Vec3x, Vec3x> aaBB;
  Vec3x &min = aaBB.first;
  Vec3x &max = aaBB.second;

  const Vec3x start = transformPoint(transform, getVertex(elementID));

  if (1 + elementID == m_numVertices) {
    min = max = start;
  } else {
    const Vec3x end = transformPoint(transform, getVertex(elementID + 1));

    min = start.cwiseMin(end);
    max = start.cwiseMax(end);
  }

  return aaBB;
}

}  // namespace strandsim
