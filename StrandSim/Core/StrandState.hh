/**
 * \copyright 2011 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef STRANDSTATE_HH_
#define STRANDSTATE_HH_

#include <memory>

#include "../Dependencies/BendingProducts.hh"
#include "../Dependencies/Twists.hh"
#include "BandMatrixFwd.hh"
#include "Definitions.hh"

namespace strandsim {
class ElasticStrandParameters;

class SimulatedObjectState {
 public:
  SimulatedObjectState() : m_totalEnergy(0) {}
  virtual ~SimulatedObjectState() {}

  // Energy, force
  Scalar m_totalEnergy;
  VecXx m_totalForce;

  // The Jacobian is only needed once. To avoid memory duplication, share a
  // pointer between current and future geometries.
  std::shared_ptr<JacobianMatrixType> m_totalJacobian;
};

/**
 * \brief This class contains all the strand's variable data (i.e. changed by
 * simulation).
 */
class StrandState : public SimulatedObjectState {
  friend class ElasticStrand;

 public:
  explicit StrandState(const VecXx& dofs, const VecXx& height_dofs,
                       BendingMatrixBase& bendingMatrixBase,
                       const ElasticStrandParameters& parameters);

  explicit StrandState(const VecXx& dofs, BendingMatrixBase& bendingMatrixBase,
                       const ElasticStrandParameters& parameters);

  explicit StrandState(const size_t dofSize,
                       BendingMatrixBase& bendingMatrixBase,
                       const ElasticStrandParameters& parameters);

  ~StrandState();

  const VecXx& getDegreesOfFreedom() const { return m_dofs.get(); }
  void setDegreesOfFreedom(const VecXx& dof) { m_dofs.set(dof); }
  const VecXx& getAreaDegreesOfFreedom() const { return m_area_dofs.get(); }
  Scalar getAreaDegreesOfFreedom(int vtx) const {
    return m_area_dofs.get()(clamp(vtx, 0, m_numVertices - 1));
  }
  void setAreaDegreesOfFreedom(const VecXx& dof) { m_area_dofs.set(dof); }
  void setAreaDegreesOfFreedom(const IndexType vtx, Scalar value) {
    m_area_dofs.setVertex(vtx, value);
  }

  const Eigen::Map<const VecXx, Eigen::Unaligned, Eigen::InnerStride<4> >
  getThetas() const {
    return m_dofs.getThetas();
  }
  void setThetas(const VecXx& thetas, int numberOfFixedThetas = 0) {
    m_dofs.setThetas(thetas, numberOfFixedThetas);
  }

  const Vec3x getVertex(const IndexType vtx) const {
    assert(vtx < m_numVertices);

    return m_dofs.getVertex(vtx);
  }

  Scalar getHeight(const IndexType vtx) {
    int clamped_vtx = (int)clamp(vtx, 0, m_numVertices - 1);
    return m_height[clamped_vtx];
  }

  //! Gets a vertex somewhere on an edge
  const Vec3x getVertex(const IndexType vtx, Scalar localAbscissa) const {
    const Vec3x v0 = getVertex(vtx);
    if (localAbscissa > 0.) {
      return v0 + (localAbscissa) * (getVertex(vtx + 1) - v0);
    }
    return v0;
  }

  Vec3x getVertex(const IndexType vtx, const Mat4x& transform) const {
    assert(vtx < m_numVertices);

    Vec4x v;
    v.head<3>() = m_dofs.getVertex(vtx);
    v[3] = 1.;
    return Vec4x(transform * v).head<3>();
  }

  //! Gets a vertex somewhere on an edge
  const Vec3x getVertex(const IndexType vtx, Scalar localAbscissa,
                        const Mat4x& transform) const {
    assert(vtx < m_numVertices);

    Vec4x v;
    v.head<3>() = getVertex(vtx, localAbscissa);
    v[3] = 1.;
    return Vec4x(transform * v).head<3>();
  }

  Scalar getTheta(const IndexType vtx) const {
    assert(vtx < m_numVertices - 1);

    return m_dofs.getTheta(vtx);
  }

  void setTheta(const IndexType vtx, const Scalar newTheta) {
    assert(vtx < m_numVertices - 1);

    m_dofs.setTheta(vtx, newTheta);
  }

  Scalar getVertexLength(const IndexType vtx) {
    return m_dynamic_voronoi_lengths[vtx];
  }

  Vec3x getJumpEdge(const IndexType vtx) { return m_jump_edges.get()[vtx]; }

  Scalar getNeighborDistances(const IndexType vtx) {
    return m_neighbor_distances[vtx];
  }

  Scalar getEdgeLength(const IndexType vtx) { return m_lengths[vtx]; }

  const Vec3x getEdgeVector(const IndexType vtx) const {
    if (vtx < m_numVertices - 1)
      return getVertex(vtx + 1) - getVertex(vtx);
    else
      return getVertex(vtx) - getVertex(vtx - 1);
  }

  inline int numVertices() const { return m_numVertices; }

  void storeInitialFrames(const Vec3x& initRefFrame1 = Vec3x()) {
    m_referenceFrames1.storeInitialFrames(initRefFrame1);
  }

  void setVertex(const IndexType vtx, const Vec3x& coordinates) {
    assert(vtx < m_numVertices);

    m_dofs.setVertex(vtx, coordinates);
  }

  void resizeSelf();
  void freeCachedQuantities();

  //! Computes H and uFree such that v = H q + uFree
  /*! where v is the world velocitiy of the point at (edge, alpha)
   and q the non-fixed degrees of freedom .
   Allocates H if pH ( pointer to H ) is NULL
   */
  void computeDeformationGradient(const unsigned edge, const Scalar alpha,
                                  const VecXx& velocities,
                                  SparseRowMatx*& pH) const;

  const Vec3x getReferenceFrame1(const IndexType vtx) const {
    assert(vtx < m_numVertices - 1);

    return m_referenceFrames1[vtx];
  }

  void setReferenceFrame1(const IndexType vtx, const Vec3x& vec) {
    assert(vtx < m_numVertices - 1);

    m_referenceFrames1.set(vtx, vec);
  }

  const Vec3x getReferenceFrame2(const IndexType vtx) const {
    assert(vtx < m_numVertices - 1);

    return m_referenceFrames2[vtx];
  }

  void setReferenceFrame2(const IndexType vtx, const Vec3x& vec) {
    assert(vtx < m_numVertices - 1);

    m_referenceFrames2.set(vtx, vec);
  }

  const Vec3x getMaterialFrame1(const IndexType vtx) const {
    assert(vtx < m_numVertices - 1);

    return m_materialFrames1[vtx];
  }

  void setMaterialFrame1(const IndexType vtx, const Vec3x& vec) {
    assert(vtx < m_numVertices - 1);

    m_materialFrames1.set(vtx, vec);
  }

  const Vec3x getMaterialFrame2(const IndexType vtx) const {
    assert(vtx < m_numVertices - 1);

    return m_materialFrames2[vtx];
  }

  void setMaterialFrame2(const IndexType vtx, const Vec3x& vec) {
    assert(vtx < m_numVertices - 1);

    m_materialFrames2.set(vtx, vec);
  }

  Vec3x closestPoint(const Vec3x& x) const;

  bool hasSmallForces(const Scalar lTwoTol, const Scalar lInfTol) const;

  std::pair<Vec3x, Vec3x> getAABB() const;
  std::pair<Vec3x, Vec3x> getAABB(unsigned elementID) const;
  std::pair<Vec3x, Vec3x> getAABB(const Mat4x& transform) const;
  std::pair<Vec3x, Vec3x> getAABB(unsigned elementID,
                                  const Mat4x& transform) const;

  /**
   * Member variables
   */
 public:
  // Convenience copy of the original number of vertices (owned by the strand).
  const IndexType m_numVertices;

  ////////////////////////////////////////

  AreaDOFs m_area_dofs;
  DOFs m_dofs;
  Edges m_edges;
  JumpEdges m_jump_edges;
  NeighborDistances m_neighbor_distances;
  Lengths m_lengths;
  DynamicVoronoiLengths m_dynamic_voronoi_lengths;
  Height m_height;
  LaplaceHeight m_laplace_height;
  Tangents m_tangents;
  mutable ReferenceFrames1 m_referenceFrames1;
  mutable ReferenceFrames2 m_referenceFrames2;
  ReferenceTwists m_referenceTwists;
  Twists m_twists;
  CurvatureBinormals m_curvatureBinormals;
  TrigThetas m_trigThetas;
  mutable MaterialFrames<1> m_materialFrames1;
  mutable MaterialFrames<2> m_materialFrames2;
  Kappas m_kappas;
  GradKappas m_gradKappas;
  GradTwists m_gradTwists;
  GradTwistsSquared m_gradTwistsSquared;
  HessKappas m_hessKappas;
  HessTwists m_hessTwists;
  BendingProducts m_bendingProducts;
};

}  // namespace strandsim

#endif /* STRANDSTATE_HH_ */
