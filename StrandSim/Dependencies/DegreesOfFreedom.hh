/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef DEGREESOFFREEDOM_HH_
#define DEGREESOFFREEDOM_HH_

#include "DependencyNode.hh"

namespace strandsim {
class ElasticStrandParameters;
/**
 * Unit: cm^2 for area dofs
 */
class AreaDOFs : public DependencyNode<VecXx> {
 public:
  AreaDOFs(const VecXx& dofValues) : DependencyNode<VecXx>(dofValues) {
    m_numVertices = dofValues.size();
    m_numEdges = m_numVertices - 1;
    setClean();
  }

  // As this class is meant to be the root of the dependency tree, we provide a
  // const getter This means that DOFs are never "computed" but updated by the
  // outer loop algorithm
  const VecXx& get() const { return m_value; }

  Scalar getVertex(IndexType vtx) const {
    assert(vtx < m_numVertices);

    return get()[vtx];
  }

  void setVertex(IndexType vtx, const Scalar& point) {
    m_value[vtx] = point;
    setDependentsDirty();
  }

  IndexType getNumEdges() const { return m_numEdges; }

  virtual const char* name() const { return "AreaDOFs"; }

 protected:
  virtual void compute()  // Not implemented as this is an pure input node
  {
    ErrorStream(g_log, "")
        << "DegreesOfFreedom::compute() should never be called";
  }

  IndexType m_numVertices, m_numEdges;
};

/**
 * Unit: cm for position dofs, no dimension for theta
 */
class DOFs : public DependencyNode<VecXx> {
 public:
  DOFs(const VecXx& dofValues) : DependencyNode<VecXx>(dofValues) {
    assert(dofValues.size() % 4 == 3);
    m_numEdges = dofValues.size() / 4;
    m_numVertices = m_numEdges + 1;
    setClean();
  }

  // As this class is meant to be the root of the dependency tree, we provide a
  // const getter This means that DOFs are never "computed" but updated by the
  // outer loop algorithm
  const VecXx& get() const { return m_value; }

  Vec3x getVertex(IndexType vtx) const {
    assert(vtx < m_numVertices);

    return get().segment<3>(4 * vtx);
  }

  void setVertex(IndexType vtx, const Vec3x& point) {
    m_value.segment<3>(4 * vtx) = point;
    setDependentsDirty();
  }

  // Accessors to the theta degrees of freedom
  const Eigen::Map<const VecXx, Eigen::Unaligned, Eigen::InnerStride<4> >
  getThetas() const {
    return Eigen::Map<const VecXx, Eigen::Unaligned, Eigen::InnerStride<4> >(
        m_value.data() + 3, m_numEdges);
  }

  Scalar getTheta(IndexType vtx) const {
    assert(vtx < m_numEdges);

    return get()[4 * vtx + 3];
  }

  void setThetas(const VecXx& thetas, int numberOfFixedThetas = 0) {
    assert(thetas.size() == m_numEdges);

    Eigen::Map<VecXx, Eigen::Unaligned, Eigen::InnerStride<4> >(
        m_value.data() + 4 * numberOfFixedThetas + 3,
        m_numEdges - numberOfFixedThetas) =
        thetas.tail(m_numEdges - numberOfFixedThetas);
    setDependentsDirty();
  }

  void setTheta(IndexType vtx, Scalar theta) {
    m_value[4 * vtx + 3] = theta;
    setDependentsDirty();
  }

  IndexType getNumEdges() const { return m_numEdges; }

  IndexType getNumVertices() const { return m_numVertices; }

  virtual const char* name() const { return "DOFs"; }

 protected:
  virtual void compute()  // Not implemented as this is an pure input node
  {
    ErrorStream(g_log, "")
        << "DegreesOfFreedom::compute() should never be called";
  }

 private:
  IndexType m_numVertices, m_numEdges;
};

/**
 * Unit: cm
 */
class Edges : public DependencyNode<Vec3xArray> {
 public:
  Edges(DOFs& dofs)
      : DependencyNode<Vec3xArray>(0, dofs.getNumEdges()), m_dofs(dofs) {
    m_dofs.addDependent(this);
  }

  virtual const char* name() const { return "Edges"; }

 protected:
  virtual void compute();

  DOFs& m_dofs;
};

/**
 * Unit: cm
 */
class JumpEdges : public DependencyNode<Vec3xArray> {
 public:
  JumpEdges(DOFs& dofs)
      : DependencyNode<Vec3xArray>(0, dofs.getNumVertices()), m_dofs(dofs) {
    m_dofs.addDependent(this);
  }

  virtual const char* name() const { return "Jump Edges"; }

 protected:
  virtual void compute();

  DOFs& m_dofs;
};

/**
 * Unit: cm
 */
class Lengths : public DependencyNode<std::vector<Scalar> > {
 public:
  Lengths(Edges& edges)
      : DependencyNode<std::vector<Scalar> >(0, edges.size()), m_edges(edges) {
    m_edges.addDependent(this);
  }

  virtual const char* name() const { return "Lengths"; }

 protected:
  virtual void compute();

  Edges& m_edges;
};

/**
 * Unit: cm
 */
class NeighborDistances : public DependencyNode<std::vector<Scalar> > {
 public:
  NeighborDistances(JumpEdges& jump_edges)
      : DependencyNode<std::vector<Scalar> >(0, jump_edges.size()),
        m_jump_edges(jump_edges) {
    m_jump_edges.addDependent(this);
  }

  virtual const char* name() const { return "Neighbor Distance"; }

 protected:
  virtual void compute();

  JumpEdges& m_jump_edges;
};

/**
 * Unit: cm
 */
class DynamicVoronoiLengths : public DependencyNode<std::vector<Scalar> > {
 public:
  DynamicVoronoiLengths(Lengths& lengths)
      : DependencyNode<std::vector<Scalar> >(0, lengths.size() + 1),
        m_lengths(lengths) {
    m_lengths.addDependent(this);
  }

  virtual const char* name() const { return "Dynamic Voronoi Lengths"; }

 protected:
  virtual void compute();

  Lengths& m_lengths;
};

/**
 * Unit: unitless
 */
class Height : public DependencyNode<std::vector<Scalar> > {
 public:
  Height(AreaDOFs& dofs, const ElasticStrandParameters& parameters);

  virtual const char* name() const { return "Height"; }

 protected:
  virtual void compute();

  AreaDOFs& m_area_dofs;
  const ElasticStrandParameters& m_parameters;
};

/**
 * Unit: cm^{-1}
 */
class LaplaceHeight : public DependencyNode<std::vector<Scalar> > {
 public:
  LaplaceHeight(Height& height, DynamicVoronoiLengths& lengths)
      : DependencyNode<std::vector<Scalar> >(0, lengths.size()),
        m_voronoi_lengths(lengths),
        m_heights(height) {
    m_voronoi_lengths.addDependent(this);
    m_heights.addDependent(this);
  }

  virtual const char* name() const { return "Laplace Height"; }

 protected:
  virtual void compute();

  Height& m_heights;
  DynamicVoronoiLengths& m_voronoi_lengths;
};

/**
 * Unit: no dimension
 */
class Tangents : public DependencyNode<Vec3xArray> {
 public:
  Tangents(Edges& edges, Lengths& lengths)
      : DependencyNode<Vec3xArray>(0, edges.size()),
        m_edges(edges),
        m_lengths(lengths) {
    m_edges.addDependent(this);
    m_lengths.addDependent(this);
  }

  virtual const char* name() const { return "Tangents"; }

 protected:
  virtual void compute();

  Edges& m_edges;
  Lengths& m_lengths;
};

/**
 * Unit: no dimension
 */
class CurvatureBinormals : public DependencyNode<Vec3xArray> {
 public:
  CurvatureBinormals(Tangents& tangents)
      : DependencyNode<Vec3xArray>(1, tangents.size()), m_tangents(tangents) {
    m_tangents.addDependent(this);
  }

  virtual const char* name() const { return "CurvatureBinormals"; }

 protected:
  virtual void compute();

  Tangents& m_tangents;
};

/**
 * Unit: no dimension
 */
class TrigThetas : public DependencyNode<std::pair<VecXx, VecXx> > {
 public:
  TrigThetas(DOFs& dofs)
      : DependencyNode<std::pair<VecXx, VecXx> >(
            std::make_pair(VecXx(), VecXx())),
        m_dofs(dofs) {
    m_dofs.addDependent(this);
  }

  virtual const char* name() const { return "TrigThetas"; }

  const VecXx& getSines() { return get().first; }

  const VecXx& getCosines() { return get().second; }

 protected:
  virtual void compute();

  DOFs& m_dofs;
};

}  // namespace strandsim

#endif /* DEGREESOFFREEDOM_HH_ */
