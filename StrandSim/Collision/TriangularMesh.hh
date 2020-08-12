/**
 * \copyright 2012 Jean-Marie Aubry, 2019 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef TRIANGULARMESH_HH_
#define TRIANGULARMESH_HH_

#include <Eigen/Sparse>

#include "../Core/Definitions.hh"
#include "../Utils/ThreadUtils.hh"

namespace strandsim {

class MeshScriptingController;
class TriangularMeshFlow;

struct TriangularFace {
  int idx[3];
  int boundaryEdges;

  TriangularFace() {
    idx[0] = idx[1] = idx[2] = 0;
    boundaryEdges = 0;
  }

  TriangularFace(const TriangularFace& f) {
    idx[0] = f.idx[0];
    idx[1] = f.idx[1];
    idx[2] = f.idx[2];
    boundaryEdges = f.boundaryEdges;
  }

  TriangularFace(int un, int deux, int trois) : boundaryEdges(0) {
    idx[0] = un;
    idx[1] = deux;
    idx[2] = trois;
  }

  bool edgeOnBoundary(const short edge) const {
    return (boundaryEdges & (1 << edge));
  }

  bool hasBoundaryEdges() const { return boundaryEdges; }

  int operator()(int i) const { return idx[i]; }

  void invert() { std::swap(idx[1], idx[2]); }
};

typedef std::vector<TriangularFace> FaceArray;

class TriangularMesh {
 public:
  TriangularMesh(
      const std::shared_ptr<MeshScriptingController>& controller = nullptr);
  virtual ~TriangularMesh();

  void buildTriangleIndexer(const Scalar& bucket_size);

  void getRigidTransform(int face_idx, Vec3x& center, Vec3x& translate,
                         Eigen::Quaternion<Scalar>& rot);

  int getClosestTriangle(const Vec3x& pos);

  void resizeVertices(size_t size) {
    m_vertices.resize(size);
    m_displacements.resize(size);
    m_previous_displacements.resize(size);
  }

  void updateFaceNormalArea();

  void resizeFaces(size_t size) { m_faces.resize(size); }

  void addVertex(const Vec3x& point) {
    m_vertices.push_back(point);
    m_displacements.push_back(Vec3x());
    m_previous_displacements.push_back(Vec3x());
  }

  Vec3x getFaceNormal(int fidx) {
    return m_face_frames.block<3, 1>(fidx * 3, 1);
  }

  void scale(const Vec3x& s) {
    for (Vec3x& v : m_vertices) {
      v.array() *= s.array();
    }
  }

  MutexType& getGeometryMutex() { return m_geometry_mutex; }

  Vec3x getVertex(size_t vtx) const { return m_vertices[vtx]; }

  void setVertex(size_t vtx, const Vec3x& point) { m_vertices[vtx] = point; }

  Vec3x getDisplacement(size_t vtx) const { return m_displacements[vtx]; }

  Vec3x getAcceleration(size_t vtx, const Scalar& dt) {
    return (m_displacements[vtx] - m_previous_displacements[vtx]) / (dt * dt);
  }

  void setDisplacement(size_t vtx, const Vec3x& disp) {
    m_displacements[vtx] = disp;
  }

  void saveDisplacement() { m_previous_displacements = m_displacements; }

  const TriangularFace& getFace(unsigned f) const { return m_faces[f]; }

  void addFace(int un, int deux, int trois) {
    m_faces.push_back(TriangularFace(un, deux, trois));
  }

  void tagBoundaryEdge(unsigned f, int firstVertex, int secondVertex);

  const std::vector<TriangularFace>& getFaces() const { return m_faces; }

  std::shared_ptr<MeshScriptingController> associatedController() const {
    return m_associatedController;
  }

  void setAssociatedController(
      const std::shared_ptr<MeshScriptingController>& controller) {
    m_associatedController = controller;
  }

  size_t nv() const { return m_vertices.size(); }

  size_t nf() const { return m_faces.size(); }

  Scalar getVertexArea(int vidx) const { return m_vertex_areas(vidx); }

  const Vec3xArray& vertices() const { return m_vertices; }

  Vec3xArray& vertices() { return m_vertices; }

  Vec3x getBBoxCenter() const { return m_bbox_center; }

  Vec3x getBBoxExtent() const { return m_bbox_extent; }

  void translate(const Vec3x& t) {
    for (Vec3x& v : m_vertices) {
      v += t;
    }
  }

  void transform(const Eigen::Quaternion<double>& transformation,
                 const Vec3x& center, const Vec3x& translate,
                 const Vec3x& scaling) {
    for (unsigned i = 0; i < nv(); ++i) {
      Vec3x vert = getVertex(i);
      Vec3x vertNext =
          transformation * Vec3x((vert - center).array() * scaling.array()) +
          center + translate;
      setVertex(i, vertNext);
      setDisplacement(i, (vertNext - m_stored_vertices[i]));
    }
  }

  void zeroDisplacement() {
    m_displacements.assign(m_displacements.size(), Vec3x::Zero());
  }

  void getClosestBarycentric(int iFaceCenter, const Vec3x& offset, int& target,
                             Vec3x& coord);

  void mapFaceVecToVertex(const VecXx& face_vec, VecXx& vertex_vec, int nrow);

  int getNeighborFaceAtDir(int i_vert, const Vec3x& dir);

  const std::vector<int>& getNeighborFacesAtVertex(int i_vert) const {
    return m_vertex_triangles[i_vert];
  }

  Vec3x getCenter() const {
    Vec3x ret = Vec3x::Zero();
    for (const Vec3x& v : m_vertices) {
      ret += v;
    }

    ret /= (Scalar)m_vertices.size();
    return ret;
  }

  void clear();

  void invert() {
    const int num_faces = nf();
    for (int i = 0; i < num_faces; ++i) {
      m_faces[i].invert();
    }
  }

  void computeBBox() {
    Vec3x bbmin = Vec3x(std::numeric_limits<Scalar>::infinity(),
                        std::numeric_limits<Scalar>::infinity(),
                        std::numeric_limits<Scalar>::infinity());

    Vec3x bbmax = Vec3x(-std::numeric_limits<Scalar>::infinity(),
                        -std::numeric_limits<Scalar>::infinity(),
                        -std::numeric_limits<Scalar>::infinity());

    int N = m_vertices.size();
    for (int i = 0; i < N; ++i) {
      const Vec3x& v = m_vertices[i];
      bbmin[0] = std::min(bbmin[0], v[0]);
      bbmin[1] = std::min(bbmin[1], v[1]);
      bbmin[2] = std::min(bbmin[2], v[2]);
      bbmax[0] = std::max(bbmax[0], v[0]);
      bbmax[1] = std::max(bbmax[1], v[1]);
      bbmax[2] = std::max(bbmax[2], v[2]);
    }

    m_bbox_center = (bbmax + bbmin) * 0.5;
    m_bbox_extent = (bbmax - bbmin) * 0.5;
  }

  void computeFlattenedNeighbors();
  void computeAdjacency(int num_rings);

  // Used for surface flow
  std::vector<MatXx> m_flattened_neighbors;

  std::vector<std::vector<int> > m_vertex_triangles;
  std::vector<std::vector<int> > m_neighbors_vertex_global_indices;
  std::vector<std::vector<int> > m_neighbors_face_global_indices;
  std::vector<MatXi> m_neighbor_local_faces;
  std::vector<int> m_global_to_local_faces;

  // Used for dynamics
  Vec3xArray m_stored_vertices;
  Vec3xArray m_vertices;
  Vec3xArray m_displacements;
  Vec3xArray m_previous_displacements;
  Vec3xArray m_vertex_normals;

  MatXx m_face_frames;
  std::vector<TriangularFace> m_faces;

  // DDG operators
  VecXx m_face_areas;
  VecXx m_vertex_areas;
  Eigen::SparseMatrix<Scalar> m_G_f;
  Eigen::SparseMatrix<Scalar> m_G_v;
  Eigen::SparseMatrix<Scalar> m_inv_G_f;
  Eigen::SparseMatrix<Scalar> m_inv_G_v;
  Eigen::SparseMatrix<Scalar> m_weight_face_vertex;
  Eigen::SparseMatrix<Scalar> m_weight_vertex_face;
  Eigen::SparseMatrix<Scalar> m_grad_operator;
  Eigen::SparseMatrix<Scalar> m_div_operator;
  Eigen::SparseMatrix<Scalar> m_laplacian_operator;
  Eigen::SparseMatrix<Scalar> m_reorder_operator;

  std::shared_ptr<MeshScriptingController> m_associatedController;
  std::shared_ptr<TriangularMeshFlow> m_associatedFlow;

  std::vector<std::vector<int> > m_vertex_adjacency;

  Vec3x m_bbox_center;
  Vec3x m_bbox_extent;

  std::vector<std::vector<int> > m_triangle_indexer;
  Vec3x m_indexer_origin;
  int m_indexer_ni;
  int m_indexer_nj;
  int m_indexer_nk;
  Scalar m_indexer_bucket_size;

  MutexType m_geometry_mutex;
};

} /* namespace strandsim */
#endif /* TRIANGULARMESH_HH_ */
