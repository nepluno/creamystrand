/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "TriangularMesh.hh"

#include <igl/adjacency_list.h>
#include <igl/adjacency_matrix.h>
#include <igl/doublearea.h>
#include <igl/exact_geodesic.h>
#include <igl/grad.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/point_simplex_squared_distance.h>
#include <igl/principal_curvature.h>
#include <igl/vertex_triangle_adjacency.h>

#include <set>
#include <unordered_map>
#include <unordered_set>

#include "../Utils/MathUtilities.hh"

namespace strandsim {

TriangularMesh::TriangularMesh(
    const std::shared_ptr<MeshScriptingController>& controller)
    : m_associatedController(controller), m_associatedFlow(nullptr) {}

void TriangularMesh::getRigidTransform(int face_idx, Vec3x& center,
                                       Vec3x& translate,
                                       Eigen::Quaternion<Scalar>& rot) {
  const int p = m_faces[face_idx](0), q = m_faces[face_idx](1),
            r = m_faces[face_idx](2);

  Vec3x old_center =
      (m_stored_vertices[p] + m_stored_vertices[q] + m_stored_vertices[r]) /
      3.0;
  Vec3x new_center = (m_vertices[p] + m_vertices[q] + m_vertices[r]) / 3.0;

  center = old_center;

  translate = new_center - old_center;
  Vec3x old_n = ((m_stored_vertices[r] - m_stored_vertices[p])
                     .cross(m_stored_vertices[q] - m_stored_vertices[p]))
                    .normalized();
  Vec3x new_n =
      ((m_vertices[r] - m_vertices[p]).cross(m_vertices[q] - m_vertices[p]))
          .normalized();

  Vec3x a = old_n.cross(new_n);
  rot = Eigen::Quaternion<Scalar>(1.0 + old_n.dot(new_n), a(0), a(1), a(2));
  rot.normalize();
}

int TriangularMesh::getClosestTriangle(const Vec3x& pos) {
  int faceidx = -1;
  double mindist = 1e+20;

  Vec3i ipos =
      Vec3i(clamp((int)((pos(0) - m_indexer_origin(0)) / m_indexer_bucket_size),
                  0, m_indexer_ni - 1),
            clamp((int)((pos(1) - m_indexer_origin(1)) / m_indexer_bucket_size),
                  0, m_indexer_nj - 1),
            clamp((int)((pos(2) - m_indexer_origin(2)) / m_indexer_bucket_size),
                  0, m_indexer_nk - 1));

  for (int k = ipos(2) - 1; k < ipos(2) + 2; ++k)
    for (int j = ipos(1) - 1; j < ipos(1) + 2; ++j)
      for (int i = ipos(0) - 1; i < ipos(0) + 2; ++i) {
        if (k < 0 || k >= m_indexer_nk || j < 0 || j >= m_indexer_nj || i < 0 ||
            i >= m_indexer_ni)
          continue;

        const std::vector<int>& cell =
            m_triangle_indexer[k * m_indexer_nj * m_indexer_ni +
                               j * m_indexer_ni + i];
        for (int f : cell) {
          const int p = m_faces[f](0), q = m_faces[f](1), r = m_faces[f](2);

          Scalar t1, t2, t3;
          double dist = point_triangle_distance(
              pos, m_vertices[p], m_vertices[q], m_vertices[r], t1, t2, t3);
          if (dist < mindist) {
            mindist = dist;
            faceidx = f;
          }
        }
      }

  return faceidx;
}

void TriangularMesh::computeAdjacency(int num_rings) {
  int num_faces = nf();
  m_neighbor_local_faces.resize(num_faces);
  m_neighbors_face_global_indices.resize(num_faces);
  m_neighbors_vertex_global_indices.resize(num_faces);
  m_global_to_local_faces.resize(num_faces);

  // compute adj matrix
  MatXi F(num_faces, 3);
  for (int i = 0; i < num_faces; ++i) {
    for (int j = 0; j < 3; ++j) F(i, j) = m_faces[i](j);
  }

  Eigen::SparseMatrix<Scalar> mat_adj;
  igl::adjacency_matrix(F, mat_adj);

  igl::adjacency_list(F, m_vertex_adjacency);

  int num_verts = nv();

  // build vertex-triangle adjacency
  std::vector<std::vector<int> > vertex_triangles_incidency;

  igl::vertex_triangle_adjacency(num_verts, F, m_vertex_triangles,
                                 vertex_triangles_incidency);

  // for each face, for all 3 vertices, find adj vertices in # rings
  for_each(0, num_faces, [&](int i) {
    m_global_to_local_faces[i] = -1;

    VecXx test_vec(num_verts);
    test_vec.setZero();

    for (int j = 0; j < 3; ++j) test_vec(m_faces[i](j)) = 1;

    for (int k = 0; k < num_rings; ++k) test_vec = mat_adj * test_vec;

    std::vector<int> adj_verts;
    std::unordered_set<int> adj_verts_set;

    // for each face, check if adj faces are covered by adj vertices
    for (int j = 0; j < num_verts; ++j) {
      if (test_vec[j] > 0.) {
        adj_verts.push_back(j);
        adj_verts_set.insert(j);
      }
    }

    std::unordered_map<int, int> passed_adj_verts_set;

    m_neighbors_face_global_indices[i].clear();

    std::unordered_set<int> unique_neighbor_faces;

    for (int i_vert : adj_verts) {
      const std::vector<int>& vf = m_vertex_triangles[i_vert];
      for (int f : vf) {
        if (adj_verts_set.find(m_faces[f](0)) != adj_verts_set.end() &&
            adj_verts_set.find(m_faces[f](1)) != adj_verts_set.end() &&
            adj_verts_set.find(m_faces[f](2)) != adj_verts_set.end()) {
          // for each face, put passed adj faces into indices array
          unique_neighbor_faces.insert(f);
          for (int j = 0; j < 3; ++j) {
            passed_adj_verts_set[m_faces[f](j)] = -1;
          }
        }
      }
    }

    for (int f : unique_neighbor_faces) {
      m_neighbors_face_global_indices[i].push_back(f);
    }

    // for each face, compute local set of vertices from all passed adj faces
    std::vector<int>& passed_adj_verts = m_neighbors_vertex_global_indices[i];
    passed_adj_verts.resize(passed_adj_verts_set.size());

    int k = 0;
    for (auto& p : passed_adj_verts_set) {
      passed_adj_verts[k] = p.first;
      p.second = k++;
    }

    const int num_passed_faces = m_neighbors_face_global_indices[i].size();
    m_neighbor_local_faces[i].resize(num_passed_faces, 3);

    for (int j = 0; j < num_passed_faces; ++j) {
      int f = m_neighbors_face_global_indices[i][j];
      if (f == i) {
        m_global_to_local_faces[i] = j;
      }
      m_neighbor_local_faces[i].row(j) =
          Vec3i(passed_adj_verts_set[m_faces[f](0)],
                passed_adj_verts_set[m_faces[f](1)],
                passed_adj_verts_set[m_faces[f](2)])
              .transpose();
    }
  });

  m_reorder_operator.resize(num_faces * 3, num_faces * 3);
  m_reorder_operator.reserve(num_faces * 3);

  std::vector<Eigen::Triplet<Scalar> > tri_reorder;
  tri_reorder.reserve(num_faces * 3);
  for (int i = 0; i < num_faces; ++i) {
    for (int r = 0; r < 3; ++r) {
      tri_reorder.push_back(
          Eigen::Triplet<Scalar>(i * 3 + r, r * num_faces + i, 1.0));
    }
  }
  m_reorder_operator.setFromTriplets(tri_reorder.begin(), tri_reorder.end());
}

void TriangularMesh::updateFaceNormalArea() {
  const int num_verts = nv();
  const int num_faces = nf();

  MatXx V(num_verts, 3);
  MatXi F(num_faces, 3);

  for_each(0, num_verts, [&](int i) { V.row(i) = m_vertices[i].transpose(); });

  for_each(0, num_faces, [&](int i) {
    F.row(i) = Vec3i(m_faces[i](0), m_faces[i](1), m_faces[i](2)).transpose();
  });

  MatXx N;

  // normal
  igl::per_face_normals(V, F, N);

  m_face_frames.resize(num_faces * 3, 3);

  for_each(0, num_faces, [&](int i) {
    Vec3x n = N.row(i).transpose();
    Vec3x r0 = ((Mat3x::Identity() - n * n.transpose()) *
                (m_vertices[m_faces[i](1)] - m_vertices[m_faces[i](0)]))
                   .normalized();
    Vec3x r1 = r0.cross(n);

    m_face_frames.block<3, 1>(i * 3, 0) = r0;
    m_face_frames.block<3, 1>(i * 3, 1) = n;
    m_face_frames.block<3, 1>(i * 3, 2) = r1;
  });

  // areas
  igl::doublearea(V, F, m_face_areas);
  m_face_areas *= 0.5;

  m_vertex_areas.resize(num_verts);
  m_vertex_normals.resize(num_verts);

  MatXx VN;

  igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,
                          VN);

  for_each(0, num_verts, [&](int i) {
    const std::vector<int>& neighbor_faces = m_vertex_triangles[i];
    Scalar vert_area = 0.0;

    for (int f : neighbor_faces) {
      vert_area += m_face_areas[f];
    }

    m_vertex_areas[i] = vert_area / 3.0;
    m_vertex_normals[i] = VN.block<1, 3>(i, 0).transpose();
  });

  m_G_f.resize(num_faces * 3, num_faces * 3);
  m_G_v.resize(num_verts, num_verts);
  m_inv_G_f.resize(num_faces * 3, num_faces * 3);
  m_inv_G_v.resize(num_verts, num_verts);

  std::vector<Eigen::Triplet<Scalar> > tri_G_f;
  std::vector<Eigen::Triplet<Scalar> > tri_G_v;
  std::vector<Eigen::Triplet<Scalar> > tri_inv_G_f;
  std::vector<Eigen::Triplet<Scalar> > tri_inv_G_v;

  tri_G_f.reserve(num_faces * 3);
  tri_G_v.reserve(num_verts);
  tri_inv_G_f.reserve(num_faces * 3);
  tri_inv_G_v.reserve(num_verts);

  for (int i = 0; i < num_faces; ++i) {
    for (int r = 0; r < 3; ++r) {
      tri_G_f.push_back(
          Eigen::Triplet<Scalar>(i * 3 + r, i * 3 + r, m_face_areas[i]));
      tri_inv_G_f.push_back(
          Eigen::Triplet<Scalar>(i * 3 + r, i * 3 + r, 1.0 / m_face_areas[i]));
    }
  }

  for (int i = 0; i < num_verts; ++i) {
    tri_G_v.push_back(Eigen::Triplet<Scalar>(i, i, m_vertex_areas[i]));
    tri_inv_G_v.push_back(
        Eigen::Triplet<Scalar>(i, i, 1.0 / m_vertex_areas[i]));
  }

  m_G_f.setFromTriplets(tri_G_f.begin(), tri_G_f.end());
  m_G_v.setFromTriplets(tri_G_v.begin(), tri_G_v.end());
  m_inv_G_f.setFromTriplets(tri_inv_G_f.begin(), tri_inv_G_f.end());
  m_inv_G_v.setFromTriplets(tri_inv_G_v.begin(), tri_inv_G_v.end());

  // interpolants
  m_weight_face_vertex.resize(num_verts, num_faces);
  std::vector<Eigen::Triplet<Scalar> > tri_weights;
  tri_weights.reserve(num_verts * 6);
  for (int i = 0; i < num_verts; ++i) {
    const std::vector<int>& neighbor_faces = m_vertex_triangles[i];

    for (int f : neighbor_faces) {
      tri_weights.push_back(Eigen::Triplet<Scalar>(
          i, f, m_face_areas[f] / (m_vertex_areas[i] * 3.0)));
    }
  }
  m_weight_face_vertex.setFromTriplets(tri_weights.begin(), tri_weights.end());

  m_weight_vertex_face.resize(num_faces, num_verts);
  tri_weights.clear();
  tri_weights.reserve(num_faces * 3);
  for (int i = 0; i < num_faces; ++i) {
    const TriangularFace& f = m_faces[i];
    for (int r = 0; r < 3; ++r) {
      tri_weights.push_back(Eigen::Triplet<Scalar>(i, f(r), 1.0 / 3.0));
    }
  }
  m_weight_vertex_face.setFromTriplets(tri_weights.begin(), tri_weights.end());

  // gradient
  igl::grad(V, F, m_grad_operator);
  m_grad_operator = m_reorder_operator * m_grad_operator;

  // divergence
  m_div_operator = -m_inv_G_v * m_grad_operator.transpose() * m_G_f;

  m_laplacian_operator = -m_div_operator * m_grad_operator;
}

int TriangularMesh::getNeighborFaceAtDir(int i_vert, const Vec3x& dir) {
  Vec3x pos = m_vertices[i_vert] + dir;

  const std::vector<int>& neighbor_faces = m_vertex_triangles[i_vert];
  Scalar min_dist2 = 1e+40;
  int ret = -1;

  for (int f : neighbor_faces) {
    Vec3x center = (m_vertices[m_faces[f](0)] + m_vertices[m_faces[f](1)] +
                    m_vertices[m_faces[f](2)]) /
                   3.0;
    Scalar dist = (pos - center).squaredNorm();
    if (dist < min_dist2) {
      min_dist2 = dist;
      ret = f;
    }
  }

  return ret;
}

void TriangularMesh::mapFaceVecToVertex(const VecXx& face_vec,
                                        VecXx& vertex_vec, int nrow) {
  vertex_vec.resize(m_weight_face_vertex.rows() * nrow);
  const MatXx& face_mat =
      Eigen::Map<const MatXx>(face_vec.data(), nrow, face_vec.size() / nrow);
  Eigen::Map<MatXx>(vertex_vec.data(), nrow, vertex_vec.size() / nrow) =
      face_mat * m_weight_face_vertex.transpose();
}

void TriangularMesh::computeFlattenedNeighbors() {
  // for each face, grab the current vertices
  const int num_faces = nf();
  m_flattened_neighbors.resize(num_faces);

  for (int i = 0; i < num_faces; ++i) {
    const std::vector<int>& neigh_verts_indices =
        m_neighbors_vertex_global_indices[i];
    const int num_local_verts = neigh_verts_indices.size();

    MatXx vert_pos_local(num_local_verts, 3);

    Vec3x center = (m_vertices[m_faces[i](0)] + m_vertices[m_faces[i](1)] +
                    m_vertices[m_faces[i](2)]) /
                   3.0;

    for (int j = 0; j < num_local_verts; ++j) {
      vert_pos_local.row(j) =
          (m_vertices[neigh_verts_indices[j]] - center).transpose();
    }

    const int num_local_faces = m_neighbor_local_faces[i].rows();

    VecXi FS(1);
    FS(0) = m_global_to_local_faces[i];

    assert(FS(0) >= 0);

    VecXi vert_indices_target(num_local_verts);

    for (int j = 0; j < num_local_verts; ++j) {
      vert_indices_target(j) = j;
    }

    VecXx D;

    igl::exact_geodesic(vert_pos_local, m_neighbor_local_faces[i], VecXi(), FS,
                        vert_indices_target, VecXi(), D);

    // project to local space
    const Mat3x& R = m_face_frames.block<3, 3>(i * 3, 0);

    m_flattened_neighbors[i].resize(num_local_verts, 2);

    for (int j = 0; j < num_local_verts; ++j) {
      Vec3x lp = (vert_pos_local.row(j) * R).transpose();

      // exp projection
      m_flattened_neighbors[i].row(j) =
          Vec2x(lp(0), lp(2)).normalized().transpose() * D(j);
    }
  }
}

void TriangularMesh::getClosestBarycentric(int iFaceCenter, const Vec3x& offset,
                                           int& target, Vec3x& coord) {
  Vec3x local_offset =
      m_face_frames.block<3, 3>(iFaceCenter * 3, 0).transpose() * offset;
  MatXx local_offset_proj(1, 2);
  local_offset_proj.block<1, 2>(0, 0) =
      Vec2x(local_offset(0), local_offset(2)).transpose();

  VecXx sqrD;
  VecXi indices;
  MatXx C;

  igl::point_mesh_squared_distance(
      local_offset_proj, m_flattened_neighbors[iFaceCenter],
      m_neighbor_local_faces[iFaceCenter], sqrD, indices, C);

  int local_target = indices(0);
  Vec2x qp = C.row(0).transpose();
  Scalar sqrdist;
  Vec2x cp;

  igl::point_simplex_squared_distance<2>(qp, m_flattened_neighbors[iFaceCenter],
                                         m_neighbor_local_faces[iFaceCenter],
                                         local_target, sqrdist, cp, coord);

  target = m_neighbors_face_global_indices[iFaceCenter][local_target];
}

void TriangularMesh::buildTriangleIndexer(const Scalar& bucket_size) {
  m_indexer_bucket_size = bucket_size;

  computeBBox();

  m_indexer_origin = m_bbox_center - m_bbox_extent;

  m_indexer_ni = (int)ceil(m_bbox_extent(0) * 2.0 / bucket_size);
  m_indexer_nj = (int)ceil(m_bbox_extent(1) * 2.0 / bucket_size);
  m_indexer_nk = (int)ceil(m_bbox_extent(2) * 2.0 / bucket_size);

  m_triangle_indexer.resize(m_indexer_ni * m_indexer_nj * m_indexer_nk);

  auto triangle_bbx = [](const Vec3x& v0, const Vec3x& v1, const Vec3x& v2,
                         Vec3x& bbx_min, Vec3x& bbx_max) {
    bbx_min(0) = std::min(std::min(v0(0), v1(0)), v2(0));
    bbx_min(1) = std::min(std::min(v0(1), v1(1)), v2(1));
    bbx_min(2) = std::min(std::min(v0(2), v1(2)), v2(2));
    bbx_max(0) = std::max(std::max(v0(0), v1(0)), v2(0));
    bbx_max(1) = std::max(std::max(v0(1), v1(1)), v2(1));
    bbx_max(2) = std::max(std::max(v0(2), v1(2)), v2(2));
  };

  for (int i = 0; i < (int)m_triangle_indexer.size(); ++i) {
    m_triangle_indexer[i].clear();
  }

  const int num_faces = nf();
  for (int i = 0; i < num_faces; ++i) {
    const int p = m_faces[i](0), q = m_faces[i](1), r = m_faces[i](2);
    Vec3x bbx_min, bbx_max;

    triangle_bbx(m_vertices[p], m_vertices[q], m_vertices[r], bbx_min, bbx_max);

    Vec3i ibbx_min = Vec3i(
        clamp((int)((bbx_min(0) - m_indexer_origin(0)) / m_indexer_bucket_size),
              0, m_indexer_ni - 1),
        clamp((int)((bbx_min(1) - m_indexer_origin(1)) / m_indexer_bucket_size),
              0, m_indexer_nj - 1),
        clamp((int)((bbx_min(2) - m_indexer_origin(2)) / m_indexer_bucket_size),
              0, m_indexer_nk - 1));

    Vec3i ibbx_max = Vec3i(
        clamp((int)((bbx_max(0) - m_indexer_origin(0)) / m_indexer_bucket_size),
              0, m_indexer_ni - 1),
        clamp((int)((bbx_max(1) - m_indexer_origin(1)) / m_indexer_bucket_size),
              0, m_indexer_nj - 1),
        clamp((int)((bbx_max(2) - m_indexer_origin(2)) / m_indexer_bucket_size),
              0, m_indexer_nk - 1));

    for (int t = ibbx_min(2); t <= ibbx_max(2); ++t)
      for (int s = ibbx_min(1); s <= ibbx_max(1); ++s)
        for (int r = ibbx_min(0); r <= ibbx_max(0); ++r) {
          m_triangle_indexer[t * m_indexer_nj * m_indexer_ni +
                             s * m_indexer_ni + r]
              .push_back(i);
        }
  }
}

void TriangularMesh::tagBoundaryEdge(unsigned f, int firstVertex,
                                     int secondVertex) {
  TriangularFace& face = m_faces[f];

  for (short k = 0; k < 3; ++k) {
    const short next = (k + 1) % 3;
    const short prev = (k + 2) % 3;
    if (firstVertex == face.idx[k]) {
      if (secondVertex == face.idx[next]) {
        face.boundaryEdges |= (1 << k);
      } else if (secondVertex == face.idx[prev]) {
        face.boundaryEdges |= (1 << prev);
      }
    }
  }
}

void TriangularMesh::clear() {
  m_vertices.clear();
  m_displacements.clear();
  m_previous_displacements.clear();
  m_faces.clear();
}

TriangularMesh::~TriangularMesh() {
  // TODO Auto-generated destructor stub
}

} /* namespace strandsim */
