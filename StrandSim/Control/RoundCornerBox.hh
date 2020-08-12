/**
 * \copyright 2018 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef RoundCornerBox_hpp
#define RoundCornerBox_hpp

#include <fstream>
#include <string>
#include <vector>

#include "../Collision/TriangularMesh.hh"

using namespace strandsim;

namespace strandsim {
class RoundCornerBox {
  strandsim::TriangularMesh* m_mesh;

 protected:
  std::vector<int> m_index_to_verts;
  double m_radius;
  int m_N_edge;
  inline void AddVertex(int i, int j, int k, const Vec3x& pos,
                        const Vec3x& base_pos) {
    int pidx = k * m_N_edge * m_N_edge + j * m_N_edge + i;
    if (m_index_to_verts[pidx] < 0) {
      int next_idx = (int)m_mesh->m_vertices.size();
      m_index_to_verts[pidx] = next_idx;

      Vec3x dir = pos - base_pos;
      if (dir.norm() > 0.0) {
        dir.normalize();
        m_mesh->addVertex(base_pos + dir * m_radius);
      } else {
        m_mesh->addVertex(pos);
      }
    }
  }
  inline int translateIndices(int i, int j, int k) {
    int pidx = k * m_N_edge * m_N_edge + j * m_N_edge + i;
    return m_index_to_verts[pidx];
  }
  inline void AddFace(int i, int j, int k, bool inversed) {
    if (inversed) {
      m_mesh->addFace(i, k, j);
    } else {
      m_mesh->addFace(i, j, k);
    }
  }

 public:
  RoundCornerBox(int N, const Vec3x& b, const double& radius,
                 strandsim::TriangularMesh* mesh, bool inverted);
};
};  // namespace strandsim

#endif /* CornerBox_hpp */
