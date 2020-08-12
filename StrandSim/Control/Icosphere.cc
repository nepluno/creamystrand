/**
 * \copyright 2018 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Icosphere.hh"

using namespace strandsim;

namespace strandsim {
Icosphere::Icosphere(int recursionLevel, const double& radius,
                     strandsim::TriangularMesh* mesh, bool inverted)
    : index(0), m_mesh(mesh) {
  middlePointIndexCache.clear();
  m_mesh->m_vertices.clear();
  m_mesh->m_displacements.clear();
  m_mesh->m_faces.clear();

  const double t = (1.0 + sqrt(5.0)) / 2.0;

  addVertexWithIndices(Vec3x(-1, t, 0));
  addVertexWithIndices(Vec3x(1, t, 0));
  addVertexWithIndices(Vec3x(-1, -t, 0));
  addVertexWithIndices(Vec3x(1, -t, 0));

  addVertexWithIndices(Vec3x(0, -1, t));
  addVertexWithIndices(Vec3x(0, 1, t));
  addVertexWithIndices(Vec3x(0, -1, -t));
  addVertexWithIndices(Vec3x(0, 1, -t));

  addVertexWithIndices(Vec3x(t, 0, -1));
  addVertexWithIndices(Vec3x(t, 0, 1));
  addVertexWithIndices(Vec3x(-t, 0, -1));
  addVertexWithIndices(Vec3x(-t, 0, 1));

  // 5 faces around point 0
  m_mesh->addFace(0, 11, 5);
  m_mesh->addFace(0, 5, 1);
  m_mesh->addFace(0, 1, 7);
  m_mesh->addFace(0, 7, 10);
  m_mesh->addFace(0, 10, 11);

  m_mesh->addFace(1, 5, 9);
  m_mesh->addFace(5, 11, 4);
  m_mesh->addFace(11, 10, 2);
  m_mesh->addFace(10, 7, 6);
  m_mesh->addFace(7, 1, 8);

  m_mesh->addFace(3, 9, 4);
  m_mesh->addFace(3, 4, 2);
  m_mesh->addFace(3, 2, 6);
  m_mesh->addFace(3, 6, 8);
  m_mesh->addFace(3, 8, 9);

  m_mesh->addFace(4, 9, 5);
  m_mesh->addFace(2, 4, 11);
  m_mesh->addFace(6, 2, 10);
  m_mesh->addFace(8, 6, 7);
  m_mesh->addFace(9, 8, 1);

  std::vector<TriangularFace>& indices = m_mesh->m_faces;
  // refine triangles
  for (int i = 0; i < recursionLevel; i++) {
    std::vector<TriangularFace> indices2;
    for (const TriangularFace& tri : indices) {
      // replace triangle by 4 triangles
      int a = getMiddlePoint(tri(0), tri(1));
      int b = getMiddlePoint(tri(1), tri(2));
      int c = getMiddlePoint(tri(2), tri(0));

      indices2.push_back(TriangularFace(tri(0), a, c));
      indices2.push_back(TriangularFace(tri(1), b, a));
      indices2.push_back(TriangularFace(tri(2), c, b));
      indices2.push_back(TriangularFace(a, b, c));
    }
    m_mesh->m_faces = indices2;
  }

  for (auto& v : m_mesh->m_vertices) {
    v *= radius;
  }

  if (inverted) {
    for (TriangularFace& f : mesh->m_faces) {
      f.invert();
    }
  }
}

int Icosphere::addVertexWithIndices(const Vec3x& p) {
  m_mesh->m_vertices.push_back(p.normalized());
  m_mesh->m_displacements.push_back(Vec3x());
  return index++;
}

int Icosphere::getMiddlePoint(int p1, int p2) {
  int smallerIndex = std::min(p1, p2);
  int greaterIndex = std::max(p1, p2);

  uint64 key = ((uint64)(smallerIndex) << 32UL) | (uint64)greaterIndex;
  auto itr = middlePointIndexCache.find(key);
  if (itr != middlePointIndexCache.end()) {
    return itr->second;
  }

  Vec3x middle = (m_mesh->m_vertices[p1] + m_mesh->m_vertices[p2]) * 0.5;
  int i = addVertexWithIndices(middle);

  middlePointIndexCache[key] = i;
  return i;
}

}  // namespace strandsim
