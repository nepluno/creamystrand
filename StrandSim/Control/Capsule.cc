/**
 * \copyright 2018 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Capsule.hh"
#ifndef M_PI
#define M_PI 3.141592653589793238462650288
#endif
#include <cmath>

using namespace strandsim;

namespace strandsim {
Capsule::Capsule(int N, const double& radius, const double& halfheight,
                 strandsim::TriangularMesh* mesh, bool inverted) {
  mesh->m_vertices.clear();
  mesh->m_faces.clear();

  const int num_verts = (N + 1) * (N / 2 + 2);
  mesh->m_vertices.resize(num_verts);
  mesh->m_displacements.resize(num_verts);

  int N_4 = N / 4;
  int N_2 = N / 2;

  int index = 0;
  for (int j = 0; j <= N_4; ++j) {
    for (int i = 0; i <= N; ++i) {
      double theta = (double)i * M_PI * 2.0 / (double)N;
      double phi = -M_PI / 2.0 + M_PI * (double)j / (double)N_2;
      mesh->m_vertices[index](0) = radius * sin(phi) - halfheight;
      mesh->m_vertices[index](1) = radius * cos(phi) * sin(theta);
      mesh->m_vertices[index](2) = radius * cos(phi) * cos(theta);
      ++index;
    }
  }

  for (int j = N_4; j <= N_2; ++j) {
    for (int i = 0; i <= N; ++i) {
      double theta = (double)i * M_PI * 2.0 / (double)N;
      double phi = -M_PI / 2.0 + M_PI * (double)j / (double)N_2;
      mesh->m_vertices[index](0) = radius * sin(phi) + halfheight;
      mesh->m_vertices[index](1) = radius * cos(phi) * sin(theta);
      mesh->m_vertices[index](2) = radius * cos(phi) * cos(theta);
      ++index;
    }
  }

  for (int i = 0; i < num_verts; ++i) {
    mesh->m_displacements[i].setZero();
  }

  mesh->m_faces.resize((N_2 + 1) * N * 2);
  index = 0;
  for (int j = 0; j <= N_2; ++j) {
    for (int i = 0; i < N; ++i) {
      int i1 = j * (N + 1) + i;
      int i2 = j * (N + 1) + (i + 1);
      int i3 = (j + 1) * (N + 1) + (i + 1);
      int i4 = (j + 1) * (N + 1) + i;
      mesh->m_faces[index++] = TriangularFace(i1, i3, i2);
      mesh->m_faces[index++] = TriangularFace(i1, i4, i3);
    }
  }

  if (inverted) {
    for (TriangularFace& f : mesh->m_faces) {
      f.invert();
    }
  }
}
};  // namespace strandsim
