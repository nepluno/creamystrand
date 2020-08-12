/**
 * \copyright 2018 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "RoundCylinder.hh"
#ifndef M_PI
#define M_PI 3.141592653589793238462650288
#endif
#ifndef M_PI_2
#define M_PI_2 (M_PI / 2.)
#endif
#include <cmath>

using namespace strandsim;

namespace strandsim {
RoundCylinder::RoundCylinder(int N, int M, const double& ra, const double& rb,
                             const double& h, strandsim::TriangularMesh* mesh,
                             bool inverted) {
  mesh->m_vertices.clear();
  mesh->m_faces.clear();

  const int num_verts = 2 + N * (M + 1) * 2;
  mesh->m_vertices.resize(num_verts);
  mesh->m_displacements.resize(num_verts);

  // top and bottom centers
  mesh->m_vertices[0] = Vec3x(0, h + rb, 0);
  mesh->m_vertices[1] = Vec3x(0, -h - rb, 0);

  const Scalar dtheta = M_PI * 2.0 / (Scalar)N;
  const Scalar dphi = M_PI_2 / (Scalar)M;

  for (int i = 0; i < N; ++i) {
    double theta = dtheta * (Scalar)i;
    double ct = cos(theta);
    double st = sin(theta);

    for (int j = 0; j <= M; ++j) {
      double phi = dphi * (Scalar)j;
      double cp = cos(phi);
      double sp = sin(phi);

      const int index_u = 2 + i * (M + 1) + j;
      const int index_d = index_u + N * (M + 1);

      mesh->m_vertices[index_u](0) = mesh->m_vertices[index_d](0) =
          (ra + rb * cp) * ct;
      mesh->m_vertices[index_u](1) = h + rb * sp;
      mesh->m_vertices[index_d](1) = -(h + rb * sp);
      mesh->m_vertices[index_u](2) = mesh->m_vertices[index_d](2) =
          (ra + rb * cp) * st;
    }
  }

  mesh->m_faces.resize(N * (M + 1) * 4);

  int base_idx = 0;
  // top and bottom cap
  for (int i = 0; i < N; ++i) {
    int j = (i + 1) % N;
    mesh->m_faces[i] =
        TriangularFace(0, 2 + j * (M + 1) + M, 2 + i * (M + 1) + M);
    mesh->m_faces[i + N] =
        TriangularFace(1, 2 + (i + N) * (M + 1) + M, 2 + (j + N) * (M + 1) + M);
  }

  base_idx += 2 * N;

  // side pads
  for (int i = 0; i < N; ++i) {
    int j = (i + 1) % N;
    mesh->m_faces[base_idx + i * 2 + 0] =
        TriangularFace(2 + j * (M + 1), 2 + (i + N) * (M + 1), 2 + i * (M + 1));
    mesh->m_faces[base_idx + i * 2 + 1] = TriangularFace(
        2 + j * (M + 1), 2 + (j + N) * (M + 1), 2 + (i + N) * (M + 1));
  }

  base_idx += 2 * N;

  // upper bevel
  for (int i = 0; i < N; ++i) {
    int j = (i + 1) % N;
    for (int k = 0; k < M; ++k) {
      mesh->m_faces[base_idx + (i * M + k) * 2 + 0] =
          TriangularFace(2 + j * (M + 1) + (k + 1), 2 + i * (M + 1) + k,
                         2 + i * (M + 1) + (k + 1));
      mesh->m_faces[base_idx + (i * M + k) * 2 + 1] = TriangularFace(
          2 + j * (M + 1) + (k + 1), 2 + j * (M + 1) + k, 2 + i * (M + 1) + k);
    }
  }

  base_idx += M * N * 2;

  // lower bevel
  for (int i = 0; i < N; ++i) {
    int j = (i + 1) % N;
    for (int k = 0; k < M; ++k) {
      mesh->m_faces[base_idx + (i * M + k) * 2 + 0] = TriangularFace(
          2 + (j + N) * (M + 1) + k, 2 + (i + N) * (M + 1) + (k + 1),
          2 + (i + N) * (M + 1) + k);
      mesh->m_faces[base_idx + (i * M + k) * 2 + 1] = TriangularFace(
          2 + (j + N) * (M + 1) + k, 2 + (j + N) * (M + 1) + (k + 1),
          2 + (i + N) * (M + 1) + (k + 1));
    }
  }

  if (inverted) {
    for (TriangularFace& f : mesh->m_faces) {
      f.invert();
    }
  }
}
};  // namespace strandsim
