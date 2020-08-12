/**
 * \copyright 2018 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ICOSPHERE_H__
#define ICOSPHERE_H__

/* Icosphere generation algorithm.
 * Adapted from Andreas Kahler's C# implementation found here:
 * http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
 */

#include <unordered_map>
#include <vector>

#include "../Collision/TriangularMesh.hh"

using namespace strandsim;

namespace strandsim {
class Icosphere {
  int index;
  std::unordered_map<uint64, int> middlePointIndexCache;
  strandsim::TriangularMesh* m_mesh;

  int addVertexWithIndices(const Vec3x& p);
  int getMiddlePoint(int p1, int p2);

 public:
  Icosphere(int recursionLevel, const double& radius,
            strandsim::TriangularMesh* mesh, bool inverted);
};
};  // namespace strandsim

#endif
