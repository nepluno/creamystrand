/**
 * \copyright 2018 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef CAPSULE_H__
#define CAPSULE_H__

/* Capsule generation algorithm.
 * Adapted from Paul Bourke's C implementation found here:
 * http://paulbourke.net/geometry/capsule/
 */

#include <unordered_map>
#include <vector>

#include "../Collision/TriangularMesh.hh"

namespace strandsim {
class Capsule {
 public:
  Capsule(int N, const double& radius, const double& halfheight,
          strandsim::TriangularMesh* mesh, bool inverted);
};
};  // namespace strandsim

#endif
