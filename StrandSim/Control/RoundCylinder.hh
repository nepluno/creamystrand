/**
 * \copyright 2018 Yun (Raymond) Fei
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ROUNDCYLINDER_H__
#define ROUNDCYLINDER_H__

/* Rounded Cylinder generation algorithm.
 */

#include "../Collision/TriangularMesh.hh"
#include <vector>
#include <unordered_map>

namespace strandsim {
class RoundCylinder
{
public:
	RoundCylinder(int N, int M, const double& ra, const double& rb, const double& h, strandsim::TriangularMesh* mesh, bool inverted);
};
};


#endif
