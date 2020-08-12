/**
 * \copyright 2013 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ORIENTEDBOUNDINGBOX_HH_
#define ORIENTEDBOUNDINGBOX_HH_

#include "../Core/CollisionParameters.hh"

namespace strandsim {

class StrandState;
class ElasticStrand;

class OrientedBoundingBox {
 public:
  OrientedBoundingBox() {}
  OrientedBoundingBox(
      const StrandState& geometry, const unsigned edge,
      const CollisionParameters& params,
      CollisionParameters::CollisionType type = CollisionParameters::SELF);
  OrientedBoundingBox(
      const ElasticStrand& strand, const unsigned edge,
      CollisionParameters::CollisionType type = CollisionParameters::SELF);

  Vec3x m_center;
  Vec3x m_axis[3];
  Scalar m_extents[3];

  void getAABB(Vec3x& min, Vec3x& max);

  bool intersects(const OrientedBoundingBox& other) const {
    return areIntersecting(*this, other);
  }

  static bool areIntersecting(const OrientedBoundingBox& obb1,
                              const OrientedBoundingBox& obb2);

 private:
  void init(const StrandState& geometry, const unsigned edge,
            const CollisionParameters& params,
            CollisionParameters::CollisionType type);
};

} /* namespace strandsim */
#endif /* ORIENTEDBOUNDINGBOX_HH_ */
