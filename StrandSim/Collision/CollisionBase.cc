/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "CollisionBase.hh"

#include "ElementProxy.hh"

namespace strandsim {

CollisionBase::CollisionBase() {
  // TODO Auto-generated constructor stub
}

CollisionBase::~CollisionBase() {
  // TODO Auto-generated destructor stub
}

void CollisionBase::print(std::ostream& os) const {
  os << "Print function not implemented yet";
}

Scalar FaceCollision::faceFrictionCoefficient() const {
  const FaceProxy* faceProxy = face();

  if (faceProxy) {
    auto controller = faceProxy->getMesh()->associatedController();
    if (controller) {
      return faceProxy->getFrictionCoefficient(m_u, m_v, m_w);
    }
  }

  return 0.;
}

std::ostream& operator<<(std::ostream& os, const CollisionBase& collision) {
  collision.print(os);

  return os;
}

bool sameCT(const CollisionBase* c1, const CollisionBase* c2) {
  return !compareCT(c1, c2) && !compareCT(c2, c1);
}

} /* namespace strandsim */
