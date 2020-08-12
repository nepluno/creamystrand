/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ContinuousTimeCollision.hh"

#include "../Core/ElasticStrand.hh"
#include "CollisionUtils.hh"
#include "EdgeFaceCollision.hh"
#include "ElementProxy.hh"
#include "VertexFaceCollision.hh"

namespace strandsim {

static const double SQ_TOLERANCE = 1e-12;
static const double EXTRA_RADIUS = 1.e-3;

bool compareTimes(const CollisionBase* cb1, const CollisionBase* cb2) {
  const ContinuousTimeCollision* ct1 =
      dynamic_cast<const ContinuousTimeCollision*>(cb1);
  const ContinuousTimeCollision* ct2 =
      dynamic_cast<const ContinuousTimeCollision*>(cb2);

  if (ct1 == NULL || ct2 == NULL)
    return false;
  else
    return ct1->m_time < ct2->m_time;
}

bool compareCT(const CollisionBase* ct1, const CollisionBase* ct2) {
  {
    const VertexFaceCollision* vf1 =
        dynamic_cast<const VertexFaceCollision*>(ct1);
    const VertexFaceCollision* vf2 =
        dynamic_cast<const VertexFaceCollision*>(ct2);

    if (vf1 != NULL && vf2 == NULL) return true;

    if (vf1 == NULL && vf2 != NULL

    )
      return false;

    if (vf1 != NULL && vf2 != NULL) {
      return compare(vf1, vf2);
    }
  }

  {
    const EdgeFaceCollision* ef1 = dynamic_cast<const EdgeFaceCollision*>(ct1);
    const EdgeFaceCollision* ef2 = dynamic_cast<const EdgeFaceCollision*>(ct2);

    if (ef1 != NULL && ef2 == NULL) return true;

    if (ef1 == NULL && ef2 != NULL) return false;

    if (ef1 != NULL && ef2 != NULL) {
      return compare(ef1, ef2);
    }
  }

  std::cout << "We should never arrive here!\n";
  return false;
}

void ContinuousTimeCollision::postAnalyse(const Vec3x& relativeDisplacement) {
  m_normalRelativeDisplacement = relativeDisplacement.dot(m_normal);
  if (m_normalRelativeDisplacement > 0.0) {
    m_normal = -m_normal;
    m_normalRelativeDisplacement = -m_normalRelativeDisplacement;
  }
  m_tangentialRelativeDisplacement =
      relativeDisplacement - m_normalRelativeDisplacement * m_normal;

  m_normalRelativeDisplacement -= EXTRA_RADIUS;
}

Vec3x ContinuousTimeCollision::offset() const {
  return (m_offset.dot(m_normal) + EXTRA_RADIUS) * m_normal;
}

} /* namespace strandsim */
