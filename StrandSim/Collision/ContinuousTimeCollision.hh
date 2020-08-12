/**
 * \copyright 2012 Jean-Marie Aubry
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef CONTINUOUSTIMECOLLISION_HH_
#define CONTINUOUSTIMECOLLISION_HH_

#include "CollisionBase.hh"

namespace strandsim {

class ElasticStrand;
class FaceProxy;
class TriangularMesh;

bool compareTimes(const CollisionBase* ct1, const CollisionBase* ct2);

class ContinuousTimeCollision : public CollisionBase {
 public:
  ContinuousTimeCollision(ElasticStrand* firstStrand, int firstVertex)
      : m_firstStrand(firstStrand), m_firstVertex(firstVertex) {}
  virtual ~ContinuousTimeCollision() {}

  ElasticStrand* getFirstStrand() const { return m_firstStrand; }
  int getFirstVertex() const { return m_firstVertex; }

  virtual bool notIn(const std::set<const ElasticStrand*>& alreadySeenStrands) {
    return find(alreadySeenStrands.begin(), alreadySeenStrands.end(),
                m_firstStrand) == alreadySeenStrands.end();
  }

  virtual void putIn(std::set<const ElasticStrand*>& alreadySeenStrands) {
    alreadySeenStrands.insert(m_firstStrand);
  }

  friend bool compareTimes(const CollisionBase* ct1, const CollisionBase* ct2);

  Scalar time() const { return m_time; }
  const Vec3x& normal() const { return m_normal; }
  virtual Vec3x offset() const;

  // protected:
 public:
  void postAnalyse(const Vec3x& relativeDisplacement);

  ElasticStrand* m_firstStrand;
  int m_firstVertex;

  Scalar m_time;
  Vec3x m_normal;
  Vec3x m_offset;
  Scalar m_normalRelativeDisplacement;
  Vec3x m_tangentialRelativeDisplacement;
};

class EdgeCollision : public ContinuousTimeCollision {
 public:
  EdgeCollision(ElasticStrand* firstStrand, int firstVertex)
      : ContinuousTimeCollision(firstStrand, firstVertex) {}

  Scalar abscissa() const { return m_s; }

 public:  // FIXME please
  Scalar m_s;
};

} /* namespace strandsim */
#endif /* CONTINUOUSTIMECOLLISION_HH_ */
